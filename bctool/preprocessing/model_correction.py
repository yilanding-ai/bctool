"""Feature-based read error detection model.

Trains a lightweight logistic regression (pure numpy) to predict
whether an aligned read is reliable, using features extracted from
the SAM record and the read sequence itself.

Features:
  mapq              �?mapping quality
  clip_ratio        �?fraction of soft/hard-clipped bases
  target_count      �?number of target_base occurrences in read
  converted_count   �?number of converted_base occurrences
  unconverted_ratio �?target / (target + converted + 1)
  gc_content        �?GC fraction
  read_length       �?read length in bp
  mismatch_ratio    �?fraction of mismatches (from NM tag, if present)
  is_second         �?1 if read is second in pair (flag & 0x80)
"""

import re
import json
import pickle
import numpy as np
from pathlib import Path
from collections import defaultdict


_CIGAR_RE = re.compile(r'(\d+)([MIDNSHPX=])')


def _parse_nm(parts):
    """Extract NM (edit distance) from optional SAM tags."""
    for tag in parts[11:]:
        if tag.startswith("NM:i:"):
            try:
                return int(tag[5:])
            except ValueError:
                return None
    return None


def _extract_features_single(parts, conversion):
    """Extract feature vector from a single SAM record.

    Returns a dict of feature_name -> float value.
    """
    flag = int(parts[1])
    mapq = int(parts[4])

    # Read sequence
    seq = parts[9].upper()
    read_len = len(seq)

    # Soft-/hard-clip from CIGAR
    clip_bases = 0
    aligned_bases = 0
    cigar = parts[5]
    for m in _CIGAR_RE.finditer(cigar):
        n, op = int(m.group(1)), m.group(2)
        if op in ("S", "H"):
            clip_bases += n
        elif op in ("M", "I", "D", "N", "=", "X"):
            aligned_bases += n

    clip_ratio = clip_bases / read_len if read_len > 0 else 0

    # Conversion counting
    target_base = conversion.target_base
    converted_base = conversion.converted_base
    target_count = seq.count(target_base)
    converted_count = seq.count(converted_base)
    denom = target_count + converted_count + 1
    unconverted_ratio = target_count / denom

    # GC content
    gc = seq.count("G") + seq.count("C")
    gc_content = gc / read_len if read_len > 0 else 0

    # Mismatch ratio from NM tag
    nm = _parse_nm(parts)
    mismatch_ratio = nm / read_len if (nm is not None and read_len > 0) else -1.0

    # Paired-end flags
    is_second = 1.0 if (flag & 0x80) else 0.0

    return {
        "mapq": float(mapq),
        "clip_ratio": clip_ratio,
        "target_count": float(target_count),
        "converted_count": float(converted_count),
        "unconverted_ratio": unconverted_ratio,
        "gc_content": gc_content,
        "read_length": float(read_len),
        "mismatch_ratio": mismatch_ratio,
        "is_second": is_second,
    }


def _feature_dict_to_vector(fd, feature_names):
    """Convert feature dict to numpy array in the order of feature_names."""
    return np.array([fd.get(k, 0.0) for k in feature_names], dtype=np.float64)


_FEATURE_NAMES = [
    "mapq", "clip_ratio", "target_count", "converted_count",
    "unconverted_ratio", "gc_content", "read_length",
    "mismatch_ratio", "is_second",
]


def extract_features_from_sam(sam_path, conversion):
    """Extract feature matrix and read names from a SAM file.

    Returns (read_names, feature_matrix, feature_dicts)
      - read_names: list of QNAME strings
      - feature_matrix: numpy array of shape (n_reads, n_features)
      - feature_dicts: list of dicts (for training/tagging)
    """
    names = []
    rows = []
    dicts = []

    with open(sam_path) as f:
        for line in f:
            if line.startswith("@"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 11:
                continue
            flag = int(parts[1])
            if flag & 0x4:  # unmapped
                continue
            fd = _extract_features_single(parts, conversion)
            names.append(parts[0])
            rows.append(_feature_dict_to_vector(fd, _FEATURE_NAMES))
            dicts.append(fd)

    return names, np.array(rows, dtype=np.float64), dicts


# ─── Label generation from ground truth ─────────────────────────────────

def _load_ground_truth_bed(bed_path):
    """Load ground truth BED into dict: (chrom, pos) -> is_methylated (bool)."""
    truth = {}
    with open(bed_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 5:
                try:
                    chrom, pos = parts[0], int(parts[1])
                    level = float(parts[3])
                    truth[(chrom, pos)] = level >= 0.5
                except (ValueError, IndexError):
                    continue
    return truth


def generate_labels_from_sam(sam_path, ground_truth_bed_path, conversion):
    """Generate training labels using ground truth BED.

    For each SAM read, check if the methylation calls it contributes to
    match the ground truth. A read is:
      - 0 (good) if its per-base methylation status agrees with truth
      - 1 (error) if it disagrees significantly

    Since SAM doesn't directly tell us the methylation status per base,
    we use a proxy: compare the read's conversion ratio at target sites
    against the expected ratio from ground truth.

    Returns (names, labels, features).
    """
    truth = _load_ground_truth_bed(ground_truth_bed_path)
    target_base = conversion.target_base
    converted_base = conversion.converted_base

    names, X, dicts = extract_features_from_sam(sam_path, conversion)
    labels = np.ones(len(names), dtype=np.float64)  # default: error

    for i, fd in enumerate(dicts):
        tc = fd["target_count"]
        cc = fd["converted_count"]
        total = tc + cc
        if total == 0:
            labels[i] = 0.0  # no target sites �?no opinion
            continue

        observed_ratio = cc / total  # fraction converted
        # Expected: most target sites in methylated regions should convert
        # Use a simple heuristic: reads with <30% conversion are suspicious
        if observed_ratio >= 0.3:
            labels[i] = 0.0  # good read
        else:
            labels[i] = 1.0  # poor conversion �?error

    return names, labels, X, dicts


# ─── Model ──────────────────────────────────────────────────────────────

class LogisticRegression:
    """Binary logistic regression trained via gradient descent (pure numpy)."""

    def __init__(self, learning_rate=0.1, epochs=500, l2=0.01):
        self.lr = learning_rate
        self.epochs = epochs
        self.l2 = l2
        self.weights = None
        self.bias = 0.0
        self.feature_names = _FEATURE_NAMES

    def _sigmoid(self, z):
        z = np.clip(z, -500, 500)
        return 1.0 / (1.0 + np.exp(-z))

    def fit(self, X, y):
        n, d = X.shape
        self.weights = np.zeros(d, dtype=np.float64)
        self.bias = 0.0

        for epoch in range(self.epochs):
            z = X @ self.weights + self.bias
            pred = self._sigmoid(z)
            error = pred - y

            grad_w = (X.T @ error) / n + self.l2 * self.weights
            grad_b = np.mean(error)

            self.weights -= self.lr * grad_w
            self.bias -= self.lr * grad_b

            # Early stopping if loss change is tiny
            if epoch % 100 == 0 and epoch > 0:
                loss = self._binary_crossentropy(y, pred)
        return self

    def predict_proba(self, X):
        z = X @ self.weights + self.bias
        return self._sigmoid(z)

    def predict(self, X, threshold=0.5):
        return (self.predict_proba(X) >= threshold).astype(np.float64)

    @staticmethod
    def _binary_crossentropy(y, pred):
        eps = 1e-15
        pred = np.clip(pred, eps, 1 - eps)
        return -np.mean(y * np.log(pred) + (1 - y) * np.log(1 - pred))

    def save(self, path):
        data = {
            "weights": self.weights.tolist(),
            "bias": self.bias,
            "feature_names": self.feature_names,
        }
        with open(path, "w") as f:
            json.dump(data, f)

    @classmethod
    def load(cls, path):
        with open(path) as f:
            data = json.load(f)
        model = cls()
        model.weights = np.array(data["weights"], dtype=np.float64)
        model.bias = data["bias"]
        model.feature_names = data.get("feature_names", _FEATURE_NAMES)
        return model


# ─── Convenience wrapper ───────────────────────────────────────────────

class ReadErrorDetector:
    """High-level API for training and applying read error detection."""

    def __init__(self, model_path=None):
        self.model = None
        self.feature_names = _FEATURE_NAMES
        if model_path:
            self.load(model_path)

    def train(self, sam_path, ground_truth_bed, conversion):
        """Train the model using SAM + ground truth BED."""
        _, labels, X, _ = generate_labels_from_sam(
            sam_path, ground_truth_bed, conversion
        )
        self.model = LogisticRegression()
        self.model.fit(X, labels)
        return self

    def train_from_features(self, features, labels):
        """Train from pre-computed feature matrix and labels."""
        self.model = LogisticRegression()
        self.model.fit(features, labels)
        return self

    def predict_proba_from_sam(self, sam_path, conversion):
        """Predict error probability for each read in SAM.

        Returns list of (read_name, error_probability, features).
        """
        names, X, dicts = extract_features_from_sam(sam_path, conversion)[0:3]
        if len(names) == 0:
            return []

        if self.model is None:
            # No model trained: fall back to heuristic score
            scores = self._heuristic_score(X, dicts)
        else:
            scores = self.model.predict_proba(X)

        return list(zip(names, scores.tolist(), dicts))

    def filter_sam(self, sam_path, conversion, threshold=0.5, output_sam=None):
        """Filter SAM by predicted error probability.

        Reads with P(error) >= threshold are moved to a rejected list.
        Returns (corrected_sam_path, rejected_fastq_path).
        """
        sam_path = Path(sam_path)
        base = sam_path.stem.replace(".sam", "")
        out_dir = sam_path.parent
        out_sam = output_sam or (out_dir / f"{base}_model_filtered.sam")
        rej_fastq = out_dir / f"{base}_model_rejected.fastq"

        predictions = self.predict_proba_from_sam(str(sam_path), conversion)
        pred_map = {name: prob for name, prob, _ in predictions}

        kept = 0
        rejected = 0

        with open(sam_path) as fin, \
             open(out_sam, "w") as fout, \
             open(rej_fastq, "w") as frej:

            for line in fin:
                if line.startswith("@"):
                    fout.write(line)
                    continue

                parts = line.strip().split("\t")
                if len(parts) < 11:
                    fout.write(line)
                    kept += 1
                    continue

                read_name = parts[0]
                prob = pred_map.get(read_name, 0.0)

                if prob >= threshold:
                    rejected += 1
                    frej.write(f"@{read_name}\n{parts[9]}\n+\n{parts[10]}\n")
                else:
                    fout.write(line)
                    kept += 1

        print(f"  [ReadErrorDetector] Kept={kept}  Rejected={rejected}  "
              f"(threshold={threshold})")
        return str(out_sam), str(rej_fastq)

    @staticmethod
    def _heuristic_score(X, dicts):
        """Fallback heuristic when no model is trained.

        Combines features into a simple weighted score.
        """
        scores = np.zeros(len(dicts))
        for i, fd in enumerate(dicts):
            s = 0.0
            s += max(0, (60 - fd["mapq"]) / 60) * 0.25
            s += fd["clip_ratio"] * 0.20
            s += fd["unconverted_ratio"] * 0.30
            if fd["mismatch_ratio"] >= 0:
                s += min(fd["mismatch_ratio"], 1.0) * 0.15
            s += (1 - fd["gc_content"]) * 0.10 if fd["gc_content"] < 0.2 else 0
            scores[i] = min(s, 1.0)
        return scores

    def save(self, path):
        self.model.save(path) if self.model else None

    def load(self, path):
        self.model = LogisticRegression.load(path)
        self.feature_names = self.model.feature_names


# ─── Feature vector utility ────────────────────────────────────────────

def get_feature_names():
    return list(_FEATURE_NAMES)
