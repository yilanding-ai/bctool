import numpy as np
from pathlib import Path
from collections import defaultdict
from .base import AlignerBase
from ..methylation.converter import ConversionType


_ENCODE_2BIT = {"A": 0b00, "C": 0b01, "G": 0b10, "T": 0b11, "N": 0b00}
_DECODE_2BIT = {0b00: "A", 0b01: "C", 0b10: "G", 0b11: "T"}
_BASE2INT = {"A": 0, "C": 1, "G": 2, "T": 3}

_CONV3_MAP_WIDE = {
    "C>T": {"A": 0, "C": 1, "G": 2, "T": 1},
    "T>C": {"A": 0, "C": 1, "G": 2, "T": 1},
    "A>G": {"A": 1, "C": 0, "G": 1, "T": 2},
    "G>A": {"A": 1, "C": 0, "G": 1, "T": 2},
    "A>C": {"A": 1, "C": 1, "G": 0, "T": 2},
    "C>A": {"A": 1, "C": 1, "G": 0, "T": 2},
    "G>T": {"A": 2, "C": 0, "G": 1, "T": 1},
    "T>G": {"A": 2, "C": 0, "G": 1, "T": 1},
    "A>T": {"A": 1, "C": 0, "G": 2, "T": 1},
    "T>A": {"A": 1, "C": 0, "G": 2, "T": 1},
    "C>G": {"A": 0, "C": 1, "G": 1, "T": 2},
    "G>C": {"A": 0, "C": 1, "G": 1, "T": 2},
}


def encode_2bit(seq):
    n = len(seq)
    words = (n + 31) // 32
    packed = np.zeros(words, dtype=np.uint64)
    for i, base in enumerate(seq):
        code = _ENCODE_2BIT.get(base.upper(), 0b00)
        packed[i // 32] |= np.uint64(code) << (2 * (i % 32))
    return packed


def decode_2bit(packed, length):
    chars = []
    for i in range(length):
        code = (packed[i // 32] >> (2 * (i % 32))) & 0b11
        chars.append(_DECODE_2BIT.get(code, "N"))
    return "".join(chars)


def to_3base(seq, conv_key):
    mapping = _CONV3_MAP_WIDE[conv_key]
    return np.array([mapping.get(b.upper(), 0) for b in seq], dtype=np.uint8)


def build_score_lut(conv_key):
    conv = ConversionType(conv_key)
    target = conv.target_base
    converted = conv.converted_base
    lut = np.full((4, 4), -1, dtype=np.int8)
    for i, ref_base in enumerate("ACGT"):
        for j, read_base in enumerate("ACGT"):
            if ref_base == read_base:
                lut[i, j] = 1
            elif (ref_base == target and read_base == converted) or \
                 (ref_base == converted and read_base == target):
                lut[i, j] = 0
    lut[0, 0] = 1
    return lut


def count_mismatches_2bit(ref_packed, read_packed, conv_key, length):
    """Count mismatches using 2-bit XOR with conversion-aware read re-encoding."""
    conv = ConversionType(conv_key)
    target = conv.target_base
    converted = conv.converted_base
    enc_target = _ENCODE_2BIT[target]
    enc_converted = _ENCODE_2BIT[converted]
    re_encode_map = {enc_target: enc_target, enc_converted: enc_target}
    mask_map = {}
    for base, code in _ENCODE_2BIT.items():
        if base in ("A", "C", "G", "T"):
            mask_map[code] = code
    mask_map[enc_converted] = enc_target
    n_words = (length + 31) // 32
    mismatches = 0
    for w in range(n_words):
        rw = int(ref_packed[w]) if w < len(ref_packed) else 0
        rd = int(read_packed[w]) if w < len(read_packed) else 0
        for pos in range(32):
            global_pos = w * 32 + pos
            if global_pos >= length:
                break
            ref_code = (rw >> (2 * pos)) & 0b11
            read_code = (rd >> (2 * pos)) & 0b11
            masked_read = mask_map.get(read_code, read_code)
            if ref_code != masked_read:
                mismatches += 1
    return mismatches


def _build_score_lut_fast(conv_key):
    lut = np.full((4, 4), -1, dtype=np.int8)
    for i, ref_base in enumerate("ACGT"):
        for j, read_base in enumerate("ACGT"):
            if ref_base == read_base:
                lut[i, j] = 1
    conv = ConversionType(conv_key)
    ti = "ACGT".index(conv.target_base)
    ci = "ACGT".index(conv.converted_base)
    lut[ti, ci] = 0
    lut[ci, ti] = 0
    lut[0, 0] = 1
    return lut


class ConversionAwareAligner(AlignerBase):
    """Built-in aligner combining 3-base seed index with 2-bit bitmask technology.

    Supports all 12 single-base conversion types. Uses hash-based seed lookup
    on 3-base converted reference, then refines with conversion-aware scoring
    via 2-bit encoding and banded Smith-Waterman.
    """

    name = "conversion_aware"
    binary_name = "__conv_aware__"
    requires_index = False

    _instances = 0

    def __init__(self, config, output_dir, name_suffix="", seed_length=16,
                 min_seed_matches=2, match_score=1, conversion_score=0,
                 mismatch_penalty=-1, gap_open=-4, gap_extend=-1,
                 band_width_ratio=2.0, max_candidates=10, auto_adapt=True):
        super().__init__(config, output_dir)
        ConversionAwareAligner._instances += 1
        if name_suffix:
            self.name = f"conversion_aware_{name_suffix}"
        else:
            self.name = f"conversion_aware_{ConversionAwareAligner._instances}"
        self.aligner_dir = self.output_dir / self.name
        self.aligner_dir.mkdir(parents=True, exist_ok=True)
        self._seed_length = seed_length
        self._auto_adapt = auto_adapt
        self.min_seed_matches = min_seed_matches
        self.match_score = match_score
        self.conversion_score = conversion_score
        self.mismatch_penalty = mismatch_penalty
        self.gap_open = gap_open
        self.gap_extend = gap_extend
        self.band_width_ratio = band_width_ratio
        self.max_candidates = max_candidates
        self._genome = None
        self._genome_len = 0
        self._ref_name = "chr_sim"
        self._conv_idx = {}
        self._ref_2bit = None

    @staticmethod
    def _optimal_seed_length(read_length):
        if read_length < 40:
            return 8
        elif read_length < 60:
            return 10
        elif read_length < 80:
            return 12
        elif read_length < 120:
            return 14
        else:
            return 16

    def is_available(self):
        return True

    @property
    def seed_length(self):
        return self._seed_length

    @seed_length.setter
    def seed_length(self, value):
        self._seed_length = value

    def _read_fasta(self, fa_path):
        seqs = {}
        current_name = None
        current_seq = []
        with open(fa_path) as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_name and current_seq:
                        seqs[current_name] = "".join(current_seq)
                    current_name = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line.upper())
        if current_name and current_seq:
            seqs[current_name] = "".join(current_seq)
        return seqs

    @staticmethod
    def _detect_read_length(read1):
        try:
            with open(read1) as f:
                for _ in range(2):
                    line = f.readline()
                return len(f.readline().strip())
        except Exception:
            return 100

    def build_index(self, reference_fa):
        ref_path = Path(reference_fa)
        if not ref_path.exists():
            print(f"  [{self.name}] Warning: {reference_fa} not found, skipping index")
            return
        genomes = self._read_fasta(reference_fa)
        if not genomes:
            return
        first_name = list(genomes.keys())[0]
        self._ref_name = first_name
        genome_seq = genomes[first_name]
        self._genome = genome_seq
        self._genome_len = len(genome_seq)
        self._ref_2bit = encode_2bit(genome_seq)
        conv_key = ConversionType(self.config.conversion).key
        self._build_index_for_conv(conv_key)

    def _build_index_for_conv(self, conv_key):
        if conv_key in self._conv_idx:
            return
        seq_3base = to_3base(self._genome, conv_key)
        k = self._seed_length
        pow3_k = 3 ** k
        hash_table = defaultdict(list)
        current_hash = 0
        for i in range(self._genome_len):
            val = int(seq_3base[i])
            if i < k:
                current_hash = current_hash * 3 + val
                if i == k - 1:
                    hash_table[current_hash].append(0)
            else:
                old_val = int(seq_3base[i - k])
                current_hash = current_hash * 3 + val - old_val * pow3_k
                hash_table[current_hash].append(i - k + 1)
        self._conv_idx[conv_key] = {
            "seq_3base": seq_3base,
            "hash_table": dict(hash_table),
        }
        print(f"  [{self.name}] Built 3-base index for {conv_key} "
              f"(k={k}, ref_len={self._genome_len}, seeds={len(hash_table)})")

    def run_align(self, read1, read2=None, reference_fa=None, index_dir=None, threads=8):
        if self._genome is None and reference_fa:
            if self._auto_adapt:
                read_len = self._detect_read_length(read1)
                self._seed_length = self._optimal_seed_length(read_len)
            self.build_index(reference_fa)
        if self._genome is None:
            print(f"  [{self.name}] No genome loaded, creating empty SAM")
            out_sam = str(self.aligner_dir / "aligned.sam")
            with open(out_sam, "w") as f:
                f.write("@HD\tVN:1.6\tSO:coordinate\n")
                f.write(f"@SQ\tSN:{self._ref_name}\tLN:1000000\n")
            return out_sam
        conv_key = ConversionType(self.config.conversion).key
        self._build_index_for_conv(conv_key)
        out_sam = str(self.aligner_dir / "aligned.sam")
        self._generate_alignment(read1, read2, out_sam, conv_key)
        return out_sam

    def _generate_alignment(self, read1, read2, out_sam, conv_key):
        idx_data = self._conv_idx[conv_key]
        hash_table = idx_data["hash_table"]
        ref_3base = idx_data["seq_3base"]
        read_length = None
        seqs_1, quals_1, names_1 = self._read_fastq(read1)
        if read_length is None and seqs_1:
            read_length = len(seqs_1[0])
        is_paired = read2 is not None
        if is_paired:
            seqs_2, quals_2, names_2 = self._read_fastq(read2)
        band_width = max(50, int(read_length * self.band_width_ratio))
        gap_open = self.gap_open
        gap_extend = self.gap_extend
        k = self._seed_length
        # Read-adaptive tuning
        if read_length < 100:
            seed_step = max(1, read_length // 10)
            local_max_candidates = min(self.max_candidates, 5)
        else:
            seed_step = max(1, k // 2)
            local_max_candidates = self.max_candidates
        with open(out_sam, "w") as f:
            f.write(f"@HD\tVN:1.6\tSO:coordinate\n")
            f.write(f"@SQ\tSN:{self._ref_name}\tLN:{self._genome_len}\n")
            f.write(f"@PG\tID:{self.name}\tPN:conversion_aware\tVN:0.1.0\n")
            n_reads = len(seqs_1)
            for ri in range(n_reads):
                if ri % 1000 == 0 and ri > 0:
                    print(f"  [{self.name}] Aligned {ri}/{n_reads} reads...", end="\r")
                read_seq = seqs_1[ri]
                read_qual = quals_1[ri]
                read_name = names_1[ri]
                if is_paired:
                    read_seq2 = seqs_2[ri] if ri < len(seqs_2) else None
                    read_qual2 = quals_2[ri] if ri < len(quals_2) else None
                else:
                    read_seq2 = None
                result = self._align_single(read_seq, conv_key, hash_table,
                                            ref_3base, k, band_width, gap_open, gap_extend,
                                            seed_step, local_max_candidates)
                if result:
                    ref_start, cigar_str, score, nm = result
                    flag = 0 if not is_paired else 99
                    tlen = 0
                    f.write(f"{read_name}\t{flag}\t{self._ref_name}\t{ref_start + 1}\t"
                            f"60\t{cigar_str}\t*\t0\t{tlen}\t{read_seq}\t{read_qual}\t"
                            f"NM:i:{nm}\tAS:i:{score}\n")
                    if is_paired and read_seq2:
                        r2_result = self._align_single(read_seq2, conv_key,
                                                       hash_table, ref_3base,
                                                       k, band_width, gap_open, gap_extend,
                                                       seed_step, local_max_candidates)
                        if r2_result:
                            r2_start, r2_cigar, r2_score, r2_nm = r2_result
                            r2_flag = 147
                            f.write(f"{read_name}\t{r2_flag}\t{self._ref_name}\t"
                                    f"{r2_start + 1}\t60\t{r2_cigar}\t*\t0\t0\t"
                                    f"{read_seq2}\t{read_qual2}\t"
                                    f"NM:i:{r2_nm}\tAS:i:{r2_score}\n")
                else:
                    flag = 4 if not is_paired else 77
                    f.write(f"{read_name}\t{flag}\t*\t0\t0\t*\t*\t0\t0\t"
                            f"{read_seq}\t{read_qual}\n")
            print(f"  [{self.name}] Aligned {n_reads}/{n_reads} reads.")

    def _align_single(self, read_seq, conv_key, hash_table, ref_3base,
                      k, band_width, gap_open, gap_extend,
                      seed_step=None, max_candidates=None):
        L = len(read_seq)
        if L < k:
            return None
        read_3base = to_3base(read_seq, conv_key)
        read_2bit = encode_2bit(read_seq)
        candidates = defaultdict(list)
        if seed_step is None:
            seed_step = max(1, k // 2)
        if max_candidates is None:
            max_candidates = self.max_candidates
        for seed_start in range(0, L - k + 1, seed_step):
            seed = read_3base[seed_start:seed_start + k]
            shash = 0
            for val in seed:
                shash = shash * 3 + int(val)
            if shash in hash_table:
                for rpos in hash_table[shash]:
                    offset = rpos - seed_start
                    candidates[offset].append((seed_start, rpos, shash))
        best_score = -1e9
        best_pos = -1
        best_cigar = f"{L}M"
        best_nm = L
        sorted_candidates = sorted(candidates.items(),
                                   key=lambda x: len(x[1]), reverse=True)
        for offset, hits in sorted_candidates[:max_candidates]:
            n_seeds = len(set(h[0] for h in hits))
            if n_seeds < self.min_seed_matches:
                continue
            ref_start_candidate = offset
            if ref_start_candidate < 0:
                ref_start_candidate = 0
            if ref_start_candidate + L > self._genome_len:
                ref_start_candidate = self._genome_len - L
            if ref_start_candidate < 0:
                continue
            nm = count_mismatches_2bit(
                self._ref_2bit[ref_start_candidate // 32:],
                read_2bit, conv_key, L
            )
            ref_region = self._genome[ref_start_candidate:
                                       ref_start_candidate + L + band_width]
            score = self._sw_score(read_seq, ref_region, conv_key, gap_open, gap_extend)
            if score > best_score:
                best_score = score
                best_pos = ref_start_candidate
                best_cigar = f"{L}M"
                best_nm = nm
        if best_pos < 0:
            return None
        return best_pos, best_cigar, best_score, best_nm

    def _sw_score(self, read_seq, ref_region, conv_key, gap_open, gap_extend):
        lut = _build_score_lut_fast(conv_key)
        L = len(read_seq)
        R = len(ref_region)
        M = np.zeros((L + 1, R + 1), dtype=np.int32)
        Ix = np.full((L + 1, R + 1), -10**9, dtype=np.int32)
        Iy = np.full((L + 1, R + 1), -10**9, dtype=np.int32)
        ri_arr = np.array([_BASE2INT.get(b.upper(), 0) for b in read_seq], dtype=np.int32)
        rj_arr = np.array([_BASE2INT.get(b.upper(), 0) for b in ref_region], dtype=np.int32)
        max_score = 0
        for i in range(1, L + 1):
            ri = int(ri_arr[i - 1])
            s_col = lut[rj_arr, ri]
            diag = np.maximum(M[i - 1, :-1], np.maximum(Ix[i - 1, :-1], Iy[i - 1, :-1]))
            M[i, 1:] = np.maximum(0, diag + s_col)
            Ix[i, 1:] = np.maximum(M[i, :-1] + gap_open, Ix[i, :-1] + gap_extend)
            Iy[i, 1:] = np.maximum(M[i - 1, 1:] + gap_open, Iy[i - 1, 1:] + gap_extend)
            max_score = max(max_score, int(np.max(M[i])))
        return max_score

    def call_methylation(self, sam_bam_path, reference_fa=None):
        conv = ConversionType(self.config.conversion)
        target = conv.target_base
        converted = conv.converted_base
        out_path = str(self.aligner_dir / "methylation.bed")
        if self._genome is None and reference_fa:
            self.build_index(reference_fa)
        if self._genome is None:
            print(f"  [{self.name}] No genome, generating empty BED")
            Path(out_path).write_text("")
            return out_path
        ref_positions = {}
        for pos, base in enumerate(self._genome):
            if base == target:
                ref_positions[pos] = base
        with open(sam_bam_path) as sam:
            lines = [l for l in sam if not l.startswith("@")]
        chrom = self._ref_name
        site_counts = defaultdict(lambda: {"meth": 0, "unmeth": 0, "depth": 0})
        for line in lines:
            parts = line.strip().split("\t")
            if len(parts) < 11:
                continue
            flag = int(parts[1])
            if flag & 0x4:
                continue
            read_seq = parts[9]
            pos_str = parts[3]
            if pos_str == "*":
                continue
            ref_start = int(pos_str) - 1
            cigar_str = parts[5]
            ref_offset = ref_start
            read_offset = 0
            cigar_tuples = self._parse_cigar(cigar_str)
            for op, length in cigar_tuples:
                if op == "M":
                    for ci in range(length):
                        ri_pos = ref_offset + ci
                        rl_pos = read_offset + ci
                        if ri_pos < 0 or rl_pos >= len(read_seq):
                            continue
                        if ri_pos in ref_positions:
                            read_base = read_seq[rl_pos].upper()
                            site_counts[ri_pos]["depth"] += 1
                            if read_base == target:
                                site_counts[ri_pos]["meth"] += 1
                            elif read_base == converted:
                                site_counts[ri_pos]["unmeth"] += 1
                    ref_offset += length
                    read_offset += length
                elif op == "D":
                    ref_offset += length
                elif op == "I":
                    read_offset += length
                elif op == "S":
                    read_offset += length
                elif op in ("H", "N"):
                    ref_offset += length if op == "N" else 0
        with open(out_path, "w") as f:
            for pos in sorted(site_counts.keys()):
                sc = site_counts[pos]
                depth = sc["depth"]
                meth = sc["meth"]
                unmeth = sc["unmeth"]
                if depth == 0:
                    continue
                level = meth / depth
                f.write(f"{chrom}\t{pos}\t{pos+1}\t{level:.4f}\t"
                        f"{meth}\t{unmeth}\t{depth}\n")
        print(f"  [{self.name}] Called methylation: {len(site_counts)} sites")
        return out_path

    @staticmethod
    def _parse_cigar(cigar_str):
        if not cigar_str or cigar_str == "*":
            return []
        result = []
        i = 0
        num = ""
        while i < len(cigar_str):
            if cigar_str[i].isdigit():
                num += cigar_str[i]
            else:
                if num:
                    result.append((cigar_str[i], int(num)))
                    num = ""
            i += 1
        return result

    @staticmethod
    def _read_fastq(path):
        fastq_path = Path(path)
        if not fastq_path.exists():
            return [], [], []
        seqs, quals, names = [], [], []
        with open(fastq_path) as f:
            lines = f.readlines()
        for i in range(0, len(lines), 4):
            if i + 3 >= len(lines):
                break
            names.append(lines[i].strip().lstrip("@"))
            seqs.append(lines[i + 1].strip())
            quals.append(lines[i + 3].strip())
        return seqs, quals, names
