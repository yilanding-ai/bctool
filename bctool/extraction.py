import random
import gzip
from pathlib import Path


def _open_fastq(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def _read_fastq_records(f):
    while True:
        header = f.readline()
        if not header:
            break
        seq = f.readline().rstrip()
        plus = f.readline()
        qual = f.readline().rstrip()
        yield header.rstrip(), seq, qual


def count_reads(fastq_path):
    n = 0
    with _open_fastq(fastq_path) as f:
        for _ in _read_fastq_records(f):
            n += 1
    return n


def extract_reads(read1, read2, output_dir, num_reads, seed=42):
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    rng = random.Random(seed)

    if read2:
        # paired-end reservoir sampling
        with _open_fastq(read1) as f1, _open_fastq(read2) as f2:
            reservoir1, reservoir2 = [], []
            for i, (r1, r2) in enumerate(zip(_read_fastq_records(f1),
                                              _read_fastq_records(f2))):
                if i < num_reads:
                    reservoir1.append(r1)
                    reservoir2.append(r2)
                else:
                    j = rng.randint(0, i)
                    if j < num_reads:
                        reservoir1[j] = r1
                        reservoir2[j] = r2
            _write_fastq(out / "extracted_R1.fastq", reservoir1, gzip=False)
            _write_fastq(out / "extracted_R2.fastq", reservoir2, gzip=False)
    else:
        with _open_fastq(read1) as f:
            reservoir = []
            for i, rec in enumerate(_read_fastq_records(f)):
                if i < num_reads:
                    reservoir.append(rec)
                else:
                    j = rng.randint(0, i)
                    if j < num_reads:
                        reservoir[j] = rec
            _write_fastq(out / "extracted.fastq", reservoir, gzip=False)

    print(f"Extracted {num_reads} reads to {out}")
    return out


def extract_reads_pct(read1, read2, output_dir, pct, seed=42):
    """Extract a percentage of reads from FASTQ files."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    total = count_reads(read1)
    num_reads = max(1, int(total * pct / 100))
    print(f"  Total reads: {total}, sampling {pct}% = {num_reads} reads")
    return extract_reads(read1, read2, output_dir, num_reads, seed)


def _write_fastq(path, records, gzip=False):
    mode = "wt" if gzip else "w"
    opener = gzip.open if gzip else open
    with opener(path, mode) as f:
        for header, seq, qual in records:
            f.write(f"{header}\n{seq}\n+\n{qual}\n")
