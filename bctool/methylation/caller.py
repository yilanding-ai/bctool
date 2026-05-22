import re
from pathlib import Path
from collections import defaultdict
from .converter import ConversionType


class MethylationCaller:
    """Call methylation status from aligned SAM/BAM data.
    
    This is a universal methylation caller used when the aligner
    does not provide its own, or for cross-validation.
    """

    def __init__(self, conversion: ConversionType, min_qual=20, min_cov=5):
        self.conversion = conversion
        self.min_qual = min_qual
        self.min_cov = min_cov

    def parse_sam(self, sam_path):
        """Parse SAM file and extract methylation information."""
        sites = defaultdict(lambda: {"meth": 0, "unmeth": 0})

        with open(sam_path) as f:
            for line in f:
                if line.startswith("@"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 11:
                    continue

                flag = int(parts[1])
                chrom = parts[2]
                pos = int(parts[3]) - 1
                seq = parts[9]
                qual = parts[10]
                cigar = parts[5]

                if flag & 0x4:
                    continue

                # Parse CIGAR to get aligned positions
                aligned_positions = self._cigar_to_positions(pos, cigar, len(seq))

                for read_idx, (ref_pos, read_base) in enumerate(aligned_positions):
                    if read_idx >= len(qual):
                        continue
                    if ord(qual[read_idx]) - 33 < self.min_qual:
                        continue
                    if ref_pos < 0:
                        continue

                    target = self.conversion.target_base
                    conv = self.conversion.converted_base
                    comp_target = self.conversion.complement_target
                    comp_conv = self.conversion.complement_converted

                    if read_base == target:
                        sites[(chrom, ref_pos)]["unmeth"] += 1
                    elif read_base == conv:
                        sites[(chrom, ref_pos)]["meth"] += 1
                    # Also check complement strand
                    rc_base = self._rc(read_base)
                    if rc_base == comp_target:
                        sites[(chrom, ref_pos)]["unmeth"] += 1
                    elif rc_base == comp_conv:
                        sites[(chrom, ref_pos)]["meth"] += 1

        return sites

    def to_bed(self, sites, output_path):
        with open(output_path, "w") as f:
            for (chrom, pos), counts in sorted(sites.items()):
                total = counts["meth"] + counts["unmeth"]
                if total >= self.min_cov:
                    level = counts["meth"] / total
                    f.write(f"{chrom}\t{pos}\t{pos+1}\t"
                           f"{level:.4f}\t{counts['meth']}\t{counts['unmeth']}\t{total}\n")
        return output_path

    @staticmethod
    def _cigar_to_positions(pos, cigar, seq_len):
        ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
        ref_pos = pos
        read_pos = 0
        result = []

        for length, op in ops:
            length = int(length)
            if op == "M":
                for i in range(length):
                    if read_pos < seq_len:
                        result.append((ref_pos, seq[read_pos]))
                        ref_pos += 1
                        read_pos += 1
            elif op == "D":
                ref_pos += length
            elif op == "I":
                read_pos += length
            elif op == "S":
                read_pos += length
            elif op == "N":
                ref_pos += length
            elif op == "=":
                for i in range(length):
                    if read_pos < seq_len:
                        result.append((ref_pos, seq[read_pos]))
                        ref_pos += 1
                        read_pos += 1
            elif op == "X":
                for i in range(length):
                    if read_pos < seq_len:
                        result.append((ref_pos, seq[read_pos]))
                        ref_pos += 1
                        read_pos += 1
        return result

    @staticmethod
    def _rc(base):
        comp = {"A": "T", "T": "A", "C": "G", "G": "C",
                "a": "t", "t": "a", "c": "g", "g": "c"}
        return comp.get(base, base)
