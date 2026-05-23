from collections import defaultdict
import re
from .converter import ConversionType


class PerBaseCaller:
    """Extract per-position methylation data with quality scores from SAM/BAM."""

    def __init__(self, conversion: ConversionType, min_qual=0):
        self.conversion = conversion
        self.min_qual = min_qual

    def parse_sam(self, sam_path):
        sites = defaultdict(lambda: {
            "+": {"conv_quals": [], "unconv_quals": []},
            "-": {"conv_quals": [], "unconv_quals": []},
        })

        target = self.conversion.target_base
        conv = self.conversion.converted_base
        comp_target = self.conversion.complement_target
        comp_conv = self.conversion.complement_converted

        with open(sam_path) as f:
            for line in f:
                if line.startswith("@"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 11:
                    continue

                flag = int(parts[1])
                chrom = parts[2]
                pos = int(parts[3])
                seq = parts[9]
                qual = parts[10]
                cigar = parts[5]

                if flag & 0x4:
                    continue

                strand = "-" if flag & 0x10 else "+"

                aligned_positions = self._cigar_to_positions(pos - 1, cigar, seq)

                for read_idx, (ref_pos, read_base) in enumerate(aligned_positions):
                    if read_idx >= len(qual):
                        continue
                    qscore = ord(qual[read_idx]) - 33
                    if qscore < self.min_qual:
                        continue
                    if ref_pos < 0:
                        continue

                    key = (chrom, ref_pos + 1)

                    if strand == "+":
                        if read_base == target:
                            sites[key][strand]["unconv_quals"].append(str(qscore))
                        elif read_base == conv:
                            sites[key][strand]["conv_quals"].append(str(qscore))
                    else:
                        if read_base == comp_target:
                            sites[key][strand]["unconv_quals"].append(str(qscore))
                        elif read_base == comp_conv:
                            sites[key][strand]["conv_quals"].append(str(qscore))

        return sites

    def to_csv(self, sites, output_path):
        with open(output_path, "w") as f:
            f.write("ref,pos,strand,convertedBaseQualities,convertedBaseCount,"
                    "unconvertedBaseQualities,unconvertedBaseCount,未转化率\n")
            for (chrom, pos), strands in sorted(sites.items()):
                for strand in ["+", "-"]:
                    data = strands[strand]
                    if not data["conv_quals"] and not data["unconv_quals"]:
                        continue
                    conv_quals = "".join(data["conv_quals"])
                    unconv_quals = "".join(data["unconv_quals"])
                    conv_count = len(data["conv_quals"])
                    unconv_count = len(data["unconv_quals"])
                    total = conv_count + unconv_count
                    unconv_rate = unconv_count / total if total > 0 else 0.0
                    f.write(f"{chrom},{pos},{strand},{conv_quals},{conv_count},"
                            f"{unconv_quals},{unconv_count},{unconv_rate:.6f}\n")
        return output_path

    @staticmethod
    def _cigar_to_positions(pos, cigar, seq):
        ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
        ref_pos = pos
        read_pos = 0
        result = []

        for length, op in ops:
            length = int(length)
            if op in ("M", "=", "X"):
                for i in range(length):
                    if read_pos < len(seq):
                        result.append((ref_pos, seq[read_pos]))
                        ref_pos += 1
                        read_pos += 1
            elif op in ("D", "N"):
                ref_pos += length
            elif op in ("I", "S"):
                read_pos += length
        return result

    @staticmethod
    def _rc(base):
        comp = {"A": "T", "T": "A", "C": "G", "G": "C",
                "a": "t", "t": "a", "c": "g", "g": "c"}
        return comp.get(base, base)
