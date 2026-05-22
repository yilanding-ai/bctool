import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from bctool.simulation import Simulator
from bctool.methylation.converter import ConversionType


def test_simulation_ct_se():
    """Test C-to-T single-end simulation."""
    with tempfile.TemporaryDirectory() as tmp:
        sim = Simulator(
            genome_length=10000,
            num_reads=100,
            read_length=50,
            methylation_rate=0.7,
            conversion="ct",
            mode="se",
            seed=42,
        )
        out = sim.run(tmp)
        assert (out / "reference.fa").exists()
        assert (out / "ground_truth.bed").exists()
        assert (out / "simulated_R1.fastq").exists()
        print("[PASS] C-to-T single-end simulation")


def test_simulation_ag_pe():
    """Test A-to-G paired-end simulation."""
    with tempfile.TemporaryDirectory() as tmp:
        sim = Simulator(
            genome_length=10000,
            num_reads=100,
            read_length=50,
            methylation_rate=0.6,
            conversion="ag",
            mode="pe",
            seed=42,
        )
        out = sim.run(tmp)
        assert (out / "reference.fa").exists()
        assert (out / "ground_truth.bed").exists()
        assert (out / "simulated_R1.fastq").exists()
        assert (out / "simulated_R2.fastq").exists()
        print("[PASS] A-to-G paired-end simulation")


def test_conversion_type():
    ct = ConversionType.from_str("ct")
    assert ct.target_base == "C"
    assert ct.converted_base == "T"
    ag = ConversionType.from_str("ag")
    assert ag.target_base == "A"
    assert ag.converted_base == "G"
    print("[PASS] Conversion types")


def test_ground_truth_format():
    with tempfile.TemporaryDirectory() as tmp:
        sim = Simulator(
            genome_length=5000,
            num_reads=50,
            read_length=50,
            conversion="ct",
            mode="se",
            seed=42,
        )
        out = sim.run(tmp)
        with open(out / "ground_truth.bed") as f:
            lines = f.readlines()
        assert len(lines) > 0
        for line in lines[:5]:
            parts = line.strip().split("\t")
            assert len(parts) == 4
            assert parts[3] in ("0", "1")
        print(f"[PASS] Ground truth: {len(lines)} sites")


if __name__ == "__main__":
    test_simulation_ct_se()
    test_simulation_ag_pe()
    test_conversion_type()
    test_ground_truth_format()
    print("\nAll tests passed!")
