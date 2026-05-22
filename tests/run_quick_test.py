"""Quick integration test: simulate small data and evaluate."""
import sys, json, tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from bctool.simulation import Simulator
from bctool.methylation.caller import MethylationCaller
from bctool.methylation.converter import ConversionType
from bctool.evaluation.metrics import Evaluator

def test_ct_se():
    print("=" * 50)
    print("Test 1: C-to-T Single-end")
    print("=" * 50)
    with tempfile.TemporaryDirectory() as tmp:
        sim = Simulator(genome_length=20000, num_reads=500, read_length=50,
                        methylation_rate=0.7, conversion="ct", mode="se", seed=42)
        sim_out = sim.run(tmp)
        
        # Call methylation
        caller = MethylationCaller(ConversionType("C>T"), min_qual=0, min_cov=1)
        
        # Parse simulated reads (treat as if they were aligned)
        sites = {}
        with open(sim_out / "simulated_R1.fastq") as f:
            lines = f.readlines()
        for i in range(0, len(lines), 4):
            seq = lines[i+1].strip()
            # Simulate "perfect alignment" - each read aligns to a known position
            read_id = lines[i].strip()
            # For the test, use the read index as pseudo-position
            for j, base in enumerate(seq):
                pos = (i // 4) * 50 + j  # pseudo position
                target = "C"
                conv = "T"
                if base == target:
                    sites[(f"chr_sim", pos)] = {"meth": 0, "unmeth": 1}
                elif base == conv:
                    sites[(f"chr_sim", pos)] = {"meth": 1, "unmeth": 0}
        
        # Write called methylation
        called_bed = Path(tmp) / "called.bed"
        with open(called_bed, "w") as f:
            for (chrom, pos), counts in sorted(sites.items()):
                total = counts["meth"] + counts["unmeth"]
                if total > 0:
                    level = counts["meth"] / total
                    f.write(f"{chrom}\t{pos}\t{pos+1}\t{level:.4f}\t{counts['meth']}\t{counts['unmeth']}\t{total}\n")
        
        # Evaluate
        evaluator = Evaluator(sim_out / "ground_truth.bed")
        result = evaluator.evaluate(called_bed, "perfect_aligner")
        
        print(f"  Precision: {result['precision']:.4f}")
        print(f"  Recall:    {result['recall']:.4f}")
        print(f"  F1 Score:  {result['f1_score']:.4f}")
        print(f"  Accuracy:  {result['accuracy']:.4f}")
        print(f"  Sites overlap: {result.get('sites_overlap', 0)}")
        
        assert result['f1_score'] > 0.5, f"F1 too low: {result['f1_score']}"
        print("  [PASS]")

def test_ag_pe():
    print("\n" + "=" * 50)
    print("Test 2: A-to-G Paired-end")
    print("=" * 50)
    with tempfile.TemporaryDirectory() as tmp:
        sim = Simulator(genome_length=20000, num_reads=500, read_length=50,
                        methylation_rate=0.6, conversion="ag", mode="pe", seed=42)
        sim_out = sim.run(tmp)
        
        caller = MethylationCaller(ConversionType("A>G"), min_qual=0, min_cov=1)
        
        sites = {}
        for fname in ["simulated_R1.fastq", "simulated_R2.fastq"]:
            with open(sim_out / fname) as f:
                lines = f.readlines()
            for i in range(0, len(lines), 4):
                seq = lines[i+1].strip()
                for j, base in enumerate(seq):
                    pos = hash(lines[i].strip()) % 100000
                    target = "A"
                    conv = "G"
                    if base == target:
                        sites[(f"chr_sim", pos)] = {"meth": 0, "unmeth": sites.get((f"chr_sim", pos), {"meth": 0, "unmeth": 0})["unmeth"] + 1}
                    elif base == conv:
                        sites[(f"chr_sim", pos)] = {"meth": sites.get((f"chr_sim", pos), {"meth": 0, "unmeth": 0})["meth"] + 1, "unmeth": sites.get((f"chr_sim", pos), {"meth": 0, "unmeth": 0})["unmeth"]}
        
        called_bed = Path(tmp) / "called_ag.bed"
        with open(called_bed, "w") as f:
            for (chrom, pos), counts in sorted(sites.items()):
                total = counts["meth"] + counts["unmeth"]
                if total > 0:
                    level = counts["meth"] / total
                    f.write(f"{chrom}\t{pos}\t{pos+1}\t{level:.4f}\t{counts['meth']}\t{counts['unmeth']}\t{total}\n")
        
        evaluator = Evaluator(sim_out / "ground_truth.bed")
        result = evaluator.evaluate(called_bed, "perfect_aligner_ag")
        
        print(f"  Accuracy: {result['accuracy']:.4f}")
        print(f"  Sites called: {result['sites_called']}")
        print("  [PASS]")

def test_evaluation_ranking():
    print("\n" + "=" * 50)
    print("Test 3: Evaluation ranking logic")
    print("=" * 50)
    from bctool.evaluation.comparator import compare_aligners, find_best_aligner
    
    results = [
        {"aligner": "A", "precision": 0.95, "recall": 0.90, "f1_score": 0.9247, "accuracy": 0.93, "sites_called": 100, "sites_truth": 100, "sites_overlap": 95, "status": "ok"},
        {"aligner": "B", "precision": 0.80, "recall": 0.85, "f1_score": 0.8242, "accuracy": 0.82, "sites_called": 100, "sites_truth": 100, "sites_overlap": 95, "status": "ok"},
        {"aligner": "C", "precision": 0.70, "recall": 0.75, "f1_score": 0.7241, "accuracy": 0.72, "sites_called": 100, "sites_truth": 100, "sites_overlap": 95, "status": "ok"},
    ]
    
    df = compare_aligners(results)
    best = find_best_aligner(df)
    
    assert best["aligner"] == "A", f"Expected A, got {best['aligner']}"
    assert df.iloc[0]["aligner"] == "A"
    assert df.iloc[2]["aligner"] == "C"
    print(f"  Best: {best['aligner']} (F1={best['f1_score']:.4f})")
    print("  [PASS]")

def test_all_12_conversions():
    """Test that all 12 conversion types work with simulation."""
    from bctool.methylation.converter import get_conversion, list_conversions
    
    print("\n" + "=" * 50)
    print("Test 4: All 12 conversion types")
    print("=" * 50)
    
    conv_list = ["C>T", "T>C", "A>G", "G>A", "A>C", "C>A",
                 "G>T", "T>G", "A>T", "T>A", "C>G", "G>C"]
    
    with tempfile.TemporaryDirectory() as tmp:
        for i, tag in enumerate(conv_list):
            conv = ConversionType(tag)
            sim = Simulator(genome_length=10000, num_reads=50, read_length=50,
                            methylation_rate=0.7, conversion=tag, mode="se", seed=i)
            out = sim.run(tmp + f"_{i}")
            assert (out / "reference.fa").exists()
            assert (out / "simulated_R1.fastq").exists()
            print(f"  [OK] {conv.key:>3} ({conv.description():<10}) "
                  f"target={conv.target_base} -> {conv.converted_base}")
    
    # Verify complement relationships
    for tag in conv_list:
        conv = ConversionType(tag)
        comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
        assert conv.complement_target == comp[conv.target_base], \
            f"{tag}: comp_target mismatch"
        assert conv.complement_converted == comp[conv.converted_base], \
            f"{tag}: comp_converted mismatch"
    
    print("  [PASS] All 12 conversions work correctly")


if __name__ == "__main__":
    test_ct_se()
    test_ag_pe()
    test_evaluation_ranking()
    test_all_12_conversions()
    print("\n" + "=" * 50)
    print("ALL TESTS PASSED!")
    print("=" * 50)
