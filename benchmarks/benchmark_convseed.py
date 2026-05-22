"""Benchmark ConvSeed accuracy + speed across conversion types and read lengths."""

import sys, time, json, tempfile, shutil
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from bctool.simulation import Simulator
from bctool.methylation.converter import ConversionType
from bctool.methylation.caller import MethylationCaller
from bctool.evaluation.metrics import Evaluator
from bctool.aligners.conversion_aware import ConversionAwareAligner
from bctool.aligners.mock import MockAligner
from bctool.config import Config

RESULTS = {}
GENOME_LENGTH = 50000
NUM_READS = 5000

def make_config(conversion):
    cfg = Config()
    cfg.data["conversion"] = conversion
    cfg.data["threads"] = 1
    cfg.data["output_dir"] = "."
    cfg.data["correction"] = {"enabled": False}
    cfg.data["mode"] = "se"
    return cfg

def run_benchmark():

    conversions = ["ct", "ag", "ga", "tc"]
    read_lengths = [50, 100, 150]
    error_rate = 0.01

    for conv in conversions:
        for rl in read_lengths:
            key = f"{conv}_rl{rl}"
            print(f"\n{'='*60}")
            print(f"  {key}")
            print(f"{'='*60}")

            tmp = Path(tempfile.mkdtemp(prefix=f"bm_{key}_"))
            cfg = make_config(conv)

            # Simulate
            t0 = time.time()
            sim = Simulator(
                genome_length=GENOME_LENGTH, num_reads=NUM_READS,
                read_length=rl, methylation_rate=0.7,
                error_rate=error_rate, conversion=conv, mode="se", seed=42,
            )
            sim_out = sim.run(tmp)
            sim_time = time.time() - t0
            ground_truth = sim_out / "ground_truth.bed"
            fastq = sim_out / "simulated_R1.fastq"
            ref = sim_out / "reference.fa"
            evaluator = Evaluator(ground_truth)

            # ConvSeed seed=10 (optimal for 50bp after fix)
            tag = "conv_seed10"
            cs = ConversionAwareAligner(cfg, tmp, name_suffix="seed10", seed_length=10, auto_adapt=False)
            t0 = time.time()
            cs.build_index(ref)
            t_idx = time.time() - t0
            t0 = time.time()
            bam = cs.run_align(fastq, None, ref)
            t_align = time.time() - t0
            t0 = time.time()
            bed = cs.call_methylation(bam, ref)
            t_call = time.time() - t0
            res = evaluator.evaluate(bed, tag)
            res["index_time"] = round(t_idx, 3)
            res["align_time"] = round(t_align, 3)
            res["call_time"] = round(t_call, 3)
            res["sim_time"] = round(sim_time, 3)
            RESULTS.setdefault(key, {})[tag] = res
            print(f"  {tag:15s} F1={res['f1_score']:.4f}  P={res['precision']:.4f}  "
                  f"R={res['recall']:.4f}  align={t_align:.2f}s  idx={t_idx:.2f}s  "
                  f"overlap={res.get('sites_overlap',0)}/{res.get('sites_truth',0)}")

            # ConvSeed seed=12 (shorter seed = more sensitive)
            tag = "conv_seed12"
            cs12 = ConversionAwareAligner(cfg, tmp, name_suffix="seed12", seed_length=12, auto_adapt=False)
            t0 = time.time()
            cs12.build_index(ref)
            t_idx12 = time.time() - t0
            t0 = time.time()
            bam12 = cs12.run_align(fastq, None, ref)
            t_align12 = time.time() - t0
            t0 = time.time()
            bed12 = cs12.call_methylation(bam12, ref)
            t_call12 = time.time() - t0
            res12 = evaluator.evaluate(bed12, tag)
            res12["index_time"] = round(t_idx12, 3)
            res12["align_time"] = round(t_align12, 3)
            res12["call_time"] = round(t_call12, 3)
            RESULTS.setdefault(key, {})[tag] = res12
            print(f"  {tag:15s} F1={res12['f1_score']:.4f}  P={res12['precision']:.4f}  "
                  f"R={res12['recall']:.4f}  align={t_align12:.2f}s  idx={t_idx12:.2f}s  "
                  f"overlap={res12.get('sites_overlap',0)}/{res12.get('sites_truth',0)}")

            # Mock aligners (for reference: genome-level methylation simulation)
            for level, acc in [("perfect", 0.98), ("good", 0.88),
                               ("medium", 0.75), ("poor", 0.55)]:
                mock = MockAligner(cfg, tmp, accuracy=acc, name_suffix=level)
                t0 = time.time()
                mock.run_align(fastq, None, ref)
                t_align_m = time.time() - t0
                t0 = time.time()
                mbed = mock.call_methylation(None, ref)
                t_call_m = time.time() - t0
                mr = evaluator.evaluate(mbed, f"mock_{level}")
                mr["align_time"] = round(t_align_m, 3)
                mr["call_time"] = round(t_call_m, 3)
                RESULTS.setdefault(key, {})[f"mock_{level}"] = mr
                print(f"  mock_{level:10s} F1={mr['f1_score']:.4f}  P={mr['precision']:.4f}  "
                      f"R={mr['recall']:.4f}")

            shutil.rmtree(tmp)

    return RESULTS

def print_summary(results):
    print(f"\n\n{'='*90}")
    print("ConvSeed BENCHMARK SUMMARY")
    print(f"{'='*90}")
    # Per-scenario comparison
    for scenario in sorted(results.keys()):
        conv, rl = scenario.split("_rl")
        data = results[scenario]
        print(f"\n{conv.upper()}  read_len={rl}bp")
        print(f"  {'Method':20s} {'F1':>8s} {'Prec':>8s} {'Rec':>8s} {'Overlap':>8s} {'Align(s)':>10s}")
        print(f"  {'-'*20} {'-'*8} {'-'*8} {'-'*8} {'-'*8} {'-'*10}")
        for method in sorted(data.keys()):
            d = data[method]
            print(f"  {method:20s} {d['f1_score']:8.4f} {d['precision']:8.4f} "
                  f"{d['recall']:8.4f} {d.get('sites_overlap',0):8d} {d.get('align_time',0):10.2f}")

    # Aggregate table: ConvSeed across all scenarios
    for seed_tag in ["conv_seed10", "conv_seed12"]:
        print(f"\n\nConvSeed({seed_tag}) AGGREGATE")
        print(f"  {'Scenario':15s} {'F1':>8s} {'Prec':>8s} {'Rec':>8s} {'Overlap/Tot':>12s} {'Idx(s)':>8s} {'Align(s)':>10s}")
        print(f"  {'-'*15} {'-'*8} {'-'*8} {'-'*8} {'-'*12} {'-'*8} {'-'*10}")
        for scenario in sorted(results.keys()):
            d = results[scenario].get(seed_tag, {})
            ot = f"{d.get('sites_overlap',0)}/{d.get('sites_truth',0)}"
            print(f"  {scenario:15s} {d.get('f1_score',0):8.4f} {d.get('precision',0):8.4f} "
                  f"{d.get('recall',0):8.4f} {ot:>12s} {d.get('index_time',0):8.3f} {d.get('align_time',0):10.2f}")

if __name__ == "__main__":
    results = run_benchmark()
    print_summary(results)
    out = Path(__file__).resolve().parent / "benchmark_results.json"
    with open(out, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nResults saved to {out}")
