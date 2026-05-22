"""Quick timing test for ConvSeed on a single scenario."""
import sys, time, tempfile, shutil
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from bctool.simulation import Simulator
from bctool.evaluation.metrics import Evaluator
from bctool.aligners.conversion_aware import ConversionAwareAligner
from bctool.config import Config

cfg = Config()
cfg.data["conversion"] = "ct"
cfg.data["threads"] = 1
cfg.data["output_dir"] = "."
cfg.data["correction"] = {"enabled": False}
cfg.data["mode"] = "se"

tmp = Path(tempfile.mkdtemp())
sim = Simulator(genome_length=50000, num_reads=500, read_length=50,
                methylation_rate=0.7, error_rate=0.01, conversion="ct",
                mode="se", seed=42)
sim_out = sim.run(tmp)
gt = sim_out / "ground_truth.bed"
fq = sim_out / "simulated_R1.fastq"
ref = sim_out / "reference.fa"
evaluator = Evaluator(gt)

for seed_len in [8, 10, 12, 14]:
    cs = ConversionAwareAligner(cfg, tmp, name_suffix=f"s{seed_len}", seed_length=seed_len)
    t0 = time.time()
    cs.build_index(ref)
    t_idx = time.time() - t0
    t0 = time.time()
    bam = cs.run_align(fq, None, ref)
    t_align = time.time() - t0
    t0 = time.time()
    bed = cs.call_methylation(bam, ref)
    t_call = time.time() - t0
    res = evaluator.evaluate(bed, f"seed{seed_len}")
    print(f"seed={seed_len:2d}  F1={res['f1_score']:.4f}  P={res['precision']:.4f}  "
          f"R={res['recall']:.4f}  overlap={res['sites_overlap']}/{res['sites_truth']}  "
          f"idx={t_idx:.2f}s  align={t_align:.2f}s  call={t_call:.2f}s")

shutil.rmtree(tmp)
