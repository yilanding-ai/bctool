import json
import sys
import shutil
from pathlib import Path
from datetime import datetime
from .config import Config
from .extraction import extract_reads, count_reads
from .simulation import Simulator
from .aligners import get_aligner, ALIGNER_REGISTRY
from .methylation.converter import ConversionType
from .methylation.caller import MethylationCaller
from .evaluation.metrics import Evaluator
from .evaluation.comparator import compare_aligners, find_best_aligner


class Pipeline:
    def __init__(self, config: Config):
        self.config = config
        self.output_dir = Path(config.output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.conversion = ConversionType.from_str(config.conversion)
        self.work_dir = self.output_dir / "work"
        self.work_dir.mkdir(exist_ok=True)

    def check_tools(self):
        print("Checking tool availability...")
        available = []
        missing = []
        for name in self.config.enabled_aligners:
            cls = get_aligner(name)
            inst = cls(self.config, self.work_dir)
            if inst.is_available():
                available.append(name)
                print(f"  [OK] {name}")
            else:
                missing.append(name)
                print(f"  [MISSING] {name} (binary: {inst.binary_name})")

        print(f"\nAvailable: {len(available)}/{len(self.config.enabled_aligners)}")
        if missing:
            print(f"Missing: {', '.join(missing)}")
            print("\nInstall missing tools via conda:")
            print("  conda install -c bioconda bwameth bsbolt bsmap walt abismal batmeth2 hisat-3n bismark bsseeker2")
        return available, missing

    def run(self):
        start_time = datetime.now()
        print("=" * 60)
        print(f"BCTool Pipeline")
        print(f"Conversion: {self.conversion.description()}")
        print(f"Mode: {'Paired-end' if self.config.mode == 'pe' else 'Single-end'}")
        print(f"Output: {self.output_dir}")
        print(f"Started: {start_time}")
        print("=" * 60)

        # Phase 0: QC analysis
        read1, read2, reference = self._prepare_inputs()
        self._run_qc(read1, read2)

        # Phase 1.5: Preprocess (adapter trimming + UMI dedup)
        read1, read2 = self._preprocess(read1, read2)

        # Phase 2: Run each aligner
        results = self._run_aligners(read1, read2, reference)

        # Phase 3: Evaluate
        evaluation_results = self._evaluate_all(results, reference)

        # Phase 4: Report
        self._generate_output(evaluation_results)

        elapsed = (datetime.now() - start_time).total_seconds()
        print(f"\nTotal time: {elapsed:.1f}s")
        print(f"Results saved to: {self.output_dir}")

    def _prepare_inputs(self):
        read1 = self.config.read1
        read2 = self.config.read2
        reference = self.config.reference

        # Simulation mode
        if not read1 or not Path(read1).exists():
            print("\n[Phase 1] Generating simulated data...")
            sim_out = self.work_dir / "simulated"
            sim = Simulator(
                genome_length=self.config.get("simulation", "genome_length", default=1_000_000),
                num_reads=self.config.get("simulation", "num_reads", default=50_000),
                read_length=self.config.get("simulation", "read_length", default=100),
                methylation_rate=self.config.get("simulation", "methylation_rate", default=0.7),
                conversion=self.config.conversion,
                mode=self.config.mode,
            )
            sim.run(sim_out)
            read1 = str(sim_out / "simulated_R1.fastq")
            read2 = str(sim_out / "simulated_R2.fastq") if self.config.mode == "pe" else None
            reference = str(sim_out / "reference.fa")

        # Extraction mode
        elif self.config.get("extraction", "enabled", default=True):
            print("\n[Phase 1] Extracting reads...")
            sample_pct = self.config.get("extraction", "sample_pct", default=None)
            seed = self.config.get("extraction", "seed", default=42)
            extract_dir = self.work_dir / "extracted"
            if sample_pct is not None:
                from .extraction import extract_reads_pct
                extract_reads_pct(read1, read2, extract_dir, sample_pct, seed)
            else:
                num_reads = self.config.get("extraction", "num_reads", default=100000)
                extract_reads(read1, read2, extract_dir, num_reads, seed)
            read1 = str(extract_dir / "extracted_R1.fastq")
            if read2:
                read2 = str(extract_dir / "extracted_R2.fastq")

        if not reference:
            print("ERROR: No reference genome provided. Use --reference or --simulate.")
            sys.exit(1)

        return read1, read2, reference

    def _run_qc(self, read1, read2):
        """Phase 0: Run QC analysis on input reads."""
        qc_enabled = self.config.get("qc", "enabled", default=False)
        if not qc_enabled:
            return

        print("\n[Phase 0] Quality Control analysis...")
        from .qc import QCAnalyzer, QCReport
        qc_dir = self.work_dir / "qc"
        qc_dir.mkdir(exist_ok=True)

        analyzer = QCAnalyzer(max_reads=self.config.get("qc", "max_reads", default=50000))
        analyzer.analyze(read1, read2)
        print(analyzer.summary_text())

        report = QCReport(analyzer, str(qc_dir))
        report.generate_html()

    def _preprocess(self, read1, read2):
        """Phase 1.5: Adapter trimming and UMI deduplication."""
        adapter_enabled = self.config.get("adapter", "enabled", default=False)
        umi_enabled = self.config.get("umi", "enabled", default=False)

        if not adapter_enabled and not umi_enabled:
            return read1, read2

        print("\n[Phase 1.5] Preprocessing reads...")
        preproc_dir = self.work_dir / "preprocessed"
        preproc_dir.mkdir(exist_ok=True)

        # Step 1: Adapter trimming
        if adapter_enabled:
            print("  Trimming adapters...")
            from .preprocessing import AdapterTrimmer
            trimmer = AdapterTrimmer(self.config, preproc_dir / "trimmed")

            if self.config.get("adapter", "auto_detect", default=True):
                detected = trimmer.auto_detect_adapters(read1)
                if detected:
                    trimmer.trim(read1, read2, adapters=detected)
                else:
                    print("  [Adapter] No adapters detected, skipping trimming")
            else:
                adapters = {
                    "r1_3prime": self.config.get("adapter", "r1_3prime"),
                    "r2_3prime": self.config.get("adapter", "r2_3prime"),
                }
                trimmer.trim(read1, read2, adapters=adapters)

            # Update read paths to trimmed files
            trimmed_r1 = preproc_dir / "trimmed" / "trimmed_R1.fastq"
            trimmed_r2 = preproc_dir / "trimmed" / "trimmed_R2.fastq" if read2 else None
            if trimmed_r1.exists():
                read1 = str(trimmed_r1)
            if trimmed_r2 and trimmed_r2.exists():
                read2 = str(trimmed_r2)

        # Step 2: UMI deduplication
        if umi_enabled:
            print("  Processing UMIs...")
            from .preprocessing import UmiProcessor
            umi_proc = UmiProcessor(self.config, preproc_dir / "umi")
            umi_config = {
                "location": self.config.get("umi", "location", default="header"),
                "length": self.config.get("umi", "length", default=0),
                "pattern": self.config.get("umi", "pattern", default=None),
                "strategy": self.config.get("umi", "strategy", default="cluster"),
            }
            umi_proc.extract_and_deduplicate(read1, read2, umi_config)
            umi_proc.report_stats()

            dedup_r1 = preproc_dir / "umi" / "dedup_R1.fastq"
            dedup_r2 = preproc_dir / "umi" / "dedup_R2.fastq" if read2 else None
            if dedup_r1.exists():
                read1 = str(dedup_r1)
            if dedup_r2 and dedup_r2.exists():
                read2 = str(dedup_r2)

        return read1, read2

    def _run_aligners(self, read1, read2, reference):
        print("\n[Phase 2] Running alignment tools...")
        results = {}

        real_aligners = []
        for name in self.config.enabled_aligners:
            if name.startswith("mock_"):
                continue
            cls = get_aligner(name)
            aligner = cls(self.config, self.work_dir)
            if aligner.is_available():
                real_aligners.append((name, cls))
            else:
                print(f"  SKIPPED: {name} not installed (binary: {aligner.binary_name})")

        # Fall back to mock aligners if no real tools available
        if not real_aligners:
            print("  No external tools available â€?using built-in mock aligners for demonstration.")
            from .aligners.mock import MockAlignerFactory
            mocks = MockAlignerFactory.create_mocks(self.config, self.work_dir)
            for mock in mocks:
                print(f"\n  --- {mock.name} ---")
                bam = mock.run_align(read1, read2, reference, threads=self.config.threads)
                meth_bed = mock.call_methylation(bam, reference)
                results[mock.name] = {"bam": bam, "meth_bed": meth_bed}
            return results

        correction_enabled = self.config.get("correction", "enabled", default=False)

        # Train model-based error detector if "model" strategy is enabled
        model_trained = False
        if correction_enabled:
            sam_strategies_all = self.config.get("correction", "strategies", default=[])
            if "model" in sam_strategies_all:
                ground_truth = self.work_dir / "simulated" / "ground_truth.bed"
                if ground_truth.exists():
                    from .preprocessing import ErrorCorrector
                    train_corrector = ErrorCorrector(self.config, self.work_dir / "corrected")
                    # Use the first mock aligner SAM for training
                    from .aligners.mock import MockAlignerFactory
                    mock_train = MockAlignerFactory.create_mocks(self.config, self.work_dir)[:1]
                    for mock in mock_train:
                        train_sam = mock.run_align(read1, read2, reference, threads=self.config.threads)
                        train_corrector.train_model(train_sam, ground_truth, self.conversion)
                        model_trained = True
                        break
                else:
                    print("  [Pipeline] Ground truth not found, skipping model training")

        for name, cls in real_aligners:
            print(f"\n  --- {name} ---")
            aligner = cls(self.config, self.work_dir)

            if aligner.requires_index:
                print(f"  Building index...")
                aligner.build_index(reference)

            print(f"  Aligning...")
            bam = aligner.run_align(read1, read2, reference, threads=self.config.threads)

            # SAM-level error correction (real aligners only)
            meth_bed = None
            if correction_enabled:
                sam_strategies = [s for s in self.config.get("correction", "strategies",
                                                              default=["mq", "clip"])
                                  if s in ("mq", "clip", "pair", "unconverted", "model")]
                if sam_strategies:
                    from .preprocessing import ErrorCorrector
                    corrector = ErrorCorrector(self.config, self.work_dir / "corrected")
                    if model_trained:
                        corrector.detector = train_corrector.detector
                    corrected_sam = corrector.correct_sam(
                        bam, self.conversion, strategies=sam_strategies
                    )
                    bam = corrected_sam

            print(f"  Calling methylation...")
            meth_bed = aligner.call_methylation(bam, reference)
            results[name] = {"bam": bam, "meth_bed": meth_bed}

        return results

    def _evaluate_all(self, aligner_results, reference):
        print("\n[Phase 3] Evaluating accuracy...")

        ground_truth = self.work_dir / "simulated" / "ground_truth.bed"
        if not ground_truth.exists():
            for f in self.work_dir.rglob("ground_truth.bed"):
                ground_truth = f
                break

        # Apply error correction on methylation BED files if enabled
        correction_enabled = self.config.get("correction", "enabled", default=False)
        if correction_enabled:
            print("  Applying error correction to methylation calls...")
            from .preprocessing import ErrorCorrector
            corrector = ErrorCorrector(self.config, self.work_dir / "corrected")
            corr_strategies = self.config.get("correction", "strategies",
                                              default=["consensus"])
            corr_min_depth = self.config.get("correction", "min_depth", default=5)
            corr_agree = self.config.get("correction", "min_agree_ratio", default=0.7)

            for name in list(aligner_results.keys()):
                meth_bed = aligner_results[name].get("meth_bed")
                if meth_bed and Path(meth_bed).exists():
                    corrected_bed = corrector.correct_bed(
                        meth_bed, self.conversion, reference=reference,
                        min_depth=corr_min_depth,
                        min_agree_ratio=corr_agree
                    )
                    aligner_results[name]["meth_bed"] = corrected_bed

        evaluator = Evaluator(str(ground_truth))
        all_results = []

        for name, paths in aligner_results.items():
            print(f"  Evaluating {name}...")
            result = evaluator.evaluate(paths["meth_bed"], name)
            all_results.append(result)
            if result["status"] == "ok":
                print(f"    F1={result['f1_score']:.4f}  P={result['precision']:.4f}  "
                      f"R={result['recall']:.4f}  A={result['accuracy']:.4f}")
            else:
                print(f"    Status: {result['status']}")

        comparison = compare_aligners(all_results)
        best = find_best_aligner(comparison)
        if best:
            print(f"\n{'=' * 50}")
            print(f"BEST ALIGNER: {best['aligner']}")
            print(f"  F1 Score: {best['f1_score']:.4f}")
            print(f"  Precision: {best['precision']:.4f}")
            print(f"  Recall: {best['recall']:.4f}")
            print(f"{'=' * 50}")

        return {"comparison": comparison, "best": best, "results": all_results}

    def _generate_output(self, eval_data):
        print("\n[Phase 4] Generating output...")

        # JSON output - save both ranked and raw results
        df = eval_data["comparison"]
        ranked_results = df.to_dict("records") if not df.empty else []
        output = {
            "config": {
                "conversion": self.config.conversion,
                "mode": self.config.mode,
            },
            "summary": eval_data["best"],
            "results": ranked_results,
        }

        json_path = self.output_dir / "comparison_results.json"
        with open(json_path, "w") as f:
            json.dump(output, f, indent=2)

        # CSV output
        csv_path = self.output_dir / "comparison_results.csv"
        df.to_csv(csv_path, index=False)

        # Generate report
        from .reporting.report import generate_report
        generate_report(self.output_dir)

        # Print ranking
        print(f"\nFinal Ranking:")
        print(f"{'Rank':<6} {'Aligner':<25} {'F1':<10} {'Precision':<10} {'Recall':<10} {'Accuracy':<10}")
        print("-" * 70)
        for _, row in df.iterrows():
            print(f"{int(row['rank']):<6} {row['aligner']:<25} "
                  f"{row['f1_score']:<10.4f} {row['precision']:<10.4f} "
                  f"{row['recall']:<10.4f} {row['accuracy']:<10.4f}")

        best = eval_data["best"]
        if best:
            print(f"\n{'=' * 60}")
            print(f"  RECOMMENDATION: Use {best['aligner']} for highest accuracy")
            print(f"  (F1={best['f1_score']:.4f}, Precision={best['precision']:.4f}, "
                  f"Recall={best['recall']:.4f})")
            print(f"{'=' * 60}")
