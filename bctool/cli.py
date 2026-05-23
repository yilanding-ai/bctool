import click
import sys
from pathlib import Path
from .config import Config
from .pipeline import Pipeline
from .methylation.converter import ALL_CONVERSION_TAGS

CONV_CHOICES = ALL_CONVERSION_TAGS


@click.group()
@click.version_option(version="1.0.0")
def main():
    """BCTool: Multi-algorithm methylation alignment benchmark.
    
    Supports all 12 conversion types, single-end and paired-end sequencing data.
    """


@main.command()
@click.option("-1", "--read1", required=False, default=None, help="Read 1 FASTQ file")
@click.option("-2", "--read2", help="Read 2 FASTQ file (for paired-end)")
@click.option("-r", "--reference", help="Reference genome FASTA")
@click.option("-o", "--output", default="./bctool_results", help="Output directory")
@click.option("--config", "-c", "config_path", help="Config file path")
@click.option("--conversion", type=click.Choice(CONV_CHOICES), default="ct",
              help="Conversion type")
@click.option("--num-reads", type=int, default=100000, help="Number of reads to extract from real data")
@click.option("--sample-pct", type=click.FloatRange(0, 100), default=None,
              help="Percentage of reads to use for validation (alternative to --num-reads)")
@click.option("--simulate", is_flag=True, help="Generate simulated data instead")
@click.option("--sim-reads", type=int, default=50000, help="Number of simulated reads")
@click.option("--threads", type=int, default=8, help="Number of threads")
@click.option("--aligners", help="Comma-separated list of aligners (e.g. 'bwameth,bsbolt') or 'all' for all 15")
@click.option("--dry-run", is_flag=True, help="Check tools without running")
@click.option("--demo", is_flag=True, help="Quick demo with mock aligners (no external tools needed)")
@click.option("--trim-adapters", is_flag=True, help="Enable adapter trimming")
@click.option("--adapter1", help="Adapter sequence for read 1 (3' end)")
@click.option("--adapter2", help="Adapter sequence for read 2 (3' end)")
@click.option("--umi", is_flag=True, help="Enable UMI deduplication")
@click.option("--umi-location", type=click.Choice(["header", "inline_5prime", "inline_3prime"]),
              default="header", help="UMI location")
@click.option("--umi-length", type=int, default=0, help="UMI length (0=auto)")
@click.option("--qc", is_flag=True, help="Run QC analysis before alignment")
@click.option("--error-correct", is_flag=True, default=False,
              help="Enable error correction (filter low-quality alignments + consensus)")
@click.option("--correction-strategies", default="mq,clip,consensus",
              help="Correction strategies: mq,clip,pair,consensus (comma-separated)")
@click.option("--correction-min-mq", type=int, default=20,
              help="Min mapping quality for --error-correct")
@click.option("--correction-min-depth", type=int, default=5,
              help="Min depth for consensus correction")
@click.option("--correction-max-unconverted", type=int, default=3,
              help="Max allowed unconverted target bases per read before rejection")
@click.option("--correction-motifs", default=None,
              help="Comma-separated IUPAC motifs for context-aware correction (e.g. CG,AAA)")
@click.option("--correction-no-auto-motif", is_flag=True, default=False,
              help="Disable auto-detection of motifs from data")
@click.option("--correction-motif-flank", type=int, default=2,
              help="Flanking bases around each site for motif extraction")
@click.option("--read-length", type=int, default=150,
              help="Sequencing read length (bp) for algorithm selection")
def run(read1, read2, reference, output, config_path, conversion,
        num_reads, sample_pct, simulate, sim_reads, threads, aligners, dry_run, demo,
        trim_adapters, adapter1, adapter2, umi, umi_location, umi_length, qc,
        error_correct, correction_strategies, correction_min_mq, correction_min_depth,
        correction_max_unconverted, correction_motifs, correction_no_auto_motif,
        correction_motif_flank, read_length):
    """Run the full benchmarking pipeline."""
    cfg = Config(config_path)
    
    if conversion:
        cfg.data["conversion"] = conversion
    if output:
        cfg.data["output_dir"] = output
    if threads:
        cfg.data["threads"] = threads
    if sample_pct is not None:
        cfg.data["extraction"]["sample_pct"] = sample_pct
    elif num_reads:
        cfg.data["extraction"]["num_reads"] = num_reads

    if demo or simulate:
        cfg.data["simulation"]["num_reads"] = sim_reads or 20000
        cfg.data["extraction"]["enabled"] = False

    if read1:
        cfg.data["input"]["read1"] = read1
    if read2:
        cfg.data["input"]["read2"] = read2
        cfg.data["mode"] = "pe"
    elif read1:
        cfg.data["mode"] = "se"
    if reference:
        cfg.data["input"]["reference"] = reference
    if aligners:
        if aligners.lower() == "all":
            # Use all 14 aligners from default config
            pass
        else:
            cfg.data["aligners"]["enabled"] = [a.strip() for a in aligners.split(",")]
    else:
        cfg.data["aligners"]["enabled"] = []  # fall back to mock

    # Preprocessing options
    if trim_adapters:
        cfg.data["adapter"]["enabled"] = True
    if adapter1:
        cfg.data["adapter"]["r1_3prime"] = adapter1
    if adapter2:
        cfg.data["adapter"]["r2_3prime"] = adapter2
    if umi:
        cfg.data["umi"]["enabled"] = True
    if umi_location:
        cfg.data["umi"]["location"] = umi_location
    if umi_length:
        cfg.data["umi"]["length"] = umi_length
    if qc:
        cfg.data["qc"]["enabled"] = True

    # Error correction options
    if error_correct:
        cfg.data["correction"]["enabled"] = True
        if correction_strategies:
            cfg.data["correction"]["strategies"] = [
                s.strip() for s in correction_strategies.split(",")
            ]
        if correction_min_mq:
            cfg.data["correction"]["min_mq"] = correction_min_mq
        if correction_min_depth:
            cfg.data["correction"]["min_depth"] = correction_min_depth
        if correction_max_unconverted is not None:
            cfg.data["correction"]["max_unconverted"] = correction_max_unconverted
        if correction_motifs:
            cfg.data["correction"]["motifs"] = [
                m.strip() for m in correction_motifs.split(",")
            ]
        if correction_no_auto_motif:
            cfg.data["correction"]["auto_detect_motifs"] = False
        if correction_motif_flank is not None:
            cfg.data["correction"]["motif_flank"] = correction_motif_flank
        if read_length:
            cfg.data["correction"]["read_length"] = read_length

    pipe = Pipeline(cfg)
    if dry_run:
        pipe.check_tools()
    else:
        pipe.run()


@main.command()
@click.argument("read1")
@click.argument("read2", required=False)
@click.option("-n", "--num-reads", type=int, default=100000)
@click.option("--seed", type=int, default=42)
@click.option("-o", "--output", default="./extracted")
def extract(read1, read2, num_reads, seed, output):
    """Extract a subset of reads from FASTQ files."""
    from .extraction import extract_reads
    extract_reads(read1, read2, output, num_reads, seed)
    click.echo(f"Extracted {num_reads} reads to {output}")


@main.command()
@click.option("-o", "--output", default="./simulated")
@click.option("--genome-length", type=int, default=1000000)
@click.option("--num-reads", type=int, default=50000)
@click.option("--read-length", type=int, default=100)
@click.option("--methylation-rate", type=float, default=0.7)
@click.option("--conversion", type=click.Choice(CONV_CHOICES), default="ct")
@click.option("--mode", type=click.Choice(["se", "pe"]), default="pe")
def simulate(output, genome_length, num_reads, read_length,
             methylation_rate, conversion, mode):
    """Generate simulated bisulfite/converted sequencing data."""
    from .simulation import Simulator
    sim = Simulator(
        genome_length=genome_length,
        num_reads=num_reads,
        read_length=read_length,
        methylation_rate=methylation_rate,
        conversion=conversion,
        mode=mode,
    )
    sim.run(output)
    click.echo(f"Simulated data written to {output}")


@main.command()
@click.argument("result_dir")
def report(result_dir):
    """Generate report from completed benchmark results."""
    from .reporting.report import generate_report
    generate_report(Path(result_dir))
    click.echo(f"Report generated in {result_dir}")


@main.command()
def list_aligners():
    """List all supported aligners with descriptions."""
    aligners = [
        ("bwameth", "Bwa-meth: BWA-based aligner for bisulfite reads"),
        ("bsbolt", "BSBolt: Bisulfite sequence aligner"),
        ("bsmap", "BSMAP: Bisulfite mapping program"),
        ("walt", "Walt: Fast bisulfite read aligner"),
        ("abismal", "Abismal: Bisulfite read aligner"),
        ("batmeth2", "Batmeth2: Bisulfite aligner"),
        ("basal", "Basal: Bisulfite alignment algorithm"),
        ("hisat3n", "HISAT-3n: Hierarchical indexing for spliced alignment"),
        ("hisat3n_repeat", "HISAT-3n (repeat): Repeat-aware alignment"),
        ("bismark_bwt2_e2e", "Bismark-bwt2-e2e: End-to-end bowtie2 mode"),
        ("bismark_his2", "Bismark-his2: HISAT2 mode"),
        ("bsseeker2_bwt", "BSseeker2-bwt: BWA-backtrack mode"),
        ("bsseeker2_soap2", "BSseeker2-soap2: SOAP2 mode"),
        ("bsseeker2_bwt2_e2e", "BSseeker2-bwt2-e2e: End-to-end bowtie2 mode"),
        ("bsseeker2_bwt2_local", "BSseeker2-bwt2-local: Local bowtie2 mode"),
        ("mock_*", "Built-in mock aligners for demo/testing"),
        ("conversion_aware", "Built-in: 2-base seed + 2-bit bitmask + banded SW (all 12 conversions)"),
    ]
    click.echo(f"{'Aligner':<25} Description")
    click.echo("-" * 80)
    for name, desc in aligners:
        click.echo(f"{name:<25} {desc}")


@main.command()
@click.option("-o", "--output", default="./demo_results")
@click.option("--conversion", type=click.Choice(CONV_CHOICES), default="ct")
def demo(output, conversion):
    """Run a quick demo without any external dependencies."""
    click.echo("Starting BCTool demo (no external tools required)...")
    cfg = Config()
    cfg.data["output_dir"] = output
    cfg.data["conversion"] = conversion
    cfg.data["simulation"]["num_reads"] = 10000
    cfg.data["simulation"]["genome_length"] = 200000
    cfg.data["extraction"]["enabled"] = False
    cfg.data["aligners"]["enabled"] = []
    pipe = Pipeline(cfg)
    pipe.run()


@main.command()
@click.option("-o", "--output", default="./demo_results")
@click.option("--conversion", type=click.Choice(CONV_CHOICES), default="ct")
def demo(output, conversion):
    """Run a quick demo without any external dependencies."""
    click.echo("Starting BCTool demo (no external tools required)...")
    cfg = Config()
    cfg.data["output_dir"] = output
    cfg.data["conversion"] = conversion
    cfg.data["simulation"]["num_reads"] = 10000
    cfg.data["simulation"]["genome_length"] = 200000
    cfg.data["extraction"]["enabled"] = False
    cfg.data["aligners"]["enabled"] = []
    pipe = Pipeline(cfg)
    pipe.run()


@main.command()
@click.option("-o", "--output", default="./batch_demo_results")
@click.option("--samples", type=int, default=4,
              help="Number of simulated samples")
@click.option("--workers", type=int, default=2,
              help="Max parallel workers")
def batch_demo(output, samples, workers):
    """Run a multi-sample demo with simulated data (no external tools)."""
    from .multisample import MultiSampleRunner, generate_batch_config
    import tempfile

    click.echo("Generating batch demo config...")
    tmp_config = Path(tempfile.mkdtemp()) / "batch_demo.yaml"

    demo_samples = []
    conditions = [
        ("ct", "high_meth", 0.8),
        ("ct", "low_meth", 0.3),
        ("ag", "high_meth", 0.8),
        ("ag", "low_meth", 0.3),
    ]
    for i in range(min(samples, len(conditions))):
        conv, group, meth_rate = conditions[i]
        name = f"sample_{conv}_{group}"
        demo_samples.append({
            "name": name,
            "conversion": conv,
            "mode": "pe",
            "group": group,
            "sim_reads": 15000,
        })

    generate_batch_config(tmp_config, demo_samples)

    runner = MultiSampleRunner(str(tmp_config))
    runner.global_config["use_mock"] = True
    runner.run_all(output, max_workers=workers, use_mock=True)
    click.echo(f"\nBatch demo complete! Results in: {output}")


@main.command()
@click.argument("sam")
@click.option("--conversion", type=click.Choice(CONV_CHOICES), default="ct",
              help="Conversion type")
@click.option("-o", "--output", default="./perbase_methylation.csv",
              help="Output CSV path")
@click.option("--min-qual", type=int, default=0,
              help="Minimum base quality score (Phred)")
def perbase(sam, conversion, output, min_qual):
    """Extract per-position conversion data from SAM to CSV.

    Output columns: ref, pos, strand, convertedBaseQualities,
    convertedBaseCount, unconvertedBaseQualities, unconvertedBaseCount, 未转化率
    """
    from .methylation.perbase_caller import PerBaseCaller
    from .methylation.converter import ConversionType

    conv = ConversionType(conversion)
    caller = PerBaseCaller(conv, min_qual=min_qual)
    click.echo(f"Parsing {sam}...")
    sites = caller.parse_sam(sam)
    caller.to_csv(sites, output)
    click.echo(f"Per-base methylation CSV written to {output}")


@main.command()
@click.option("--port", default=8501, help="Web server port")
@click.option("--host", default="0.0.0.0", help="Web server host")
def web(port, host):
    """Launch the interactive web interface (Streamlit)."""
    from .web import run_web
    click.echo(f"Starting BCTool web interface at http://{host}:{port}")
    click.echo("Press Ctrl+C to stop")
    run_web()


if __name__ == "__main__":
    main()
