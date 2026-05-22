"""
BCTool Web Interface (Streamlit)

Provides an interactive web UI for:
- Uploading FASTQ files
- Configuring pipeline parameters
- Running analysis with progress feedback
- Viewing interactive results dashboard
- Comparing multiple samples
"""

import os
import sys
import json
import tempfile
from pathlib import Path

try:
    import streamlit as st
    import pandas as pd
    import plotly.express as px
    import plotly.graph_objects as go
    HAS_DEPS = True
except ImportError:
    HAS_DEPS = False


def run_web():
    """Launch the Streamlit web interface."""
    if not HAS_DEPS:
        print("Web interface requires: streamlit plotly pandas")
        print("Install: pip install streamlit plotly pandas")
        sys.exit(1)

    sys.argv = ["streamlit", "run", __file__, "--browser.serverAddress=0.0.0.0", "--server.port=8501"]
    from streamlit.web import cli
    cli.main()

# ============================================================================
# Streamlit App
# ============================================================================

if HAS_DEPS and __name__ == "__main__":

    st.set_page_config(
        page_title="BCTool",
        page_icon="",
        layout="wide",
        initial_sidebar_state="expanded",
    )

    # ─── Sidebar ──────────────────────────────────────────────────────────
    st.sidebar.title("BCTool")
    st.sidebar.caption("Multi-algorithm methylation alignment benchmark")

    page = st.sidebar.radio(
        "Navigation",
        ["Pipeline", "QC Analysis", "Results", "Multi-Sample", "About"],
    )

    # ─── Session State ────────────────────────────────────────────────────
    if "output_dir" not in st.session_state:
        st.session_state.output_dir = tempfile.mkdtemp()
    if "results" not in st.session_state:
        st.session_state.results = None
    if "qc_results" not in st.session_state:
        st.session_state.qc_results = None

    # ======================================================================
    # PAGE: Pipeline
    # ======================================================================
    if page == "Pipeline":
        st.title("Benchmark Pipeline")
        st.markdown("Configure and run the 14-aligner benchmarking pipeline.")

        col1, col2 = st.columns(2)

        with col1:
            st.subheader("Input Files")
            use_sim = st.checkbox("Use simulated data", value=True,
                                  help="Generate synthetic data with known ground truth")

            read1 = None
            read2 = None
            reference = None

            if not use_sim:
                read1 = st.file_uploader("Read 1 (FASTQ)", type=["fastq", "fq", "gz"])
                read2 = st.file_uploader("Read 2 (FASTQ, optional)", type=["fastq", "fq", "gz"])
                reference = st.file_uploader("Reference Genome (FASTA)", type=["fa", "fasta", "gz"])

                if read1:
                    tmp = Path(tempfile.mkdtemp())
                    (tmp / "input").mkdir()
                    with open(tmp / "input" / read1.name, "wb") as f:
                        f.write(read1.getbuffer())
                    read1 = str(tmp / "input" / read1.name)
                    if read2:
                        with open(tmp / "input" / read2.name, "wb") as f:
                            f.write(read2.getbuffer())
                        read2 = str(tmp / "input" / read2.name)

        with col2:
            st.subheader("Settings")
            conversion = st.selectbox(
                "Conversion",
                ["ct (C>T)", "tc (T>C)", "ag (A>G)", "ga (G>A)",
                 "ac (A>C)", "ca (C>A)", "gt (G>T)", "tg (T>G)",
                 "at (A>T)", "ta (T>A)", "cg (C>G)", "gc (G>C)"],
                index=0,
            )
            conversion = conversion.split()[0]

            threads = st.slider("Threads", 1, 32, 8)

        # Data amount selection
        with col1:
            st.subheader("Data Amount for Validation")
            use_pct = st.checkbox("Use percentage of data", value=False,
                                  help="Sample a percentage of reads instead of a fixed count")
            if use_pct:
                sample_pct = st.slider("Percentage of reads", 1, 100, 10,
                                       help="% of reads to use for validation")
                num_reads = 0
            else:
                sample_pct = None
                num_reads = st.number_input("Number of reads to use",
                                            min_value=1000, value=20000, step=1000)

        st.subheader("Preprocessing")
        col3, col4, col5, col6 = st.columns(4)
        with col3:
            do_qc = st.checkbox("Run QC", value=True)
        with col4:
            trim_adapters = st.checkbox("Trim adapters")
        with col5:
            do_umi = st.checkbox("UMI dedup")
        with col6:
            do_correct = st.checkbox("Error correct",
                                     help="Filter low-quality alignments and apply consensus correction")

        # Aligner selection
        st.subheader("Aligners")
        all_aligners = [
            "bwameth", "bsbolt", "bsmap", "walt", "abismal", "batmeth2",
            "basal",
            "hisat3n", "hisat3n_repeat", "bismark_bwt2_e2e", "bismark_his2",
            "bsseeker2_bwt", "bsseeker2_soap2", "bsseeker2_bwt2_e2e", "bsseeker2_bwt2_local",
        ]
        use_mock = st.checkbox("Use mock aligners (no external tools needed)", value=True)
        if use_mock:
            selected = []
            select_all = False
        else:
            col_a, col_b = st.columns([3, 1])
            with col_b:
                select_all = st.checkbox("Select all", value=False)
            with col_a:
                if select_all:
                    selected = st.multiselect("Select aligners", all_aligners,
                                              default=all_aligners)
                else:
                    selected = st.multiselect("Select aligners", all_aligners,
                                              default=all_aligners[:3])

        if st.button("Run Pipeline", type="primary", use_container_width=True):
            with st.spinner("Running benchmark pipeline... This may take a while."):
                from ..config import Config
                from ..pipeline import Pipeline

                cfg = Config()
                out_dir = Path(st.session_state.output_dir) / f"run_{len(list(Path(st.session_state.output_dir).glob('run_*')))}"
                out_dir.mkdir(parents=True, exist_ok=True)

                cfg.data["output_dir"] = str(out_dir)
                cfg.data["conversion"] = conversion
                cfg.data["threads"] = threads
                cfg.data["qc"]["enabled"] = do_qc

                if use_sim:
                    cfg.data["simulation"]["num_reads"] = num_reads or 20000
                    cfg.data["extraction"]["enabled"] = False
                else:
                    cfg.data["input"]["read1"] = read1
                    if read2:
                        cfg.data["input"]["read2"] = read2
                    cfg.data["input"]["reference"] = str(reference) if reference else None
                    if sample_pct:
                        cfg.data["extraction"]["sample_pct"] = sample_pct
                    else:
                        cfg.data["extraction"]["num_reads"] = num_reads

                if trim_adapters:
                    cfg.data["adapter"]["enabled"] = True
                if do_umi:
                    cfg.data["umi"]["enabled"] = True
                if do_correct:
                    cfg.data["correction"]["enabled"] = True
                if use_mock:
                    cfg.data["aligners"]["enabled"] = []
                else:
                    cfg.data["aligners"]["enabled"] = selected

                pipe = Pipeline(cfg)
                pipe.run()

                # Load results
                res_file = out_dir / "comparison_results.json"
                if res_file.exists():
                    with open(res_file) as f:
                        st.session_state.results = json.load(f)
                    st.success(f"Pipeline complete! Results in {out_dir}")
                else:
                    st.error("Pipeline finished but no results found.")

        # Show cached results
        if st.session_state.results:
            st.subheader("Latest Results")
            _display_results(st.session_state.results)

    # ======================================================================
    # PAGE: QC Analysis
    # ======================================================================
    elif page == "QC Analysis":
        st.title("Quality Control")
        st.markdown("Run FastQC-style quality analysis on FASTQ files.")

        col1, col2 = st.columns([2, 1])
        with col1:
            qc_file = st.file_uploader("Upload FASTQ", type=["fastq", "fq", "gz"], key="qc_upload")
        with col2:
            qc_n = st.number_input("Max reads", min_value=1000, value=50000, step=1000)

        if qc_file and st.button("Run QC", type="primary"):
            with st.spinner("Analyzing..."):
                from ..qc import QCAnalyzer, QCReport
                tmp = tempfile.mktemp(suffix=".fastq")
                with open(tmp, "wb") as f:
                    f.write(qc_file.getbuffer())

                analyzer = QCAnalyzer(max_reads=qc_n)
                analyzer.analyze(tmp)
                st.session_state.qc_results = analyzer.metrics

        if st.session_state.qc_results:
            _display_qc(st.session_state.qc_results)

    # ======================================================================
    # PAGE: Results
    # ======================================================================
    elif page == "Results":
        st.title("Results Browser")
        st.markdown("Browse and visualize completed benchmark results.")

        # Allow loading results from directory
        result_path = st.text_input("Path to results directory (optional)",
                                     placeholder="e.g., ./bctool_results")
        if result_path and Path(result_path).exists():
            res_file = Path(result_path) / "comparison_results.json"
            if res_file.exists():
                with open(res_file) as f:
                    st.session_state.results = json.load(f)
                st.success(f"Loaded results from {res_file}")
            else:
                st.warning(f"No comparison_results.json found in {result_path}")

        if st.session_state.results:
            _display_results(st.session_state.results)
            _display_interactive_charts(st.session_state.results)
        else:
            st.info("No results loaded. Run the pipeline or load existing results.")

    # ======================================================================
    # PAGE: Multi-Sample
    # ======================================================================
    elif page == "Multi-Sample":
        st.title("Multi-Sample Comparison")
        st.markdown("Compare results across multiple samples or conditions.")

        result_dirs = st.text_area(
            "Result directories (one per line)",
            placeholder="./sample1_results\n./sample2_results\n./sample3_results"
        )
        if result_dirs and st.button("Compare", type="primary"):
            dirs = [d.strip() for d in result_dirs.split("\n") if d.strip()]
            all_results = []
            for d in dirs:
                res_file = Path(d) / "comparison_results.json"
                if res_file.exists():
                    with open(res_file) as f:
                        data = json.load(f)
                    for r in data.get("results", []):
                        r["sample"] = Path(d).name
                        all_results.append(r)

            if all_results:
                df = pd.DataFrame(all_results)
                st.subheader("Multi-Sample Comparison")

                # Heatmap
                pivot = df.pivot_table(index="sample", columns="aligner",
                                       values="f1_score", aggfunc="mean")
                fig = px.imshow(pivot, text_auto=".3f", color_continuous_scale="YlOrRd",
                                title="F1 Score by Sample and Aligner")
                st.plotly_chart(fig, use_container_width=True)

                # Bar chart
                fig2 = px.bar(df, x="sample", y="f1_score", color="aligner",
                              barmode="group", title="F1 Score Comparison")
                st.plotly_chart(fig2, use_container_width=True)

                # Best per sample
                best = df.loc[df.groupby("sample")["f1_score"].idxmax()]
                st.subheader("Best Aligner Per Sample")
                st.dataframe(best[["sample", "aligner", "f1_score", "precision", "recall"]],
                            use_container_width=True)
            else:
                st.warning("No results found.")

    # ======================================================================
    # PAGE: About
    # ======================================================================
    elif page == "About":
        st.title("About BCTool")
        st.markdown("""
        ### Multi-Algorithm Methylation Alignment Benchmark

        **Supported Aligners:**
        1. Bwa-meth
        2. BSBolt
        3. BSMAP
        4. Walt
        5. Abismal
        6. Batmeth2
        7. Basal
        8. HISAT-3n
        9. HISAT-3n (repeat)
        10. Bismark-bwt2-e2e
        11. Bismark-his2
        12. BSseeker2-bwt
        13. BSseeker2-soap2
        14. BSseeker2-bwt2-e2e
        15. BSseeker2-bwt2-local

        **Features:**
        - 12 conversion types: C>T, T>C, A>G, G>A, A>C, C>A, G>T, T>G, A>T, T>A, C>G, G>C
        - Single-end and paired-end
        - Built-in data simulation with ground truth
        - QC analysis (FastQC-style)
        - Adapter trimming (auto-detect + cutadapt)
        - UMI deduplication
        - Multi-sample batch analysis
        - Docker containerization
        """)

        st.code("BCTool --help", language="bash")


# ============================================================================
# Display Helpers
# ============================================================================

def _display_results(results):
    summary = results.get("summary")
    if summary:
        cols = st.columns(4)
        cols[0].metric("Best Aligner", summary["aligner"])
        cols[1].metric("F1 Score", f"{summary['f1_score']:.4f}")
        cols[2].metric("Precision", f"{summary['precision']:.4f}")
        cols[3].metric("Recall", f"{summary['recall']:.4f}")

    res = results.get("results", [])
    if res:
        df = pd.DataFrame(res)
        if "rank" not in df.columns and not df.empty:
            df = df.sort_values("f1_score", ascending=False)
            df["rank"] = range(1, len(df) + 1)
        st.dataframe(df, use_container_width=True)


def _display_interactive_charts(results):
    res = results.get("results", [])
    if not res:
        return

    df = pd.DataFrame(res)

    st.subheader("Interactive Charts")

    tab1, tab2, tab3 = st.tabs(["F1 Score", "Metrics Comparison", "Radar Chart"])

    with tab1:
        fig = px.bar(df, x="aligner", y="f1_score",
                     color="f1_score", color_continuous_scale="Viridis",
                     title="F1 Score by Aligner",
                     text_auto=".4f")
        fig.update_layout(xaxis_tickangle=-45)
        st.plotly_chart(fig, use_container_width=True)

    with tab2:
        metrics_df = df.melt(id_vars=["aligner"],
                             value_vars=["precision", "recall", "accuracy", "f1_score"],
                             var_name="metric", value_name="score")
        fig = px.bar(metrics_df, x="aligner", y="score", color="metric",
                     barmode="group", title="All Metrics by Aligner")
        fig.update_layout(xaxis_tickangle=-45)
        st.plotly_chart(fig, use_container_width=True)

    with tab3:
        if len(df) > 0:
            # Radar chart for top aligners
            top = df.nlargest(5, "f1_score")
            fig = go.Figure()
            for _, row in top.iterrows():
                fig.add_trace(go.Scatterpolar(
                    r=[row["precision"], row["recall"], row["f1_score"], row["accuracy"],
                       row.get("level_correlation", 0) or row["precision"]],
                    theta=["Precision", "Recall", "F1", "Accuracy", "Level Corr"],
                    fill="toself",
                    name=row["aligner"],
                ))
            fig.update_layout(polar=dict(radialaxis=dict(visible=True, range=[0, 1])))
            st.plotly_chart(fig, use_container_width=True)


def _display_qc(metrics):
    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Total Reads", metrics.get("total_reads", 0))
    col2.metric("Total Bases", f"{metrics.get('total_bases', 0):,}")
    col3.metric("Mean Quality", f"{metrics.get('qual_summary', {}).get('mean', 0):.1f}")
    gc = metrics.get("gc_summary", {}).get("mean", 0)
    col4.metric("GC Content", f"{gc:.1f}%")

    # Per-base quality heatmap-style chart
    sq = metrics.get("per_base_qual_summary", {})
    if sq:
        st.subheader("Per-Base Quality")
        df_q = pd.DataFrame(sq).T
        df_q = df_q.reset_index()
        df_q.columns = ["position", "mean", "med", "q1", "q3", "lower", "upper", "min", "max"]
        fig = px.line(df_q, x="position", y=["mean", "med", "q1", "q3"],
                      title="Quality Scores Across Read Positions")
        st.plotly_chart(fig, use_container_width=True)

    # Base content
    content = metrics.get("per_base_content", {})
    if content:
        st.subheader("Per-Base Sequence Content")
        df_c = pd.DataFrame(content).T
        df_c = df_c.reset_index()
        cols = {"index": "position"}
        for b in "ACGTN":
            if b in df_c.columns:
                df_c[b] = df_c[b] / df_c[list("ACGTN")].sum(axis=1) * 100
                cols[b] = b
        df_c = df_c.rename(columns=cols)
        fig = px.line(df_c, x="position", y=["A", "C", "G", "T", "N"],
                      title="Base Composition Across Positions")
        st.plotly_chart(fig, use_container_width=True)

    # Overrepresented
    over = metrics.get("overrepresented", [])
    if over:
        st.subheader("Overrepresented Sequences")
        st.dataframe(pd.DataFrame(over[:10]), use_container_width=True)

    # Adapters
    adapters = metrics.get("detected_adapters", [])
    if adapters:
        st.warning(f"Adapter detected: {adapters[0]['adapter']} ({adapters[0]['pct']:.1f}%)")
