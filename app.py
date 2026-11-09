"""
BioNexus - Multi-Omics Bioinformatics Platform
==============================================
Phase: Version 0.1 - Genomics Foundation
Architecture: Modular Industry-Grade Bioinformatics Platform
"""

import os
from pathlib import Path
from typing import Optional
import pandas as pd
import streamlit as st

from src.fasta_parser import FastaParseError, parse_fasta
from src.sequence_analysis import DatasetSummary, analyze_fasta_records
from src.visualization import (
    plot_gc_content_distribution,
    plot_gc_vs_at_comparison,
    plot_gc_vs_length_scatter,
    plot_overall_nucleotide_composition,
    plot_per_sequence_composition,
    plot_sequence_lengths,
)

# ------------------------------------------------------------------------------
# 1. Page Configuration & Custom Styling
# ------------------------------------------------------------------------------
st.set_page_config(
    page_title="BioNexus | Multi-Omics Platform",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

CUSTOM_CSS = """
<style>
    /* Global BioNexus Theme Adjustments */
    .stApp {
        background-color: #0b0f19;
        color: #e2e8f0;
        font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
    }
    
    /* Brand Header Card */
    .brand-container {
        background: linear-gradient(135deg, rgba(30, 41, 59, 0.7) 0%, rgba(15, 23, 42, 0.9) 100%);
        border: 1px solid rgba(56, 189, 248, 0.2);
        border-radius: 12px;
        padding: 24px 28px;
        margin-bottom: 24px;
        box-shadow: 0 4px 20px -2px rgba(0, 0, 0, 0.5);
    }
    
    .brand-title {
        font-size: 2.2rem;
        font-weight: 800;
        letter-spacing: -0.025em;
        background: linear-gradient(90deg, #38bdf8 0%, #818cf8 50%, #34d399 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        margin: 0;
        display: flex;
        align-items: center;
        gap: 12px;
    }
    
    .brand-subtitle {
        font-size: 0.98rem;
        color: #94a3b8;
        margin-top: 6px;
        font-weight: 400;
    }
    
    .badge-phase {
        display: inline-block;
        background-color: rgba(56, 189, 248, 0.15);
        color: #38bdf8;
        border: 1px solid rgba(56, 189, 248, 0.35);
        padding: 2px 10px;
        border-radius: 9999px;
        font-size: 0.75rem;
        font-weight: 600;
        letter-spacing: 0.05em;
        text-transform: uppercase;
        vertical-align: middle;
        margin-left: 10px;
    }
    
    /* Scientific Metric Cards */
    .metric-card {
        background: rgba(30, 41, 59, 0.5);
        border: 1px solid rgba(148, 163, 184, 0.12);
        border-radius: 10px;
        padding: 16px 20px;
        transition: transform 0.2s ease, border-color 0.2s ease;
    }
    .metric-card:hover {
        border-color: rgba(56, 189, 248, 0.4);
    }
    .metric-card-valid {
        background: rgba(16, 185, 129, 0.08);
        border: 1px solid rgba(16, 185, 129, 0.3);
        border-radius: 10px;
        padding: 16px 20px;
    }
    .metric-card-valid:hover {
        border-color: rgba(52, 211, 153, 0.6);
    }
    .metric-label {
        font-size: 0.8rem;
        text-transform: uppercase;
        letter-spacing: 0.05em;
        color: #94a3b8;
        font-weight: 600;
    }
    .metric-label-valid {
        font-size: 0.8rem;
        text-transform: uppercase;
        letter-spacing: 0.05em;
        color: #6ee7b7;
        font-weight: 600;
    }
    .metric-value {
        font-size: 1.8rem;
        font-weight: 700;
        color: #f8fafc;
        margin-top: 4px;
    }
    .metric-footnote {
        font-size: 0.75rem;
        color: #64748b;
        margin-top: 2px;
    }

    /* Section Cards */
    .section-box {
        background: rgba(15, 23, 42, 0.6);
        border: 1px solid rgba(148, 163, 184, 0.1);
        border-radius: 12px;
        padding: 20px 24px;
        margin-bottom: 24px;
    }

    /* Module Roadmap Card */
    .module-card {
        background: rgba(30, 41, 59, 0.35);
        border: 1px solid rgba(148, 163, 184, 0.1);
        border-radius: 10px;
        padding: 18px 20px;
        margin-bottom: 16px;
    }
    .module-badge-ready {
        color: #34d399;
        background: rgba(52, 211, 153, 0.15);
        border: 1px solid rgba(52, 211, 153, 0.3);
        padding: 2px 8px;
        border-radius: 6px;
        font-size: 0.72rem;
        font-weight: 600;
    }
    .module-badge-upcoming {
        color: #fbbf24;
        background: rgba(251, 191, 36, 0.12);
        border: 1px solid rgba(251, 191, 36, 0.25);
        padding: 2px 8px;
        border-radius: 6px;
        font-size: 0.72rem;
        font-weight: 600;
    }
</style>
"""
st.markdown(CUSTOM_CSS, unsafe_allow_html=True)


# ------------------------------------------------------------------------------
# 2. Sidebar Navigation & Global Controls
# ------------------------------------------------------------------------------
def render_sidebar():
    with st.sidebar:
        st.markdown(
            """
            <div style="display: flex; align-items: center; gap: 10px; margin-bottom: 8px;">
                <span style="font-size: 2rem;">🧬</span>
                <span style="font-size: 1.4rem; font-weight: 800; color: #f8fafc; letter-spacing: -0.02em;">BioNexus</span>
            </div>
            <div style="font-size: 0.78rem; color: #94a3b8; margin-bottom: 20px;">
                Industry Multi-Omics Platform • v0.1
            </div>
            """,
            unsafe_allow_html=True
        )

        st.subheader("Platform Navigation")
        module_choice = st.radio(
            "Select Module",
            options=[
                "Genomics (FASTA Analysis)",
                "Platform Overview",
                "Transcriptomics (Coming Soon)",
                "Proteomics (Coming Soon)",
                "Metagenomics (Coming Soon)",
                "Epigenomics (Coming Soon)",
                "Multi-Omics Integration (Coming Soon)"
            ],
            index=0,
            label_visibility="collapsed"
        )

        st.divider()

        st.subheader("Quick Dataset Presets")
        st.caption("Load curated sample sequences directly into the Genomics engine:")
        col1, col2 = st.columns(2)
        with col1:
            if st.button("E. coli Genes", use_container_width=True):
                st.session_state["sample_choice"] = "ecoli"
        with col2:
            if st.button("Edge Cases", use_container_width=True):
                st.session_state["sample_choice"] = "edge"

        st.divider()
        st.markdown(
            """
            <div style="font-size: 0.72rem; color: #64748b; line-height: 1.4;">
                <b>Local Processing Engine</b><br>
                • Zero Cloud Dependency<br>
                • Offline & Private<br>
                • Windows 64-bit Verified
            </div>
            """,
            unsafe_allow_html=True
        )

    return module_choice


# ------------------------------------------------------------------------------
# 3. View: Platform Overview
# ------------------------------------------------------------------------------
def render_platform_overview():
    st.markdown(
        """
        <div class="brand-container">
            <h1 class="brand-title">
                BioNexus Platform Overview
                <span class="badge-phase">v0.1 Roadmap</span>
            </h1>
            <p class="brand-subtitle">
                An industry-oriented, modular multi-omics bioinformatics platform engineered for reproducible 
                genomic, transcriptomic, proteomic, and systems biology research.
            </p>
        </div>
        """,
        unsafe_allow_html=True
    )

    st.markdown("### 🗺️ Multi-Omics Architecture Roadmap")

    modules_info = [
        {
            "name": "Genomics Foundation",
            "phase": "v0.1 — ACTIVE",
            "status": "ready",
            "icon": "🧬",
            "desc": "High-performance FASTA parser, nucleotide stoichiometry, IUPAC degenerate base tracking, strict GC% profiling, and quality control metrics."
        },
        {
            "name": "Transcriptomics Engine",
            "phase": "v0.2 — IN DESIGN",
            "status": "upcoming",
            "icon": "🧪",
            "desc": "RNA-Seq quantification, differential gene expression (DESeq2/edgeR equivalents), volcano plots, TPM/FPKM normalization, and pathway enrichment."
        },
        {
            "name": "Proteomics Suite",
            "phase": "v0.3 — IN DESIGN",
            "status": "upcoming",
            "icon": "🔬",
            "desc": "Mass spectrometry peptide identification, protein-protein interaction networks, post-translational modification analysis, and isobaric labeling quant."
        },
        {
            "name": "Metagenomics & Microbiome",
            "phase": "v0.4 — PLANNED",
            "status": "upcoming",
            "icon": "🦠",
            "desc": "16S rRNA taxonomy classification, shotgun metagenome binning, taxonomic abundance heatmaps, and alpha/beta diversity computation."
        },
        {
            "name": "Epigenomics Profiler",
            "phase": "v0.5 — PLANNED",
            "status": "upcoming",
            "icon": "⚛️",
            "desc": "DNA methylation calling, histone modification peak exploration, chromatin accessibility (ATAC-Seq), and epigenetic biomarker mapping."
        },
        {
            "name": "Cross-Omics Integrative Analysis",
            "phase": "v1.0 — LONG TERM",
            "status": "upcoming",
            "icon": "🌐",
            "desc": "Multi-omics factor analysis (MOFA), joint dimensionality reduction, cross-layer biological interpretation, and automated PDF research reporting."
        }
    ]

    for m in modules_info:
        badge_class = "module-badge-ready" if m["status"] == "ready" else "module-badge-upcoming"
        st.markdown(
            f"""
            <div class="module-card">
                <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 6px;">
                    <div style="font-size: 1.1rem; font-weight: 700; color: #f1f5f9;">
                        {m['icon']} {m['name']}
                    </div>
                    <span class="{badge_class}">{m['phase']}</span>
                </div>
                <div style="font-size: 0.88rem; color: #94a3b8; line-height: 1.4;">
                    {m['desc']}
                </div>
            </div>
            """,
            unsafe_allow_html=True
        )


# ------------------------------------------------------------------------------
# 4. View: Coming Soon Placeholder Modules
# ------------------------------------------------------------------------------
def render_coming_soon_module(module_name: str, icon: str, description: str, planned_version: str):
    st.markdown(
        f"""
        <div class="brand-container">
            <h1 class="brand-title">
                {icon} {module_name}
                <span class="badge-phase">Roadmap: {planned_version}</span>
            </h1>
            <p class="brand-subtitle">
                This module is part of the BioNexus modular roadmap and is currently under active architectural design.
            </p>
        </div>
        """,
        unsafe_allow_html=True
    )

    st.info(f"💡 **Upcoming Capability:** {description}")

    col1, col2 = st.columns([2, 1])
    with col1:
        st.markdown(
            f"""
            #### Planned Core Features for {module_name}:
            - **End-to-End Pipeline:** Automated ingestion, quality control, normalization, and statistical modeling.
            - **Local Scientific Compute:** High-efficiency vectorized processing with Pandas and NumPy.
            - **Interactive Visualizations:** Publication-ready Plotly charts, heatmaps, volcano plots, and correlation matrices.
            - **Exportable Outputs:** Formatted scientific CSV summaries and automated PDF research dossiers.
            """
        )
    with col2:
        st.markdown(
            """
            #### BioNexus Development Phasing
            1. **v0.1:** Genomics Foundation *(Current)*
            2. **v0.2:** Transcriptomics Engine
            3. **v0.3:** Proteomics Suite
            4. **v0.4:** Metagenomics & Microbiome
            5. **v0.5:** Epigenomics Profiler
            6. **v1.0:** Multi-Omics Integration & Reports
            """
        )


# ------------------------------------------------------------------------------
# 5. View: Genomics Module (DNA FASTA Analysis & QC)
# ------------------------------------------------------------------------------
def render_genomics_module():
    st.markdown(
        """
        <div class="brand-container">
            <h1 class="brand-title">
                🧬 Genomics Analysis Engine
                <span class="badge-phase">v0.1 Active</span>
            </h1>
            <p class="brand-subtitle">
                High-throughput DNA FASTA sequence parsing, stoichiometry breakdown, strict GC% calculation, 
                IUPAC degenerate base audit, and quality control profiling.
            </p>
        </div>
        """,
        unsafe_allow_html=True
    )

    # --- Dataset Upload & Preset Selection Area ---
    st.subheader("1. FASTA Dataset Ingestion")
    st.caption("Upload raw sequence records formatted as `.fa`, `.fasta`, `.fna`, or `.txt`:")

    col_upload, col_preset = st.columns([3, 1])

    with col_upload:
        uploaded_file = st.file_uploader(
            "Choose a FASTA file",
            type=["fasta", "fa", "fna", "txt"],
            help="Upload single or multi-sequence FASTA files with standard headers (>id description)",
            label_visibility="collapsed"
        )

    with col_preset:
        st.markdown("<div style='height: 4px;'></div>", unsafe_allow_html=True)
        preset_action = st.selectbox(
            "Quick Load Presets",
            options=["Select Preset...", "E. coli Reference Genes", "Degenerate & Edge Cases"],
            index=0
        )

    # Determine input data source
    records = None
    source_filename = None

    sample_choice = st.session_state.get("sample_choice", None)
    if preset_action == "E. coli Reference Genes" or sample_choice == "ecoli":
        sample_path = Path("data/uploads/sample_ecoli_genes.fasta")
        if sample_path.exists():
            source_filename = "sample_ecoli_genes.fasta (Preset)"
            with open(sample_path, "r", encoding="utf-8") as f:
                records = parse_fasta(f)
            if "sample_choice" in st.session_state:
                del st.session_state["sample_choice"]

    elif preset_action == "Degenerate & Edge Cases" or sample_choice == "edge":
        sample_path = Path("data/uploads/sample_ambiguous_and_edge.fasta")
        if sample_path.exists():
            source_filename = "sample_ambiguous_and_edge.fasta (Preset)"
            with open(sample_path, "r", encoding="utf-8") as f:
                records = parse_fasta(f)
            if "sample_choice" in st.session_state:
                del st.session_state["sample_choice"]

    elif uploaded_file is not None:
        source_filename = uploaded_file.name
        try:
            records = parse_fasta(uploaded_file)
        except FastaParseError as e:
            st.error(f"⚠️ FASTA Parsing Error: {str(e)}")
            st.info("Ensure the file starts with a valid header line beginning with `>` and contains valid nucleotide sequences.")
            return
        except Exception as ex:
            st.error(f"Unexpected error while reading file: {str(ex)}")
            return

    if records is None:
        st.info("👆 Upload a FASTA file or select a preset dataset above to initiate analysis.")
        return

    # --- Run Analysis ---
    with st.spinner(f"Analyzing {len(records)} sequence records..."):
        summary: DatasetSummary = analyze_fasta_records(records, filename=source_filename)
        df_results = summary.to_dataframe()
        qc_df = summary.generate_qc_report_dataframe()

    # --- Section: High-Level KPI Metrics (All Records Overview) ---
    st.markdown("### 2. Dataset Executive Summary")

    st.markdown("#### 🌐 Overall Dataset Overview (All Records)")
    m_col1, m_col2, m_col3, m_col4, m_col5, m_col6 = st.columns(6)

    with m_col1:
        st.markdown(
            f"""
            <div class="metric-card">
                <div class="metric-label">Filename</div>
                <div style="font-size: 1.05rem; font-weight: 700; color: #38bdf8; word-break: break-all; margin-top: 6px;">
                    {summary.filename}
                </div>
                <div class="metric-footnote">Source identifier</div>
            </div>
            """,
            unsafe_allow_html=True
        )

    with m_col2:
        st.markdown(
            f"""
            <div class="metric-card">
                <div class="metric-label">Total Records</div>
                <div class="metric-value">{summary.total_sequences:,}</div>
                <div class="metric-footnote">{summary.valid_sequences} valid • {summary.invalid_sequences} invalid</div>
            </div>
            """,
            unsafe_allow_html=True
        )

    with m_col3:
        st.markdown(
            f"""
            <div class="metric-card">
                <div class="metric-label">Total Length</div>
                <div class="metric-value">{summary.total_length:,} <span style="font-size: 1rem; color:#94a3b8">bp</span></div>
                <div class="metric-footnote">Range: {summary.min_length:,} – {summary.max_length:,} bp</div>
            </div>
            """,
            unsafe_allow_html=True
        )

    with m_col4:
        st.markdown(
            f"""
            <div class="metric-card">
                <div class="metric-label">All-Records Avg Length</div>
                <div class="metric-value">{summary.all_avg_length:,.1f} <span style="font-size: 1rem; color:#94a3b8">bp</span></div>
                <div class="metric-footnote">Median: {summary.median_length:,.1f} bp (incl. all)</div>
            </div>
            """,
            unsafe_allow_html=True
        )

    with m_col5:
        st.markdown(
            f"""
            <div class="metric-card">
                <div class="metric-label">Overall GC%</div>
                <div class="metric-value" style="color: #38bdf8;">{summary.overall_gc_percentage:.2f}%</div>
                <div class="metric-footnote">Excl. ambiguous & invalid</div>
            </div>
            """,
            unsafe_allow_html=True
        )

    with m_col6:
        st.markdown(
            f"""
            <div class="metric-card">
                <div class="metric-label">Overall AT%</div>
                <div class="metric-value" style="color: #f59e0b;">{summary.overall_at_percentage:.2f}%</div>
                <div class="metric-footnote">Complement ratio</div>
            </div>
            """,
            unsafe_allow_html=True
        )

    st.markdown("<div style='height: 8px;'></div>", unsafe_allow_html=True)

    # --- Section: Separate Valid Records Summary (Requirements 1 & 2) ---
    st.markdown("#### ✅ Valid Records Quality Baseline")
    valid_pct = (summary.valid_sequences / summary.total_sequences * 100) if summary.total_sequences > 0 else 0.0

    v_col1, v_col2, v_col3, v_col4, v_col5 = st.columns(5)

    with v_col1:
        st.markdown(
            f"""
            <div class="metric-card-valid">
                <div class="metric-label-valid">Valid Sequences</div>
                <div class="metric-value" style="color: #34d399;">{summary.valid_sequences:,} <span style="font-size: 1rem; color:#94a3b8">/ {summary.total_sequences:,}</span></div>
                <div class="metric-footnote">{valid_pct:.1f}% of uploaded dataset</div>
            </div>
            """,
            unsafe_allow_html=True
        )

    with v_col2:
        st.markdown(
            f"""
            <div class="metric-card-valid">
                <div class="metric-label-valid">Valid Total Length</div>
                <div class="metric-value" style="color: #34d399;">{summary.valid_total_length:,} <span style="font-size: 1rem; color:#94a3b8">bp</span></div>
                <div class="metric-footnote">Range: {summary.valid_min_length:,} – {summary.valid_max_length:,} bp</div>
            </div>
            """,
            unsafe_allow_html=True
        )

    with v_col3:
        st.markdown(
            f"""
            <div class="metric-card-valid">
                <div class="metric-label-valid">Valid-Records Avg Length</div>
                <div class="metric-value" style="color: #34d399;">{summary.valid_avg_length:,.1f} <span style="font-size: 1rem; color:#94a3b8">bp</span></div>
                <div class="metric-footnote">Median: {summary.valid_median_length:,.1f} bp (valid only)</div>
            </div>
            """,
            unsafe_allow_html=True
        )

    with v_col4:
        st.markdown(
            f"""
            <div class="metric-card-valid">
                <div class="metric-label-valid">Valid GC%</div>
                <div class="metric-value" style="color: #38bdf8;">{summary.valid_gc_percentage:.2f}%</div>
                <div class="metric-footnote">Valid canonical bases only</div>
            </div>
            """,
            unsafe_allow_html=True
        )

    with v_col5:
        st.markdown(
            f"""
            <div class="metric-card-valid">
                <div class="metric-label-valid">Valid AT%</div>
                <div class="metric-value" style="color: #f59e0b;">{summary.valid_at_percentage:.2f}%</div>
                <div class="metric-footnote">Valid canonical bases only</div>
            </div>
            """,
            unsafe_allow_html=True
        )

    st.markdown("<div style='height: 12px;'></div>", unsafe_allow_html=True)

    # --- Section: Stoichiometry & GC% Mathematical Specification (Requirements 7, 8, 9, 10) ---
    with st.expander("🔬 Nucleotide Stoichiometry & GC% Formula Specification", expanded=True):
        doc_col1, doc_col2 = st.columns([1, 1])

        with doc_col1:
            st.markdown("##### 🧬 Canonical vs Ambiguous Nucleotide Counts")
            st.markdown(
                f"""
                - **Adenine (A):** `{summary.total_A:,}` bp
                - **Thymine (T):** `{summary.total_T:,}` bp
                - **Guanine (G):** `{summary.total_G:,}` bp
                - **Cytosine (C):** `{summary.total_C:,}` bp
                - **Total Canonical Bases (A + T + G + C):** `{summary.total_canonical:,}` bp
                - **Total IUPAC Ambiguous Bases (N, R, Y, etc.):** `{summary.total_ambiguous:,}` bp
                - **Total Non-IUPAC Invalid Characters:** `{summary.total_invalid:,}` characters
                """
            )

        with doc_col2:
            st.markdown("##### 📐 Mathematical Definition of GC% Denominator")
            st.latex(r"\text{GC\%} = \left(\frac{G + C}{A + T + G + C}\right) \times 100")
            st.markdown(
                f"""
                - **Canonical Denominator Value:** `{summary.total_canonical:,}` bp
                - **G + C Numerator Value:** `{(summary.total_G + summary.total_C):,}` bp
                - **Strict Exclusion Rule:** Ambiguous IUPAC bases (*N, R, Y, S, W, K, M, B, D, H, V*) 
                  and non-IUPAC invalid characters are **strictly excluded** from the GC% denominator 
                  to preserve exact stoichiometry.
                """
            )

    st.markdown("<div style='height: 16px;'></div>", unsafe_allow_html=True)

    # --- Section: Quality Control & Anomaly Inspection (Requirements 4, 5, 6, 11) ---
    st.markdown("### 3. Quality Control (QC) & Anomaly Inspection")

    has_anomalies = (summary.invalid_sequences > 0 or summary.empty_sequences > 0 or summary.total_ambiguous > 0)

    if not has_anomalies:
        st.success(f"✅ **Certified Quality Control Passed:** All {summary.total_sequences} sequences are canonical, non-empty, and free of illegal characters.")
    else:
        if summary.invalid_sequences > 0:
            st.error(
                f"🚨 **Invalid Sequences Detected ({summary.invalid_sequences} record(s)):** Sequence IDs: "
                f"`{', '.join(summary.invalid_sequence_ids)}`"
            )

        if summary.empty_sequences > 0:
            st.warning(
                f"⚠️ **Empty Sequences Detected ({summary.empty_sequences} record(s)):** Sequence IDs: "
                f"`{', '.join(summary.empty_sequence_ids)}`"
            )

        if summary.total_ambiguous > 0:
            st.info(
                f"ℹ️ **Ambiguous IUPAC Sequences ({len(summary.ambiguous_sequence_ids)} record(s)):** Sequence IDs: "
                f"`{', '.join(summary.ambiguous_sequence_ids)}`"
            )

    # Display QC Report Table
    st.markdown("##### 📋 Quality Control Audit Table")
    st.dataframe(
        qc_df,
        use_container_width=True,
        hide_index=True,
        column_config={
            "sequence_id": st.column_config.TextColumn("Sequence ID", width="medium"),
            "qc_status": st.column_config.TextColumn("QC Status", width="small"),
            "error_types": st.column_config.TextColumn("Error Category", width="medium"),
            "issue_details": st.column_config.TextColumn("Issue Details", width="large"),
            "invalid_character_count": st.column_config.NumberColumn("Invalid Chars", format="%d"),
            "invalid_characters": st.column_config.TextColumn("Illegal Chars"),
            "invalid_positions": st.column_config.TextColumn("Positions (1-based)", width="large"),
            "ambiguous_count": st.column_config.NumberColumn("Ambiguous Bases", format="%d"),
            "recommended_action": st.column_config.TextColumn("Recommended Action", width="large")
        }
    )

    # Download button for QC Error Report CSV (Requirement 11)
    qc_csv_data = qc_df.to_csv(index=False).encode("utf-8")
    qc_export_filename = f"bionexus_{Path(source_filename).stem}_qc_error_report.csv"

    qc_down_col1, qc_down_col2 = st.columns([1, 3])
    with qc_down_col1:
        st.download_button(
            label="📥 Download QC Error Report (CSV)",
            data=qc_csv_data,
            file_name=qc_export_filename,
            mime="text/csv",
            use_container_width=True
        )

    with qc_down_col2:
        # Save snapshot of QC error report to outputs directory
        output_dir = Path("outputs")
        output_dir.mkdir(exist_ok=True)
        local_qc_path = output_dir / qc_export_filename
        try:
            with open(local_qc_path, "wb") as f_qc:
                f_qc.write(qc_csv_data)
            st.caption(f"QC snapshot automatically synchronized to: `{local_qc_path}`")
        except Exception:
            pass

    st.markdown("<div style='height: 20px;'></div>", unsafe_allow_html=True)

    # --- Section: Interactive Visualizations ---
    st.markdown("### 4. Interactive Genomics Visualizations")

    tab1, tab2, tab3, tab4 = st.tabs([
        "📊 Nucleotide Composition",
        "⚖️ GC Content Analysis",
        "📏 Sequence Lengths",
        "🎯 QC Scatter (Length vs GC%)"
    ])

    with tab1:
        col_c1, col_c2 = st.columns([1, 1])
        with col_c1:
            fig_overall = plot_overall_nucleotide_composition(summary)
            st.plotly_chart(fig_overall, use_container_width=True)
        with col_c2:
            fig_per_seq = plot_per_sequence_composition(df_results, max_records=25)
            st.plotly_chart(fig_per_seq, use_container_width=True)

    with tab2:
        col_gc1, col_gc2 = st.columns([2, 1])
        with col_gc1:
            fig_gc_dist = plot_gc_content_distribution(df_results, summary)
            st.plotly_chart(fig_gc_dist, use_container_width=True)
        with col_gc2:
            fig_gc_pie = plot_gc_vs_at_comparison(summary)
            st.plotly_chart(fig_gc_pie, use_container_width=True)

    with tab3:
        fig_lengths = plot_sequence_lengths(df_results, summary)
        st.plotly_chart(fig_lengths, use_container_width=True)

    with tab4:
        fig_scatter = plot_gc_vs_length_scatter(df_results)
        st.plotly_chart(fig_scatter, use_container_width=True)

    st.markdown("<div style='height: 16px;'></div>", unsafe_allow_html=True)

    # --- Section: Individual Sequence Statistics Table ---
    st.markdown("### 5. Individual Sequence Records & Verification")

    filter_col1, filter_col2, filter_col3 = st.columns([2, 1, 1])
    with filter_col1:
        search_query = st.text_input("Filter by Sequence ID or Description", placeholder="e.g. dnaA, b0001, rec1...")
    with filter_col2:
        status_filter = st.selectbox(
            "Filter by Status",
            options=["All Records", "Valid Only", "Contains Ambiguous Bases", "Invalid Characters", "Empty Sequences"]
        )
    with filter_col3:
        sort_by = st.selectbox(
            "Sort by",
            options=["Original Order", "Length (Descending)", "GC% (Descending)", "GC% (Ascending)"]
        )

    # Apply filters
    filtered_df = df_results.copy()

    if search_query:
        query = search_query.strip().lower()
        filtered_df = filtered_df[
            filtered_df["id"].str.lower().str.contains(query) |
            filtered_df["description"].str.lower().str.contains(query)
        ]

    if status_filter == "Valid Only":
        filtered_df = filtered_df[filtered_df["is_valid"] == True]
    elif status_filter == "Contains Ambiguous Bases":
        filtered_df = filtered_df[filtered_df["count_ambiguous"] > 0]
    elif status_filter == "Invalid Characters":
        filtered_df = filtered_df[filtered_df["count_invalid"] > 0]
    elif status_filter == "Empty Sequences":
        filtered_df = filtered_df[filtered_df["length"] == 0]

    # Apply sorting
    if sort_by == "Length (Descending)":
        filtered_df = filtered_df.sort_values(by="length", ascending=False)
    elif sort_by == "GC% (Descending)":
        filtered_df = filtered_df.sort_values(by="gc_percentage", ascending=False)
    elif sort_by == "GC% (Ascending)":
        filtered_df = filtered_df.sort_values(by="gc_percentage", ascending=True)

    st.dataframe(
        filtered_df,
        use_container_width=True,
        hide_index=True,
        column_config={
            "id": st.column_config.TextColumn("Sequence ID", width="medium"),
            "description": st.column_config.TextColumn("Description", width="large"),
            "length": st.column_config.NumberColumn("Length (bp)", format="%d"),
            "gc_percentage": st.column_config.NumberColumn("GC %", format="%.2f%%"),
            "at_percentage": st.column_config.NumberColumn("AT %", format="%.2f%%"),
            "count_A": st.column_config.NumberColumn("A", format="%d"),
            "count_T": st.column_config.NumberColumn("T", format="%d"),
            "count_G": st.column_config.NumberColumn("G", format="%d"),
            "count_C": st.column_config.NumberColumn("C", format="%d"),
            "count_canonical": st.column_config.NumberColumn("Canonical (bp)", format="%d"),
            "count_ambiguous": st.column_config.NumberColumn("Ambiguous", format="%d"),
            "ambiguous_breakdown": st.column_config.TextColumn("Ambiguous Breakdown"),
            "count_invalid": st.column_config.NumberColumn("Invalid", format="%d"),
            "invalid_characters": st.column_config.TextColumn("Illegal Chars"),
            "invalid_positions_summary": st.column_config.TextColumn("Illegal Char Positions (1-based)", width="large"),
            "is_valid": st.column_config.CheckboxColumn("Valid?"),
            "status": st.column_config.TextColumn("Status", width="medium")
        }
    )
    st.caption(f"Displaying {len(filtered_df)} of {len(df_results)} sequence records.")

    # --- Section: Export Options ---
    st.markdown("### 6. Export Full Sequence Analysis")
    down_col1, down_col2 = st.columns([1, 3])

    csv_data = df_results.to_csv(index=False).encode("utf-8")
    export_filename = f"bionexus_{Path(source_filename).stem}_analysis.csv"

    with down_col1:
        st.download_button(
            label="📥 Download Full Analysis CSV",
            data=csv_data,
            file_name=export_filename,
            mime="text/csv",
            use_container_width=True
        )

    with down_col2:
        output_dir = Path("outputs")
        output_dir.mkdir(exist_ok=True)
        local_saved_path = output_dir / export_filename
        try:
            with open(local_saved_path, "wb") as f_out:
                f_out.write(csv_data)
            st.success(f"💾 Full dataset snapshot saved to: `{local_saved_path}`")
        except Exception:
            pass


# ------------------------------------------------------------------------------
# 6. Main Application Dispatcher
# ------------------------------------------------------------------------------
def main():
    selected_module = render_sidebar()

    if selected_module == "Genomics (FASTA Analysis)":
        render_genomics_module()
    elif selected_module == "Platform Overview":
        render_platform_overview()
    elif selected_module == "Transcriptomics (Coming Soon)":
        render_coming_soon_module(
            module_name="Transcriptomics Engine",
            icon="🧪",
            description="High-throughput RNA-Seq quantification, differential expression modeling, volcano plots, and Gene Ontology (GO) pathway enrichment.",
            planned_version="v0.2"
        )
    elif selected_module == "Proteomics (Coming Soon)":
        render_coming_soon_module(
            module_name="Proteomics Suite",
            icon="🔬",
            description="Mass spectrometry spectra matching, protein identification, isobaric quantitation (TMT/iTRAQ), and post-translational modification profiling.",
            planned_version="v0.3"
        )
    elif selected_module == "Metagenomics (Coming Soon)":
        render_coming_soon_module(
            module_name="Metagenomics & Microbiome",
            icon="🦠",
            description="Taxonomic classification of 16S/18S rRNA amplicons, metagenomic binning, functional pathway abundance, and diversity index computation.",
            planned_version="v0.4"
        )
    elif selected_module == "Epigenomics (Coming Soon)":
        render_coming_soon_module(
            module_name="Epigenomics Profiler",
            icon="⚛️",
            description="Bisulfite sequencing DNA methylation analysis, ChIP-seq and ATAC-seq peak calling, and epigenetic regulatory element annotation.",
            planned_version="v0.5"
        )
    elif selected_module == "Multi-Omics Integration (Coming Soon)":
        render_coming_soon_module(
            module_name="Multi-Omics Integration",
            icon="🌐",
            description="Joint latent-factor modeling (MOFA), cross-omic regulatory networks, biomarker discovery, and automated executive research reporting.",
            planned_version="v1.0"
        )


if __name__ == "__main__":
    main()
