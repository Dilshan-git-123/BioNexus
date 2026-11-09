"""
BioNexus Interactive Visualization Module
=========================================
Generates publication-quality, interactive Plotly charts and tables
for genomics sequence inspection and quality control.

Visualizations:
- Nucleotide composition breakdown (Overall & Per-sequence)
- GC% distribution and GC vs AT comparative visualization
- Sequence length comparison with mean/median benchmarks
- GC% vs Sequence Length QC scatter plot
- Formatted summary statistics table
"""

from typing import Optional
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from src.sequence_analysis import DatasetSummary

# BioNexus Color Palette for Standard DNA Bases
BASE_COLORS = {
    "A": "#2ECC71",  # Emerald Green
    "T": "#E74C3C",  # Ruby Red
    "G": "#3498DB",  # Cobalt Blue
    "C": "#F39C12",  # Amber Orange
    "Ambiguous": "#9B59B6",  # Purple
    "Invalid": "#95A5A6"   # Slate Grey
}

THEME_LAYOUT = {
    "font": {"family": "Inter, Segoe UI, sans-serif", "size": 12, "color": "#E0E6ED"},
    "paper_bgcolor": "rgba(17, 24, 39, 0.8)",  # Dark slate background
    "plot_bgcolor": "rgba(17, 24, 39, 0.4)",
    "margin": {"l": 50, "r": 30, "t": 60, "b": 50},
    "legend": {"orientation": "h", "yanchor": "bottom", "y": 1.02, "xanchor": "right", "x": 1}
}


def plot_overall_nucleotide_composition(summary: DatasetSummary) -> go.Figure:
    """
    Creates an interactive bar chart of overall nucleotide composition across the dataset.
    """
    labels = ["Adenine (A)", "Thymine (T)", "Guanine (G)", "Cytosine (C)"]
    counts = [summary.total_A, summary.total_T, summary.total_G, summary.total_C]
    colors = [BASE_COLORS["A"], BASE_COLORS["T"], BASE_COLORS["G"], BASE_COLORS["C"]]

    if summary.total_ambiguous > 0:
        labels.append("Ambiguous (IUPAC)")
        counts.append(summary.total_ambiguous)
        colors.append(BASE_COLORS["Ambiguous"])

    if summary.total_invalid > 0:
        labels.append("Invalid Characters")
        counts.append(summary.total_invalid)
        colors.append(BASE_COLORS["Invalid"])

    total_bases = sum(counts)
    percentages = [(c / total_bases * 100) if total_bases > 0 else 0 for c in counts]

    fig = go.Figure(
        data=[
            go.Bar(
                x=labels,
                y=counts,
                marker=dict(
                    color=colors,
                    line=dict(color="#FFFFFF", width=1)
                ),
                text=[f"{count:,} ({pct:.1f}%)" for count, pct in zip(counts, percentages)],
                textposition="auto",
                hovertemplate="<b>%{x}</b><br>Count: %{y:,}<br>Proportion: %{text}<extra></extra>"
            )
        ]
    )

    fig.update_layout(
        title="<b>Overall Nucleotide Composition</b>",
        xaxis_title="Nucleotide Base",
        yaxis_title="Total Base Count",
        **THEME_LAYOUT
    )
    fig.update_xaxes(showgrid=False, gridcolor="#2D3748")
    fig.update_yaxes(showgrid=True, gridcolor="#2D3748")
    return fig


def plot_per_sequence_composition(df: pd.DataFrame, max_records: int = 50) -> go.Figure:
    """
    Generates a stacked horizontal bar chart showing nucleotide proportions per sequence.
    """
    if df.empty:
        return go.Figure()

    display_df = df.head(max_records).copy()
    display_df = display_df.iloc[::-1]  # Reverse for top-to-bottom reading

    fig = go.Figure()
    bases = [("count_A", "A", BASE_COLORS["A"]),
             ("count_T", "T", BASE_COLORS["T"]),
             ("count_G", "G", BASE_COLORS["G"]),
             ("count_C", "C", BASE_COLORS["C"])]

    for col, name, color in bases:
        fig.add_trace(
            go.Bar(
                y=display_df["id"],
                x=display_df[col],
                name=f"Base {name}",
                orientation="h",
                marker=dict(color=color)
            )
        )

    # Add ambiguous bases if present
    if "count_ambiguous" in display_df.columns and display_df["count_ambiguous"].sum() > 0:
        fig.add_trace(
            go.Bar(
                y=display_df["id"],
                x=display_df["count_ambiguous"],
                name="Ambiguous",
                orientation="h",
                marker=dict(color=BASE_COLORS["Ambiguous"])
            )
        )

    fig.update_layout(
        barmode="stack",
        title=f"<b>Per-Sequence Nucleotide Breakdown</b> (Top {len(display_df)} records)",
        xaxis_title="Nucleotide Count",
        yaxis_title="Sequence Identifier",
        height=max(380, 28 * len(display_df) + 120),
        **THEME_LAYOUT
    )
    fig.update_xaxes(showgrid=True, gridcolor="#2D3748")
    fig.update_yaxes(showgrid=False)
    return fig


def plot_gc_content_distribution(df: pd.DataFrame, summary: DatasetSummary) -> go.Figure:
    """
    Visualizes GC% across sequences using a histogram with a benchmark line for overall GC%.
    If dataset has <= 20 sequences, renders a clean bar chart with GC% values.
    """
    if df.empty:
        return go.Figure()

    if len(df) <= 25:
        # Direct bar chart for smaller record counts
        fig = px.bar(
            df,
            x="id",
            y="gc_percentage",
            color="gc_percentage",
            color_continuous_scale="Viridis",
            labels={"id": "Sequence ID", "gc_percentage": "GC Content (%)"},
            title="<b>GC Content per Sequence</b>"
        )
        fig.add_hline(
            y=summary.overall_gc_percentage,
            line_dash="dash",
            line_color="#E74C3C",
            annotation_text=f"Dataset Avg GC: {summary.overall_gc_percentage:.1f}%",
            annotation_position="top right"
        )
        fig.update_layout(
            yaxis_range=[0, 100],
            **THEME_LAYOUT
        )
    else:
        # Histogram distribution for larger datasets
        fig = px.histogram(
            df,
            x="gc_percentage",
            nbins=30,
            marginal="box",
            color_discrete_sequence=["#1ABC9C"],
            labels={"gc_percentage": "GC Content (%)"},
            title="<b>GC Content Distribution Across Dataset</b>"
        )
        fig.add_vline(
            x=summary.overall_gc_percentage,
            line_dash="dash",
            line_color="#E74C3C",
            annotation_text=f"Overall GC: {summary.overall_gc_percentage:.1f}%",
            annotation_position="top right"
        )
        fig.update_layout(
            xaxis_range=[0, 100],
            **THEME_LAYOUT
        )

    fig.update_xaxes(showgrid=True, gridcolor="#2D3748")
    fig.update_yaxes(showgrid=True, gridcolor="#2D3748")
    return fig


def plot_gc_vs_at_comparison(summary: DatasetSummary) -> go.Figure:
    """
    Generates a modern donut chart showing the balance between GC and AT content.
    """
    labels = ["GC Content", "AT Content"]
    values = [summary.overall_gc_percentage, summary.overall_at_percentage]
    colors = ["#3498DB", "#E67E22"]

    fig = go.Figure(
        data=[
            go.Pie(
                labels=labels,
                values=values,
                hole=0.6,
                marker=dict(colors=colors, line=dict(color="#111827", width=2)),
                textinfo="label+percent",
                hovertemplate="<b>%{label}</b><br>Proportion: %{percent}<extra></extra>"
            )
        ]
    )

    fig.add_annotation(
        text=f"<b>{summary.overall_gc_percentage:.1f}%</b><br><span style='font-size:12px;color:#9CA3AF'>Overall GC</span>",
        x=0.5, y=0.5,
        font=dict(size=20, color="#FFFFFF"),
        showarrow=False
    )

    fig.update_layout(
        title="<b>GC vs AT Content Ratio</b>",
        **THEME_LAYOUT
    )
    return fig


def plot_sequence_lengths(df: pd.DataFrame, summary: DatasetSummary) -> go.Figure:
    """
    Displays sequence length comparison across records, with an average length reference line.
    """
    if df.empty:
        return go.Figure()

    if len(df) <= 30:
        fig = px.bar(
            df,
            x="id",
            y="length",
            labels={"id": "Sequence ID", "length": "Sequence Length (bp)"},
            title="<b>Sequence Length Comparison</b>",
            color_discrete_sequence=["#6366F1"]
        )
    else:
        fig = px.histogram(
            df,
            x="length",
            nbins=35,
            marginal="rug",
            labels={"length": "Sequence Length (bp)"},
            title="<b>Sequence Length Distribution</b>",
            color_discrete_sequence=["#6366F1"]
        )

    if len(df) <= 30:
        fig.add_hline(
            y=summary.avg_length,
            line_dash="dash",
            line_color="#F59E0B",
            annotation_text=f"Mean Length: {summary.avg_length:,.1f} bp",
            annotation_position="top right"
        )
    else:
        fig.add_vline(
            x=summary.avg_length,
            line_dash="dash",
            line_color="#F59E0B",
            annotation_text=f"Mean Length: {summary.avg_length:,.1f} bp",
            annotation_position="top right"
        )

    fig.update_layout(**THEME_LAYOUT)
    fig.update_xaxes(showgrid=True, gridcolor="#2D3748")
    fig.update_yaxes(showgrid=True, gridcolor="#2D3748")
    return fig


def plot_gc_vs_length_scatter(df: pd.DataFrame) -> go.Figure:
    """
    QC scatter plot: Sequence Length vs GC Percentage.
    Crucial for identifying anomalous sequences, plasmid contamination, or divergent genomic regions.
    """
    if df.empty:
        return go.Figure()

    fig = px.scatter(
        df,
        x="length",
        y="gc_percentage",
        color="status",
        hover_name="id",
        hover_data=["description", "count_canonical", "count_ambiguous", "count_invalid"],
        labels={"length": "Sequence Length (bp)", "gc_percentage": "GC Content (%)", "status": "Validation Status"},
        title="<b>Genomic Quality Control: Sequence Length vs GC Content</b>"
    )

    fig.update_traces(marker=dict(size=9, opacity=0.85, line=dict(width=1, color="#FFFFFF")))
    fig.update_layout(
        yaxis_range=[0, 105],
        **THEME_LAYOUT
    )
    fig.update_xaxes(showgrid=True, gridcolor="#2D3748")
    fig.update_yaxes(showgrid=True, gridcolor="#2D3748")
    return fig
