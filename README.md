# 🧬 BioNexus: Industry-Oriented Modular Multi-Omics Bioinformatics Platform

**Version:** 0.1.0 (Phase 1: Genomics Foundation)  
**Status:** Active Development  
**Operating Environment:** Local / Windows 64-bit  

---

## 📌 Project Overview

**BioNexus** is a modular, high-throughput bioinformatics analysis platform designed for reproducible genomic, transcriptomic, proteomic, and systems biology research. Built entirely with open-source scientific Python technologies, BioNexus runs 100% locally with zero cloud dependencies, zero paid APIs, and complete data privacy.

The platform architecture is designed to progressively scale across the multi-omics landscape:
1. **Genomics (v0.1 — Active)**
2. **Transcriptomics (v0.2 — Coming Soon)**
3. **Proteomics (v0.3 — Coming Soon)**
4. **Metagenomics (v0.4 — Coming Soon)**
5. **Epigenomics (v0.5 — Coming Soon)**
6. **Multi-Omics Integration (v1.0 — Planned)**
7. **Biological Interpretation & Automated Reporting (v1.0 — Planned)**

---

## ✨ Features (Version 0.1 - Genomics Foundation)

### 1. Robust DNA FASTA Parsing
- Supports standard genomic FASTA formats: `.fa`, `.fasta`, `.fna`, `.txt`.
- Multi-line sequence aggregation and arbitrary whitespace tolerance.
- Multi-record dataset processing (thousands of contigs, genes, or chromosomes).
- Graceful error detection for non-FASTA files, empty files, or corrupted headers.

### 2. Comprehensive Nucleotide Stoichiometry & Quality Control
- **Dual Summary Framework:**
  - **Overall Dataset Overview:** Metrics across all records, including all-records average sequence length.
  - **Valid Records Baseline:** Dedicated quality metrics evaluating strictly valid sequences, including valid-records average length.
- **Canonical Nucleotide Counting:** Accurate counts for Adenine ($A$), Thymine ($T$), Guanine ($G$), and Cytosine ($C$), presented alongside canonical totals $(A + T + G + C)$.
- **Strict GC Percentage Calculation:**
  $$\text{GC\%} = \left(\frac{G + C}{A + T + G + C}\right) \times 100$$
  *The denominator is strictly defined as $A + T + G + C$. Ambiguous IUPAC bases ($N, R, Y, \dots$) and invalid characters are strictly excluded from the denominator.*
- **AT Percentage:** Complementary metric computed over canonical bases.
- **IUPAC Degenerate Base Auditing:** Full accounting for all IUPAC ambiguous codes ($R, Y, S, W, K, M, B, D, H, V, N$) with individual frequency breakdowns reported separately from canonical bases.
- **Exact Anomaly & Coordinate Tracking:**
  - Identifies exact sequence IDs containing non-IUPAC characters.
  - Tracks invalid characters and their exact 1-based positions in each record.
  - Identifies empty or zero-length FASTA records by sequence ID.
- **Downloadable QC Error Report:** Generates a dedicated CSV report detailing flagged sequence IDs, error categories, 1-based coordinate maps, and recommended remediation actions.

### 3. Interactive Data Visualizations (Plotly)
- **Overall Nucleotide Composition:** Proportional breakdown of canonical, ambiguous, and invalid bases.
- **Per-Sequence Nucleotide Stack:** Horizontal stacked bar chart comparing nucleotide distributions across records.
- **GC Content Distribution:** Interactive histogram / bar chart benchmarked against the dataset mean GC%.
- **GC vs AT Ratio:** High-contrast donut chart depicting overall base composition.
- **Sequence Length Comparison:** Sequence length bar chart and distribution with mean/median benchmarks.
- **QC Scatter Plot (Length vs GC%):** Scatter plot color-coded by validation status for outlier and anomaly detection.

### 4. Interactive Data Management & Reporting
- Filterable and searchable records table with column formatting.
- Instant CSV export of all computed metrics.
- Automatic snapshot persistence to `outputs/`.
- Pre-loaded biological datasets for instant one-click demonstration.

---

## 🛠️ Technology Stack

| Component | Technology | Description |
| :--- | :--- | :--- |
| **Language** | Python 3.12 | Core execution environment |
| **Interface** | Streamlit | Modern scientific dashboard |
| **Data Processing** | Pandas, NumPy | Vectorized sequence statistics |
| **Visualization** | Plotly | Publication-ready interactive visualizations |
| **Testing** | Pytest & Unittest | Automated unit and integration test suite |
| **Optional Bio** | Biopython | Available for downstream bioinformatics extensions |

---

## 📁 Project Structure

```text
BioNexus/
│
├── app.py                      # Main Streamlit application and view dispatcher
├── requirements.txt            # Project dependencies
├── README.md                   # Documentation and user guide
├── .gitignore                  # Git ignore rules
│
├── src/
│   ├── __init__.py             # Package initialization and version definition
│   ├── fasta_parser.py         # FASTA parser supporting multiline & edge cases
│   ├── sequence_analysis.py    # Nucleotide counting, IUPAC tracking & GC% formulas
│   └── visualization.py        # Interactive Plotly chart builders & styling
│
├── data/
│   ├── uploads/                # Directory for user uploads & sample presets
│   │   ├── sample_ecoli_genes.fasta
│   │   └── sample_ambiguous_and_edge.fasta
│   └── processed/              # Processed genomic datasets
│
├── outputs/                    # Output directory for exported CSV reports
│
├── tests/
│   ├── __init__.py             # Test package initialization
│   └── test_sequence_analysis.py # Automated test suite (26 tests)
│
└── docs/                       # Technical specifications & roadmap documents
    └── architecture.md
```

---

## 🚀 Installation & Setup

### 1. Prerequisites
- Windows 10/11 (or Linux/macOS)
- Python 3.10+ (Python 3.12 recommended)

### 2. Create Virtual Environment
Open PowerShell or Terminal in the project root:

```powershell
# Create virtual environment
py -3.12 -m venv .venv

# Activate virtual environment (Windows PowerShell)
.\.venv\Scripts\Activate.ps1
```

### 3. Install Dependencies
```powershell
# Install requirements
pip install -r requirements.txt
```

---

## 🖥️ Running the Application

Launch the Streamlit web interface using:

```powershell
streamlit run app.py
```

Or execute directly through Python in the virtual environment:

```powershell
.\.venv\Scripts\python.exe -m streamlit run app.py
```

Once running, navigate to `http://localhost:8501` in your browser.

---

## 🧪 Running Automated Tests

BioNexus includes a comprehensive test suite compatible with both `pytest` and Python's built-in `unittest`:

```powershell
# Run with Pytest
pytest tests/ -v

# Run with Python unittest
python -m unittest discover -s tests -p "test_*.py" -v
```

All 20 test cases validate:
- Single-line and multi-line FASTA records
- CRLF and LF newline handling
- Strict GC% formula and ambiguous base denominator exclusion
- Zero-length and empty sequences
- IUPAC ambiguous bases ($R, Y, S, W, K, M, B, D, H, V, N$)
- Non-IUPAC invalid character detection
- Dataset aggregation and DataFrame serialization

---

## 🔄 Current Development Status

- [x] **v0.1 - Genomics Foundation (Complete):** FASTA parser, GC% engine, Plotly visualizer, CSV export, Streamlit UI.
- [ ] **v0.2 - Transcriptomics (Planned):** RNA-Seq quantification, DESeq2 modeling, volcano plots.
- [ ] **v0.3 - Proteomics (Planned):** Peptide spectrum matching, protein quantitation.
- [ ] **v0.4 - Metagenomics (Planned):** 16S amplicon analysis, taxonomic abundance.
- [ ] **v0.5 - Epigenomics (Planned):** DNA methylation calling, peak exploration.
- [ ] **v1.0 - Multi-Omics Integration (Planned):** MOFA, cross-layer interpretation, PDF reports.

