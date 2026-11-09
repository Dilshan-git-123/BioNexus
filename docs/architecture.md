# BioNexus Technical Architecture & Scientific Specifications

## Overview
BioNexus is structured as an industry-standard multi-omics bioinformatics computational suite. It isolates parsing logic, biological computation, visualization routines, and interactive presentation into clean, decoupled modules.

## 1. Scientific Standards & Mathematical Specifications

### GC Percentage Formulation
According to molecular biology specifications for BioNexus:
$$GC\% = \left(\frac{G + C}{A + T + G + C}\right) \times 100$$

**Strict Mathematical Constraints:**
- **Denominator Definition:** Denominator is strictly the sum of canonical DNA bases $(A + T + G + C)$.
- **Exclusion Rule:** All IUPAC ambiguous/degenerate bases ($N, R, Y, S, W, K, M, B, D, H, V$) and non-IUPAC invalid characters are **strictly excluded from the denominator**.
- **Boundary Condition:** If $A + T + G + C = 0$, $GC\% = 0.0$ to prevent division by zero.

### IUPAC Nucleotide Codes Handled
| Code | Bases Represented | Type |
| :--- | :--- | :--- |
| **A** | Adenine | Canonical |
| **C** | Cytosine | Canonical |
| **G** | Guanine | Canonical |
| **T** | Thymine | Canonical |
| **R** | A or G (puRine) | Degenerate / Ambiguous |
| **Y** | C or T (pYrimidine) | Degenerate / Ambiguous |
| **S** | G or C (Strong) | Degenerate / Ambiguous |
| **W** | A or T (Weak) | Degenerate / Ambiguous |
| **K** | G or T (Keto) | Degenerate / Ambiguous |
| **M** | A or C (aMino) | Degenerate / Ambiguous |
| **B** | C, G, or T (not A) | Degenerate / Ambiguous |
| **D** | A, G, or T (not C) | Degenerate / Ambiguous |
| **H** | A, C, or T (not G) | Degenerate / Ambiguous |
| **V** | A, C, or G (not T) | Degenerate / Ambiguous |
| **N** | Any base (A, C, G, T) | Degenerate / Ambiguous |

Any character outside this set (such as numerals, whitespace within sequences, punctuation marks, or invalid letters like $Z, X, Q, J$) is registered as an invalid character, tracked by its exact 1-based coordinates, and isolated in the QC Error Report.

## 2. Quality Control (QC) Architecture

### Dual-Summary System
1. **Overall Dataset Overview (All Records):**
   - Captures entire input file characteristics including contig counts, total sequence length, all-records average sequence length, overall GC%, and overall AT%.
2. **Valid Records Quality Baseline:**
   - Evaluates only certified biological sequences (`is_valid == True` and `length > 0`).
   - Reports valid record count (% of dataset), valid total length, valid-records average length, valid median length, valid GC%, and valid AT%.

### Exact Anomaly Identification
- **Invalid Sequences Tracking:** Lists exact sequence IDs containing illegal characters.
- **1-Based Coordinate Mapping:** Identifies the exact 1-based position and character for every illegal token in a sequence record.
- **Empty Record Auditing:** Isolates zero-length or header-only FASTA records.
- **Downloadable QC Error Report:** Generates an industry-standard CSV audit trail detailing sequence ID, error type, invalid positions, and recommended remediation actions.

## 3. Directory Layout & Module Responsibilities

- **`app.py`**: Streamlit presentation tier, reactive state management, module dispatching, and user action handling.
- **`src/fasta_parser.py`**: Low-level stream ingestion, line aggregation, header parsing, and format validation.
- **`src/sequence_analysis.py`**: Core domain logic, mathematical formulation, nucleotide counting, IUPAC degenerate tracking, exact coordinate logging, and QC DataFrame generation.
- **`src/visualization.py`**: Interactive Plotly charts configured with scientific color scales and responsive themes.
- **`tests/test_sequence_analysis.py`**: Automated test suite (25 unit tests) covering parsing, mathematics, coordinate tracking, dual summaries, and QC reports.
- **`data/`**: Structured storage for uploads and processed datasets.
- **`outputs/`**: Local destination for exported CSV analysis reports and QC error audit trails.
