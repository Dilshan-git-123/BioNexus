"""
BioNexus Sequence Analysis Module
=================================
Performs industry-standard genomics sequence analysis on DNA sequences:
- Nucleotide counting (A, T, G, C)
- IUPAC degenerate / ambiguous base identification & counting
- Non-IUPAC invalid character detection & exact 1-based position tracking
- Strict GC% calculation excluding ambiguous bases and invalid characters from denominator
- AT% calculation
- Separate Valid Records Summary vs All Records Summary
- Dedicated Quality Control (QC) anomaly reporting and CSV export
"""

from collections import Counter
from dataclasses import asdict, dataclass, field
from typing import Dict, List, Optional, Set, Tuple
import pandas as pd

from src.fasta_parser import FastaRecord


# Standard canonical DNA nucleotides
CANONICAL_BASES: Set[str] = {"A", "T", "G", "C"}

# Standard IUPAC ambiguous / degenerate DNA nucleotide codes
# R: A or G
# Y: C or T
# S: G or C
# W: A or T
# K: G or T
# M: A or C
# B: C or G or T (not A)
# D: A or G or T (not C)
# H: A or C or T (not G)
# V: A or C or G (not T)
# N: any base (A, C, G, T)
IUPAC_AMBIGUOUS_BASES: Set[str] = {
    "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N"
}

# All recognized IUPAC DNA characters
ALL_VALID_DNA_BASES: Set[str] = CANONICAL_BASES.union(IUPAC_AMBIGUOUS_BASES)


def calculate_gc_percentage(
    count_g: int, count_c: int, count_a: int, count_t: int
) -> float:
    """
    Calculates GC percentage according to the strict formula:
        GC% = ((G + C) / (A + T + G + C)) * 100

    CRITICAL RULE:
    The denominator is strictly (A + T + G + C).
    Ambiguous bases (N, R, Y, etc.) and invalid characters are strictly
    excluded from the denominator.
    If the canonical base count (A + T + G + C) is 0, returns 0.0.

    Parameters
    ----------
    count_g : int
        Number of Guanine bases.
    count_c : int
        Number of Cytosine bases.
    count_a : int
        Number of Adenine bases.
    count_t : int
        Number of Thymine bases.

    Returns
    -------
    float
        GC percentage rounded to 2 decimal places.
    """
    canonical_sum = count_a + count_t + count_g + count_c
    if canonical_sum == 0:
        return 0.0
    gc_pct = ((count_g + count_c) / canonical_sum) * 100.0
    return round(gc_pct, 2)


def calculate_at_percentage(
    count_a: int, count_t: int, count_g: int, count_c: int
) -> float:
    """
    Calculates AT percentage:
        AT% = ((A + T) / (A + T + G + C)) * 100

    Ambiguous bases and invalid characters are strictly excluded from the denominator.
    """
    canonical_sum = count_a + count_t + count_g + count_c
    if canonical_sum == 0:
        return 0.0
    at_pct = ((count_a + count_t) / canonical_sum) * 100.0
    return round(at_pct, 2)


@dataclass
class SequenceMetrics:
    """Detailed statistical and compositional metrics for a single sequence record."""
    id: str
    description: str
    length: int
    count_A: int
    count_T: int
    count_G: int
    count_C: int
    count_canonical: int
    count_ambiguous: int
    ambiguous_breakdown: Dict[str, int]
    count_invalid: int
    invalid_characters: List[str]
    invalid_positions: List[Tuple[int, str]]
    invalid_positions_summary: str
    gc_percentage: float
    at_percentage: float
    is_valid: bool
    status: str

    def to_dict(self) -> Dict:
        """Converts metrics to a flat dictionary suitable for DataFrame serialization."""
        d = asdict(self)
        d["invalid_characters"] = ", ".join(self.invalid_characters) if self.invalid_characters else "None"
        d["invalid_positions_summary"] = self.invalid_positions_summary if self.invalid_positions_summary else "None"
        d["ambiguous_breakdown"] = (
            ", ".join(f"{k}:{v}" for k, v in sorted(self.ambiguous_breakdown.items()))
            if self.ambiguous_breakdown else "None"
        )
        # Drop raw tuple list to keep DataFrame output clean and string-serializable
        d.pop("invalid_positions", None)
        return d


@dataclass
class DatasetSummary:
    """Aggregated statistics for a multi-sequence FASTA dataset."""
    filename: str
    total_sequences: int
    valid_sequences: int
    invalid_sequences: int
    empty_sequences: int
    total_length: int
    avg_length: float           # Average sequence length across ALL records (Requirement 3)
    all_avg_length: float       # Explicit alias for clarity (Requirement 3)
    valid_avg_length: float     # Average sequence length for VALID records only (Requirement 2)
    median_length: float
    min_length: int
    max_length: int
    total_A: int
    total_T: int
    total_G: int
    total_C: int
    total_canonical: int
    total_ambiguous: int
    total_invalid: int
    overall_gc_percentage: float
    overall_at_percentage: float
    # Valid records specific baseline (Requirement 1)
    valid_total_length: int = 0
    valid_median_length: float = 0.0
    valid_min_length: int = 0
    valid_max_length: int = 0
    valid_total_A: int = 0
    valid_total_T: int = 0
    valid_total_G: int = 0
    valid_total_C: int = 0
    valid_total_canonical: int = 0
    valid_total_ambiguous: int = 0
    valid_gc_percentage: float = 0.0
    valid_at_percentage: float = 0.0
    # QC identifier tracking lists (Requirements 4, 6)
    invalid_sequence_ids: List[str] = field(default_factory=list)
    empty_sequence_ids: List[str] = field(default_factory=list)
    ambiguous_sequence_ids: List[str] = field(default_factory=list)
    sequence_metrics: List[SequenceMetrics] = field(default_factory=list)

    def to_dataframe(self) -> pd.DataFrame:
        """Converts individual sequence metrics into a pandas DataFrame."""
        if not self.sequence_metrics:
            return pd.DataFrame()
        records = [m.to_dict() for m in self.sequence_metrics]
        df = pd.DataFrame(records)
        columns_order = [
            "id", "description", "length", "gc_percentage", "at_percentage",
            "count_A", "count_T", "count_G", "count_C",
            "count_canonical", "count_ambiguous", "ambiguous_breakdown",
            "count_invalid", "invalid_characters", "invalid_positions_summary",
            "is_valid", "status"
        ]
        present_cols = [c for c in columns_order if c in df.columns]
        return df[present_cols]

    def generate_qc_report_dataframe(self) -> pd.DataFrame:
        """
        Builds a dedicated QC error report DataFrame suitable for inspection and CSV export.
        Flags records containing invalid characters, empty sequences, or ambiguous IUPAC bases.
        If all records pass QC cleanly, generates a certified clean audit report.
        """
        qc_rows = []
        for m in self.sequence_metrics:
            issues = []
            recommendations = []
            error_types = []

            if m.length == 0:
                issues.append("Empty Sequence (0 bp)")
                error_types.append("EMPTY_SEQUENCE")
                recommendations.append("Remove empty header record from FASTA")

            if m.count_invalid > 0:
                issues.append(f"{m.count_invalid} invalid character(s): {', '.join(m.invalid_characters)}")
                error_types.append("INVALID_CHARACTERS")
                recommendations.append(f"Strip or correct illegal characters at positions: {m.invalid_positions_summary}")

            if m.count_ambiguous > 0:
                breakdown_str = ", ".join(f"{k}:{v}" for k, v in sorted(m.ambiguous_breakdown.items()))
                issues.append(f"{m.count_ambiguous} IUPAC ambiguous base(s): {breakdown_str}")
                error_types.append("AMBIGUOUS_BASES")
                recommendations.append("Review degenerate primer or consensus sequencing tolerance")

            if issues:
                qc_rows.append({
                    "sequence_id": m.id,
                    "description": m.description,
                    "sequence_length": m.length,
                    "qc_status": "FLAGGED",
                    "error_types": "; ".join(error_types),
                    "issue_details": "; ".join(issues),
                    "invalid_character_count": m.count_invalid,
                    "invalid_characters": ", ".join(m.invalid_characters) if m.invalid_characters else "None",
                    "invalid_positions": m.invalid_positions_summary,
                    "ambiguous_count": m.count_ambiguous,
                    "canonical_count": m.count_canonical,
                    "gc_percentage": m.gc_percentage,
                    "recommended_action": "; ".join(recommendations)
                })

        if not qc_rows:
            return pd.DataFrame([{
                "sequence_id": "ALL_RECORDS",
                "description": f"Dataset: {self.filename}",
                "sequence_length": self.total_length,
                "qc_status": "PASSED",
                "error_types": "NONE",
                "issue_details": f"All {self.total_sequences} sequences are canonical, non-empty, and valid.",
                "invalid_character_count": 0,
                "invalid_characters": "None",
                "invalid_positions": "None",
                "ambiguous_count": 0,
                "canonical_count": self.total_canonical,
                "gc_percentage": self.overall_gc_percentage,
                "recommended_action": "Dataset is certified for downstream alignment, phylogenetics, or variant calling."
            }])

        return pd.DataFrame(qc_rows)


def analyze_sequence(record: FastaRecord) -> SequenceMetrics:
    """
    Analyzes a single FastaRecord and computes nucleotide counts,
    GC/AT percentages, IUPAC ambiguous bases, and invalid character occurrences with exact positions.

    Parameters
    ----------
    record : FastaRecord
        The sequence record to analyze.

    Returns
    -------
    SequenceMetrics
        Computed metrics, exact positions, and validation status.
    """
    seq = record.sequence.upper()
    length = len(seq)

    if length == 0:
        return SequenceMetrics(
            id=record.id,
            description=record.description,
            length=0,
            count_A=0,
            count_T=0,
            count_G=0,
            count_C=0,
            count_canonical=0,
            count_ambiguous=0,
            ambiguous_breakdown={},
            count_invalid=0,
            invalid_characters=[],
            invalid_positions=[],
            invalid_positions_summary="None",
            gc_percentage=0.0,
            at_percentage=0.0,
            is_valid=False,
            status="Empty Sequence"
        )

    counts = Counter(seq)
    count_a = counts.get("A", 0)
    count_t = counts.get("T", 0)
    count_g = counts.get("G", 0)
    count_c = counts.get("C", 0)
    canonical_count = count_a + count_t + count_g + count_c

    # Ambiguous IUPAC bases
    ambiguous_breakdown = {
        base: counts[base] for base in IUPAC_AMBIGUOUS_BASES if base in counts and counts[base] > 0
    }
    count_ambiguous = sum(ambiguous_breakdown.values())

    # Invalid characters & exact 1-based positions
    invalid_chars_set = set()
    invalid_positions: List[Tuple[int, str]] = []
    count_invalid = 0

    for idx, char in enumerate(seq):
        if char not in ALL_VALID_DNA_BASES:
            # 1-based indexing for standard biological coordinate reporting
            pos = idx + 1
            invalid_positions.append((pos, char))
            invalid_chars_set.add(char)
            count_invalid += 1

    # Format invalid positions summary string
    if invalid_positions:
        pos_entries = [f"'{char}' at pos {pos}" for pos, char in invalid_positions]
        if len(pos_entries) <= 15:
            invalid_positions_summary = "; ".join(pos_entries)
        else:
            invalid_positions_summary = (
                "; ".join(pos_entries[:15]) + f" ... (+{len(pos_entries) - 15} more)"
            )
    else:
        invalid_positions_summary = "None"

    # Strict GC and AT percentages: denominator is strictly canonical_count (A + T + G + C)
    # Both ambiguous bases and invalid characters are excluded from denominator
    gc_pct = calculate_gc_percentage(count_g, count_c, count_a, count_t)
    at_pct = calculate_at_percentage(count_a, count_t, count_g, count_c)

    # Determine validation status
    if count_invalid > 0:
        is_valid = False
        status = f"Invalid Characters Detected ({', '.join(sorted(invalid_chars_set))})"
    elif count_ambiguous > 0:
        is_valid = True
        status = "Valid (Contains Ambiguous Bases)"
    else:
        is_valid = True
        status = "Valid (Canonical DNA)"

    return SequenceMetrics(
        id=record.id,
        description=record.description,
        length=length,
        count_A=count_a,
        count_T=count_t,
        count_G=count_g,
        count_C=count_c,
        count_canonical=canonical_count,
        count_ambiguous=count_ambiguous,
        ambiguous_breakdown=ambiguous_breakdown,
        count_invalid=count_invalid,
        invalid_characters=sorted(list(invalid_chars_set)),
        invalid_positions=invalid_positions,
        invalid_positions_summary=invalid_positions_summary,
        gc_percentage=gc_pct,
        at_percentage=at_pct,
        is_valid=is_valid,
        status=status
    )


def analyze_fasta_records(
    records: List[FastaRecord], filename: str = "dataset.fasta"
) -> DatasetSummary:
    """
    Performs comprehensive batch analysis across all sequence records in a FASTA file.
    Computes both overall dataset metrics and a separate Valid Records baseline.

    Parameters
    ----------
    records : List[FastaRecord]
        Parsed sequence records.
    filename : str
        Name of the source file.

    Returns
    -------
    DatasetSummary
        Aggregated dataset metrics, separate valid records metrics, and individual sequence results.
    """
    if not records:
        return DatasetSummary(
            filename=filename,
            total_sequences=0,
            valid_sequences=0,
            invalid_sequences=0,
            empty_sequences=0,
            total_length=0,
            avg_length=0.0,
            all_avg_length=0.0,
            valid_avg_length=0.0,
            median_length=0.0,
            min_length=0,
            max_length=0,
            total_A=0,
            total_T=0,
            total_G=0,
            total_C=0,
            total_canonical=0,
            total_ambiguous=0,
            total_invalid=0,
            overall_gc_percentage=0.0,
            overall_at_percentage=0.0,
            valid_total_length=0,
            valid_median_length=0.0,
            valid_min_length=0,
            valid_max_length=0,
            valid_total_A=0,
            valid_total_T=0,
            valid_total_G=0,
            valid_total_C=0,
            valid_total_canonical=0,
            valid_total_ambiguous=0,
            valid_gc_percentage=0.0,
            valid_at_percentage=0.0,
            invalid_sequence_ids=[],
            empty_sequence_ids=[],
            ambiguous_sequence_ids=[],
            sequence_metrics=[]
        )

    metrics_list: List[SequenceMetrics] = []
    all_lengths: List[int] = []
    valid_lengths: List[int] = []

    valid_seq_count = 0
    invalid_seq_count = 0
    empty_seq_count = 0

    invalid_sequence_ids: List[str] = []
    empty_sequence_ids: List[str] = []
    ambiguous_sequence_ids: List[str] = []

    # Overall totals (all records)
    total_a = 0
    total_t = 0
    total_g = 0
    total_c = 0
    total_ambiguous = 0
    total_invalid = 0

    # Valid records totals (Requirement 1)
    valid_total_a = 0
    valid_total_t = 0
    valid_total_g = 0
    valid_total_c = 0
    valid_total_ambiguous = 0

    for rec in records:
        m = analyze_sequence(rec)
        metrics_list.append(m)
        all_lengths.append(m.length)

        if m.length == 0:
            empty_seq_count += 1
            empty_sequence_ids.append(m.id)
        elif not m.is_valid:
            invalid_seq_count += 1
            invalid_sequence_ids.append(m.id)
        else:
            valid_seq_count += 1
            valid_lengths.append(m.length)
            valid_total_a += m.count_A
            valid_total_t += m.count_T
            valid_total_g += m.count_G
            valid_total_c += m.count_C
            valid_total_ambiguous += m.count_ambiguous

        if m.count_ambiguous > 0:
            ambiguous_sequence_ids.append(m.id)

        total_a += m.count_A
        total_t += m.count_T
        total_g += m.count_G
        total_c += m.count_C
        total_ambiguous += m.count_ambiguous
        total_invalid += m.count_invalid

    total_canonical = total_a + total_t + total_g + total_c
    total_length = sum(all_lengths)
    n_seq = len(records)

    # Average sequence length across ALL records (Requirement 3)
    avg_len = round(total_length / n_seq, 2) if n_seq > 0 else 0.0

    # Average sequence length for VALID records only (Requirement 2)
    valid_total_length = sum(valid_lengths)
    valid_avg_len = round(valid_total_length / valid_seq_count, 2) if valid_seq_count > 0 else 0.0

    # Median length (all records)
    sorted_all_lens = sorted(all_lengths)
    if n_seq % 2 == 1:
        median_len = float(sorted_all_lens[n_seq // 2])
    else:
        median_len = float(sorted_all_lens[n_seq // 2 - 1] + sorted_all_lens[n_seq // 2]) / 2.0

    min_len = min(all_lengths) if all_lengths else 0
    max_len = max(all_lengths) if all_lengths else 0

    # Valid records statistics
    if valid_lengths:
        sorted_valid_lens = sorted(valid_lengths)
        n_valid = len(valid_lengths)
        if n_valid % 2 == 1:
            valid_median_len = float(sorted_valid_lens[n_valid // 2])
        else:
            valid_median_len = float(sorted_valid_lens[n_valid // 2 - 1] + sorted_valid_lens[n_valid // 2]) / 2.0
        valid_min_len = min(valid_lengths)
        valid_max_len = max(valid_lengths)
    else:
        valid_median_len = 0.0
        valid_min_len = 0
        valid_max_len = 0

    valid_total_canonical = valid_total_a + valid_total_t + valid_total_g + valid_total_c
    valid_gc = calculate_gc_percentage(valid_total_g, valid_total_c, valid_total_a, valid_total_t)
    valid_at = calculate_at_percentage(valid_total_a, valid_total_t, valid_total_g, valid_total_c)

    overall_gc = calculate_gc_percentage(total_g, total_c, total_a, total_t)
    overall_at = calculate_at_percentage(total_a, total_t, total_g, total_c)

    return DatasetSummary(
        filename=filename,
        total_sequences=n_seq,
        valid_sequences=valid_seq_count,
        invalid_sequences=invalid_seq_count,
        empty_sequences=empty_seq_count,
        total_length=total_length,
        avg_length=avg_len,
        all_avg_length=avg_len,
        valid_avg_length=valid_avg_len,
        median_length=round(median_len, 2),
        min_length=min_len,
        max_length=max_len,
        total_A=total_a,
        total_T=total_t,
        total_G=total_g,
        total_C=total_c,
        total_canonical=total_canonical,
        total_ambiguous=total_ambiguous,
        total_invalid=total_invalid,
        overall_gc_percentage=overall_gc,
        overall_at_percentage=overall_at,
        valid_total_length=valid_total_length,
        valid_median_length=round(valid_median_len, 2),
        valid_min_length=valid_min_len,
        valid_max_length=valid_max_len,
        valid_total_A=valid_total_a,
        valid_total_T=valid_total_t,
        valid_total_G=valid_total_g,
        valid_total_C=valid_total_c,
        valid_total_canonical=valid_total_canonical,
        valid_total_ambiguous=valid_total_ambiguous,
        valid_gc_percentage=valid_gc,
        valid_at_percentage=valid_at,
        invalid_sequence_ids=invalid_sequence_ids,
        empty_sequence_ids=empty_sequence_ids,
        ambiguous_sequence_ids=ambiguous_sequence_ids,
        sequence_metrics=metrics_list
    )
