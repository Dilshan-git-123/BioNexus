"""
BioNexus Automated Test Suite
=============================
Tests for FASTA parsing, nucleotide composition, strict GC percentage
calculation (with ambiguous base & invalid character exclusion from denominator),
multi-record support, empty sequences, ambiguous nucleotides, invalid character
position tracking, valid-record baseline separation, and QC error reports.

Compatible with both `pytest` and `python -m unittest`.
"""

import io
import unittest
import pandas as pd

from src.fasta_parser import FastaRecord, FastaParseError, parse_fasta
from src.sequence_analysis import (
    calculate_gc_percentage,
    calculate_at_percentage,
    analyze_sequence,
    analyze_fasta_records,
    IUPAC_AMBIGUOUS_BASES,
    CANONICAL_BASES
)


class TestFastaParser(unittest.TestCase):
    """Test suite for FASTA parsing capabilities."""

    def test_parse_single_record(self):
        """Verify parsing a basic single-record FASTA string."""
        raw = ">seq1 Human beta-globin segment\nATGCGATCGATC"
        records = parse_fasta(raw)
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].id, "seq1")
        self.assertEqual(records[0].description, "Human beta-globin segment")
        self.assertEqual(records[0].sequence, "ATGCGATCGATC")
        self.assertEqual(records[0].length, 12)

    def test_parse_multiline_sequence(self):
        """Verify sequence across multiple lines with arbitrary whitespaces is concatenated correctly."""
        raw = (
            ">gene_001 Poly-line test\n"
            "ATGC\n"
            "GATC\n"
            "  TTTT  \n"
            "CCCC\n"
        )
        records = parse_fasta(raw)
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].sequence, "ATGCGATCTTTTCCCC")
        self.assertEqual(records[0].length, 16)

    def test_parse_multiple_records(self):
        """Verify parsing a multi-record FASTA stream with blank lines."""
        raw = (
            ">rec1 First gene\n"
            "ATGC\n\n"
            ">rec2 Second gene with no description\n"
            "GGCC\n\n"
            ">rec3 Third gene\n"
            "AATT\n"
        )
        records = parse_fasta(raw)
        self.assertEqual(len(records), 3)
        self.assertEqual([r.id for r in records], ["rec1", "rec2", "rec3"])
        self.assertEqual([r.sequence for r in records], ["ATGC", "GGCC", "AATT"])

    def test_parse_crlf_newlines(self):
        """Verify Windows-style CRLF (\\r\\n) is handled seamlessly."""
        raw = ">win_seq Windows test\r\nATGCATGC\r\nGGCCAATT\r\n"
        records = parse_fasta(raw)
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].sequence, "ATGCATGCGGCCAATT")

    def test_parse_empty_input_raises_error(self):
        """Verify empty input triggers FastaParseError."""
        with self.assertRaises(FastaParseError):
            parse_fasta("   \n\n   ")

    def test_parse_invalid_header_raises_error(self):
        """Verify input without a leading '>' header triggers FastaParseError."""
        raw = "ATGCGATCGATC\nATGCATGC"
        with self.assertRaises(FastaParseError):
            parse_fasta(raw)

    def test_parse_file_like_object(self):
        """Verify parsing from StringIO or BytesIO handles file-like streams."""
        stream = io.StringIO(">io_seq Stream test\nACGTACGT\n")
        records = parse_fasta(stream)
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].id, "io_seq")
        self.assertEqual(records[0].sequence, "ACGTACGT")

    def test_parse_empty_sequence(self):
        """Verify a header followed immediately by another header or EOF produces an empty sequence record."""
        raw = ">empty_seq1 Empty sequence\n>next_seq Normal sequence\nATGC"
        records = parse_fasta(raw)
        self.assertEqual(len(records), 2)
        self.assertEqual(records[0].id, "empty_seq1")
        self.assertEqual(records[0].sequence, "")
        self.assertEqual(records[0].length, 0)
        self.assertEqual(records[1].id, "next_seq")
        self.assertEqual(records[1].sequence, "ATGC")


class TestSequenceAnalysis(unittest.TestCase):
    """Test suite for nucleotide counting, GC% calculation, and validation."""

    def test_gc_percentage_all_gc(self):
        """Verify 100% GC sequence."""
        gc = calculate_gc_percentage(count_g=50, count_c=50, count_a=0, count_t=0)
        self.assertEqual(gc, 100.0)

    def test_gc_percentage_all_at(self):
        """Verify 0% GC sequence."""
        gc = calculate_gc_percentage(count_g=0, count_c=0, count_a=50, count_t=50)
        self.assertEqual(gc, 0.0)

    def test_gc_percentage_balanced(self):
        """Verify balanced 50% GC sequence."""
        gc = calculate_gc_percentage(count_g=25, count_c=25, count_a=25, count_t=25)
        self.assertEqual(gc, 50.0)

    def test_gc_percentage_rounding(self):
        """Verify correct rounding to 2 decimal places."""
        gc = calculate_gc_percentage(count_g=1, count_c=0, count_a=1, count_t=1)
        self.assertEqual(gc, 33.33)

    def test_gc_percentage_zero_canonical(self):
        """Verify safe handling of zero canonical bases (avoid division by zero)."""
        gc = calculate_gc_percentage(count_g=0, count_c=0, count_a=0, count_t=0)
        self.assertEqual(gc, 0.0)

    def test_gc_percentage_excludes_ambiguous_bases_from_denominator(self):
        """
        CRITICAL REQUIREMENT:
        GC% = ((G + C) / (A + T + G + C)) * 100
        Ambiguous bases (such as N, R, Y) MUST be excluded from the denominator.
        """
        rec = FastaRecord(id="ambig_gc", description="Test", sequence="GGCCNNNN", raw_header="ambig_gc")
        metrics = analyze_sequence(rec)

        self.assertEqual(metrics.count_G, 2)
        self.assertEqual(metrics.count_C, 2)
        self.assertEqual(metrics.count_canonical, 4)
        self.assertEqual(metrics.count_ambiguous, 4)
        # GC% should be (4 / 4) * 100 = 100.0%, NOT (4 / 8) * 100 = 50.0%.
        self.assertEqual(metrics.gc_percentage, 100.0)
        self.assertEqual(metrics.at_percentage, 0.0)

    def test_gc_denominator_excludes_invalid_characters_and_ambiguous_bases(self):
        """
        REQUIREMENT 9 & 10:
        Exclude both invalid characters and ambiguous bases from the GC% denominator.
        Denominator = A + T + G + C.
        """
        # Sequence: "GGCC" + "NN" + "12@#" -> total length = 10
        # G=2, C=2, N=2, invalid=4 ('1','2','@','#')
        # Denominator = 4 (strictly canonical bases).
        # GC% = (4 / 4) * 100 = 100.0%
        rec = FastaRecord(id="mixed_edge", description="Test", sequence="GGCCNN12@#", raw_header="mixed_edge")
        metrics = analyze_sequence(rec)

        self.assertEqual(metrics.length, 10)
        self.assertEqual(metrics.count_canonical, 4)
        self.assertEqual(metrics.count_ambiguous, 2)
        self.assertEqual(metrics.count_invalid, 4)
        self.assertEqual(metrics.gc_percentage, 100.0)

    def test_nucleotide_counting(self):
        """Verify accurate counts for A, T, G, C."""
        rec = FastaRecord(id="count_test", description="", sequence="AAATTTGGCC", raw_header="count_test")
        m = analyze_sequence(rec)
        self.assertEqual(m.count_A, 3)
        self.assertEqual(m.count_T, 3)
        self.assertEqual(m.count_G, 2)
        self.assertEqual(m.count_C, 2)
        self.assertEqual(m.count_canonical, 10)
        self.assertEqual(m.count_ambiguous, 0)
        self.assertEqual(m.count_invalid, 0)
        self.assertTrue(m.is_valid)

    def test_case_insensitivity(self):
        """Verify lowercase input is converted and counted identically to uppercase."""
        raw = ">lower Test lowercase\naaatttggcc\n"
        records = parse_fasta(raw)
        m = analyze_sequence(records[0])
        self.assertEqual(m.count_A, 3)
        self.assertEqual(m.count_T, 3)
        self.assertEqual(m.count_G, 2)
        self.assertEqual(m.count_C, 2)
        self.assertEqual(m.gc_percentage, 40.0)

    def test_ambiguous_nucleotides_detection(self):
        """Verify all IUPAC degenerate bases (R, Y, S, W, K, M, B, D, H, V, N) are counted."""
        ambig_seq = "ACGT" + "RYSWKMBDHVN"
        rec = FastaRecord(id="all_ambig", description="", sequence=ambig_seq, raw_header="all_ambig")
        m = analyze_sequence(rec)

        self.assertEqual(m.count_canonical, 4)
        self.assertEqual(m.count_ambiguous, 11)
        self.assertEqual(m.count_invalid, 0)
        self.assertTrue(m.is_valid)
        self.assertIn("Ambiguous Bases", m.status)
        for base in "RYSWKMBDHVN":
            self.assertEqual(m.ambiguous_breakdown.get(base), 1)

    def test_invalid_sequence_characters_and_exact_positions(self):
        """
        REQUIREMENTS 4 & 5:
        Verify non-IUPAC characters are flagged and their exact 1-based positions recorded.
        """
        # "AT12GC@#" -> '1' at pos 3, '2' at pos 4, '@' at pos 7, '#' at pos 8
        seq_with_invalid = "AT12GC@#"
        rec = FastaRecord(id="invalid_rec", description="", sequence=seq_with_invalid, raw_header="invalid_rec")
        m = analyze_sequence(rec)

        self.assertEqual(m.count_canonical, 4)
        self.assertEqual(m.count_invalid, 4)
        self.assertEqual(set(m.invalid_characters), {"1", "2", "@", "#"})
        self.assertFalse(m.is_valid)
        self.assertEqual(m.invalid_positions, [(3, "1"), (4, "2"), (7, "@"), (8, "#")])
        self.assertIn("'1' at pos 3", m.invalid_positions_summary)
        self.assertIn("'@' at pos 7", m.invalid_positions_summary)

    def test_empty_sequence_metrics(self):
        """Verify empty sequence is handled cleanly without exceptions."""
        rec = FastaRecord(id="empty", description="", sequence="", raw_header="empty")
        m = analyze_sequence(rec)

        self.assertEqual(m.length, 0)
        self.assertEqual(m.count_canonical, 0)
        self.assertEqual(m.gc_percentage, 0.0)
        self.assertFalse(m.is_valid)
        self.assertEqual(m.status, "Empty Sequence")
        self.assertEqual(m.invalid_positions, [])

    def test_dataset_summary_aggregation(self):
        """Verify dataset summary aggregates metrics properly across multiple records."""
        records = [
            FastaRecord(id="rec1", description="Canonical", sequence="GGCC", raw_header="rec1"),  # GC: 100%
            FastaRecord(id="rec2", description="Canonical", sequence="AATT", raw_header="rec2"),  # GC: 0%
            FastaRecord(id="rec3", description="Ambiguous", sequence="GCNN", raw_header="rec3"),  # GC: 100%
            FastaRecord(id="rec4", description="Invalid", sequence="AT12", raw_header="rec4"),     # Invalid
            FastaRecord(id="rec5", description="Empty", sequence="", raw_header="rec5")            # Empty
        ]

        summary = analyze_fasta_records(records, filename="test_run.fasta")

        self.assertEqual(summary.filename, "test_run.fasta")
        self.assertEqual(summary.total_sequences, 5)
        self.assertEqual(summary.valid_sequences, 3)
        self.assertEqual(summary.invalid_sequences, 1)
        self.assertEqual(summary.empty_sequences, 1)

        self.assertEqual(summary.total_length, 16)
        self.assertEqual(summary.total_G, 3)
        self.assertEqual(summary.total_C, 3)
        self.assertEqual(summary.total_A, 3)
        self.assertEqual(summary.total_T, 3)
        self.assertEqual(summary.total_canonical, 12)
        self.assertEqual(summary.total_ambiguous, 2)
        self.assertEqual(summary.total_invalid, 2)

        # Overall GC: (3 G + 3 C) / 12 = 50.0%
        self.assertEqual(summary.overall_gc_percentage, 50.0)
        self.assertEqual(summary.overall_at_percentage, 50.0)

        # Test conversion to DataFrame
        df = summary.to_dataframe()
        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(len(df), 5)
        self.assertIn("gc_percentage", df.columns)
        self.assertIn("invalid_positions_summary", df.columns)
        self.assertIn("status", df.columns)

    def test_valid_avg_length_vs_all_avg_length(self):
        """
        REQUIREMENTS 1, 2, 3:
        Verify separate calculation of valid_avg_length vs all_avg_length.
        """
        # rec1: len 10 (valid)
        # rec2: len 20 (valid)
        # rec3: len 0 (empty)
        # rec4: len 10 (invalid)
        records = [
            FastaRecord(id="v1", description="", sequence="A" * 10, raw_header="v1"),
            FastaRecord(id="v2", description="", sequence="C" * 20, raw_header="v2"),
            FastaRecord(id="e1", description="", sequence="", raw_header="e1"),
            FastaRecord(id="inv1", description="", sequence="A" * 8 + "12", raw_header="inv1"),
        ]
        summary = analyze_fasta_records(records, filename="length_test.fasta")

        # Total sequences: 4. Total length: 10 + 20 + 0 + 10 = 40.
        # all_avg_length = 40 / 4 = 10.0
        self.assertEqual(summary.all_avg_length, 10.0)
        self.assertEqual(summary.avg_length, 10.0)

        # Valid sequences: 2 (v1, v2). Valid length: 10 + 20 = 30.
        # valid_avg_length = 30 / 2 = 15.0
        self.assertEqual(summary.valid_sequences, 2)
        self.assertEqual(summary.valid_total_length, 30)
        self.assertEqual(summary.valid_avg_length, 15.0)

    def test_anomaly_identification_lists(self):
        """
        REQUIREMENTS 4 & 6:
        Verify exact sequence IDs for invalid characters and empty records.
        """
        records = [
            FastaRecord(id="norm_gene", description="", sequence="ATGC", raw_header="norm_gene"),
            FastaRecord(id="empty_contig_1", description="", sequence="", raw_header="empty_contig_1"),
            FastaRecord(id="bad_read_a", description="", sequence="AT99GC", raw_header="bad_read_a"),
            FastaRecord(id="empty_contig_2", description="", sequence="", raw_header="empty_contig_2"),
            FastaRecord(id="bad_read_b", description="", sequence="GG##CC", raw_header="bad_read_b"),
            FastaRecord(id="ambig_read", description="", sequence="ATNNGC", raw_header="ambig_read"),
        ]
        summary = analyze_fasta_records(records, filename="id_test.fasta")

        self.assertEqual(summary.invalid_sequence_ids, ["bad_read_a", "bad_read_b"])
        self.assertEqual(summary.empty_sequence_ids, ["empty_contig_1", "empty_contig_2"])
        self.assertEqual(summary.ambiguous_sequence_ids, ["ambig_read"])

    def test_generate_qc_report_dataframe_with_flags(self):
        """
        REQUIREMENT 11:
        Verify generation of downloadable QC Error Report DataFrame for flagged records.
        """
        records = [
            FastaRecord(id="clean_gene", description="", sequence="ATGC", raw_header="clean_gene"),
            FastaRecord(id="corrupted_read", description="Has digits", sequence="AT12GC", raw_header="corrupted_read"),
            FastaRecord(id="empty_record", description="Blank", sequence="", raw_header="empty_record"),
        ]
        summary = analyze_fasta_records(records, filename="qc_flag_test.fasta")
        qc_df = summary.generate_qc_report_dataframe()

        self.assertIsInstance(qc_df, pd.DataFrame)
        self.assertEqual(len(qc_df), 2)  # corrupted_read and empty_record
        self.assertIn("sequence_id", qc_df.columns)
        self.assertIn("error_types", qc_df.columns)
        self.assertIn("invalid_positions", qc_df.columns)
        self.assertIn("recommended_action", qc_df.columns)

        # Check corrupted_read details
        corr_row = qc_df[qc_df["sequence_id"] == "corrupted_read"].iloc[0]
        self.assertEqual(corr_row["invalid_character_count"], 2)
        self.assertIn("pos 3", corr_row["invalid_positions"])

        # Check empty_record details
        empty_row = qc_df[qc_df["sequence_id"] == "empty_record"].iloc[0]
        self.assertIn("EMPTY_SEQUENCE", empty_row["error_types"])

    def test_generate_qc_report_dataframe_clean_dataset(self):
        """
        REQUIREMENT 11:
        Verify clean dataset produces certified PASSED audit row.
        """
        records = [
            FastaRecord(id="g1", description="", sequence="ATGC", raw_header="g1"),
            FastaRecord(id="g2", description="", sequence="CCGG", raw_header="g2"),
        ]
        summary = analyze_fasta_records(records, filename="clean.fasta")
        qc_df = summary.generate_qc_report_dataframe()

        self.assertEqual(len(qc_df), 1)
        self.assertEqual(qc_df.iloc[0]["qc_status"], "PASSED")
        self.assertEqual(qc_df.iloc[0]["error_types"], "NONE")

    def test_sample_ambiguous_and_edge_dataset_regression_audit(self):
        """
        REGRESSION AUDIT TEST:
        Validates exact numerical consistency for sample_ambiguous_and_edge.fasta.
        Ensures all-records average is 20.0 bp (120 / 6) and valid-records average
        is 25.25 bp (101 / 4).
        """
        from pathlib import Path
        sample_path = Path("data/uploads/sample_ambiguous_and_edge.fasta")
        with open(sample_path, "r", encoding="utf-8") as f:
            records = parse_fasta(f)

        self.assertEqual(len(records), 6)
        expected_lengths = [18, 19, 32, 32, 19, 0]
        actual_lengths = [r.length for r in records]
        self.assertEqual(actual_lengths, expected_lengths)

        summary = analyze_fasta_records(records, filename="sample_ambiguous_and_edge.fasta")

        # 1. Total sequence length
        self.assertEqual(summary.total_length, 120)
        # 2. Average across all 6 records: 120 / 6 = 20.0 bp
        self.assertEqual(summary.all_avg_length, 20.0)
        self.assertEqual(summary.avg_length, 20.0)
        # 3. Valid records count and valid average length: 101 / 4 = 25.25 bp
        self.assertEqual(summary.valid_sequences, 4)
        self.assertEqual(summary.valid_total_length, 101)
        self.assertEqual(summary.valid_avg_length, 25.25)
        # 4. Canonical base counts
        self.assertEqual(summary.total_A, 24)
        self.assertEqual(summary.total_T, 26)
        self.assertEqual(summary.total_G, 26)
        self.assertEqual(summary.total_C, 23)
        self.assertEqual(summary.total_canonical, 99)
        # 5. Ambiguous & invalid counts
        self.assertEqual(summary.total_ambiguous, 14)
        self.assertEqual(summary.total_invalid, 7)
        # 6. Overall & valid GC percentage
        self.assertEqual(summary.overall_gc_percentage, 49.49)
        self.assertEqual(summary.valid_gc_percentage, 49.44)


if __name__ == "__main__":
    unittest.main()
