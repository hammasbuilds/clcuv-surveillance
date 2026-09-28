"""Tests for FASTA and metadata parsing.

Every case here is a way a real submission breaks: a missing sequence, a
duplicate name, a non-nucleotide character, a date format GenBank actually
uses, a metadata sheet that spells its columns differently.
"""

from __future__ import annotations

import pytest

from clcuv.io import (
    FastaError,
    isolates_from,
    parse_fasta,
    period_from,
    read_metadata,
    write_fasta,
)


class TestParseFasta:
    def test_reads_two_records(self):
        records = parse_fasta(">A1\nACGT\n>A2\nTTTT\n")
        assert [r.name for r in records] == ["A1", "A2"]
        assert [r.sequence for r in records] == ["ACGT", "TTTT"]

    def test_sequence_lines_are_joined_and_uppercased(self):
        (record,) = parse_fasta(">A1\nacgt\nACGT\n")
        assert record.sequence == "ACGTACGT"

    def test_header_description_after_the_name_is_kept(self):
        (record,) = parse_fasta(">A1 Pakistan, 2019\nACGT\n")
        assert record.name == "A1"
        assert record.description == "Pakistan, 2019"

    def test_semicolon_comment_lines_are_skipped(self):
        records = parse_fasta("; provenance header\n>A1\nACGT\n")
        assert len(records) == 1

    def test_empty_text_is_rejected(self):
        with pytest.raises(FastaError, match="no FASTA records"):
            parse_fasta("")

    def test_header_with_no_sequence_is_rejected(self):
        with pytest.raises(FastaError, match="header but no sequence"):
            parse_fasta(">A1\n>A2\nACGT\n")

    def test_sequence_before_any_header_is_rejected(self):
        with pytest.raises(FastaError, match="before the first"):
            parse_fasta("ACGT\n>A1\nACGT\n")

    def test_empty_header_name_is_rejected(self):
        with pytest.raises(FastaError, match="empty name"):
            parse_fasta(">\nACGT\n")

    def test_duplicate_name_is_rejected(self):
        with pytest.raises(FastaError, match="appears twice"):
            parse_fasta(">A1\nACGT\n>A1\nTTTT\n")

    def test_non_nucleotide_character_is_rejected(self):
        with pytest.raises(FastaError, match="not nucleotide"):
            parse_fasta(">A1\nACGZ\n")

    def test_protein_looking_sequence_gets_a_pointed_error(self):
        with pytest.raises(FastaError, match="protein"):
            parse_fasta(">A1\nMKVL\n")

    def test_gap_characters_kept_by_default(self):
        (record,) = parse_fasta(">A1\nAC-GT.N\n")
        assert record.sequence == "AC-GT.N"

    def test_gaps_rejected_when_unaligned_input_is_required(self):
        with pytest.raises(FastaError, match="gaps"):
            parse_fasta(">A1\nAC-GT\n", allow_gaps=False)

    def test_ambiguity_codes_are_accepted(self):
        (record,) = parse_fasta(">A1\nACGTRYKMSWBDHVN\n")
        assert record.sequence == "ACGTRYKMSWBDHVN"


class TestWriteFasta:
    def test_round_trips_through_parse_fasta(self, tmp_path):
        path = tmp_path / "out.fasta"
        n = write_fasta(path, [("A1", "ACGT" * 20), ("A2", "TTTT" * 20)])
        assert n == 2
        back = parse_fasta(path.read_text(encoding="utf-8"))
        assert [r.sequence for r in back] == ["ACGT" * 20, "TTTT" * 20]

    def test_wraps_at_the_given_width(self, tmp_path):
        path = tmp_path / "out.fasta"
        write_fasta(path, [("A1", "A" * 150)], width=70)
        lines = path.read_text(encoding="utf-8").splitlines()
        seq_lines = [line for line in lines if not line.startswith((">", ";"))]
        assert [len(line) for line in seq_lines] == [70, 70, 10]

    def test_writes_lf_not_crlf(self, tmp_path):
        path = tmp_path / "out.fasta"
        write_fasta(path, [("A1", "ACGT")])
        assert b"\r\n" not in path.read_bytes()

    def test_header_lines_are_prefixed_with_semicolon(self, tmp_path):
        path = tmp_path / "out.fasta"
        write_fasta(path, [("A1", "ACGT")], header=["one per haplotype"])
        first_line = path.read_text(encoding="utf-8").splitlines()[0]
        assert first_line == "; one per haplotype"

    def test_creates_missing_parent_directories(self, tmp_path):
        path = tmp_path / "nested" / "dir" / "out.fasta"
        write_fasta(path, [("A1", "ACGT")])
        assert path.exists()


class TestPeriodFrom:
    def test_bare_year(self):
        assert period_from("2019") == "2019"

    def test_month_year(self):
        assert period_from("May-2019") != ""

    def test_day_month_year(self):
        assert period_from("01-May-2019") != ""

    def test_iso_year_month_does_not_get_truncated_to_the_wrong_period(self):
        # The bug this guards: taking the last 4 characters of "2015-01" gives
        # "5-01", a nonsense period distinct from every other 2015 record.
        result = period_from("2015-01")
        assert result != "5-01"
        assert "2015" in result

    def test_blank_value_is_empty(self):
        assert period_from("") == ""
        assert period_from("   ") == ""

    def test_unparseable_text_with_no_year_is_empty(self):
        assert period_from("unknown") == ""


class TestReadMetadata:
    def test_reads_the_name_column_by_any_accepted_alias(self, tmp_path):
        path = tmp_path / "meta.csv"
        path.write_text("accession,date,country\nA1,2019,Pakistan\n", encoding="utf-8")
        meta = read_metadata(path)
        assert meta["A1"]["period"] == "2019"
        assert meta["A1"]["location"] == "Pakistan"

    def test_missing_name_column_is_a_clear_error(self, tmp_path):
        path = tmp_path / "meta.csv"
        path.write_text("date,country\n2019,Pakistan\n", encoding="utf-8")
        with pytest.raises(ValueError, match="no name column"):
            read_metadata(path)

    def test_duplicate_name_is_rejected(self, tmp_path):
        path = tmp_path / "meta.csv"
        path.write_text("name,date\nA1,2019\nA1,2020\n", encoding="utf-8")
        with pytest.raises(ValueError, match="appears twice"):
            read_metadata(path)

    def test_blank_name_rows_are_skipped_not_crashed_on(self, tmp_path):
        path = tmp_path / "meta.csv"
        path.write_text("name,date\n,2019\nA1,2020\n", encoding="utf-8")
        meta = read_metadata(path)
        assert list(meta) == ["A1"]

    def test_utf8_bom_is_handled(self, tmp_path):
        path = tmp_path / "meta.csv"
        path.write_bytes(b"\xef\xbb\xbfname,date\nA1,2019\n")
        meta = read_metadata(path)
        assert meta["A1"]["period"] == "2019"


class TestIsolatesFrom:
    def test_joins_records_to_metadata_by_name(self):
        records = parse_fasta(">A1\nACGT\n")
        metadata = {"A1": {"period": "2019", "location": "Pakistan: Punjab"}}
        isolates, skipped = isolates_from(records, metadata)
        assert len(isolates) == 1
        assert skipped == []
        assert isolates[0].period == "2019"

    def test_joins_by_the_first_token_when_the_header_has_a_description(self):
        records = parse_fasta(">MW183409.1 Punjab isolate\nACGT\n")
        metadata = {"MW183409.1": {"period": "2019", "location": "Pakistan"}}
        isolates, skipped = isolates_from(records, metadata)
        assert len(isolates) == 1

    def test_a_record_with_no_metadata_row_is_skipped_not_dropped_silently(self):
        records = parse_fasta(">A1\nACGT\n>A2\nTTTT\n")
        metadata = {"A1": {"period": "2019", "location": "Pakistan"}}
        isolates, skipped = isolates_from(records, metadata)
        assert len(isolates) == 1
        assert skipped == ["A2"]

    def test_missing_location_or_period_is_skipped(self):
        records = parse_fasta(">A1\nACGT\n")
        metadata = {"A1": {"period": "2019", "location": ""}}
        isolates, skipped = isolates_from(records, metadata)
        assert isolates == []
        assert skipped == ["A1"]
