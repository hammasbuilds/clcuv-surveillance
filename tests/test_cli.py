"""Tests for the `clcuv` command line: a user's own FASTA and metadata sheet in.

These are the same paths a real user hits first - a file that does not exist,
a FASTA that is not FASTA, an unusable metadata join - and the whole point of
wrapping the library in a CLI was that those no longer crash with a traceback.
"""

from __future__ import annotations

import json

import pytest

from clcuv.cli import main

FASTA = """\
>A1
ATGCGTACGTTAGCATGCATCGATCGTAGCTAGCTAGCATCGATCGATCGTAGCATCGA
>A2
ATGCGTACGTTAGCATGCATCGATCGTAGCTAGCTAGCATCGATCGATCGTAGCATCGA
>A3
ATGCGTACGTTAGCATGCATCGATCGTAGCTAGCTAGCATCGATCGATCGTTGCATCGA
>B1
ATGCGTACGTTAGCATGCATCGATCGTAGCTAGCTAGCATCGATCGATCGTAGCATCGT
>B2
ATGCGTACGTTAGCATGCATCGATCGTAGCTAGCTAGCATCGATCGATCGTAGCATCGT
"""

METADATA = """\
name,date,location
A1,2019-03-01,Pakistan: Punjab
A2,2019-05-01,Pakistan: Punjab
A3,2019-08-01,Pakistan: Punjab
B1,2022-01-01,Pakistan: Sindh
B2,2022-06-01,Pakistan: Sindh
"""


@pytest.fixture
def fasta(tmp_path):
    path = tmp_path / "seqs.fasta"
    path.write_text(FASTA, encoding="utf-8")
    return path


@pytest.fixture
def metadata(tmp_path):
    path = tmp_path / "meta.csv"
    path.write_text(METADATA, encoding="utf-8")
    return path


class TestAnalyse:
    def test_runs_end_to_end_on_a_users_own_files(self, fasta, metadata, capsys):
        code = main(["analyse", "--fasta", str(fasta), "--metadata", str(metadata)])
        assert code == 0
        out = capsys.readouterr().out
        assert "5 sequences read, 5 surveillable" in out

    def test_counts_of_one_are_not_printed_as_plurals(self, fasta, metadata, capsys):
        """A stratum of one read `1 seqs -> 1 haplotypes`, which reads like a bug."""
        code = main(["analyse", "--fasta", str(fasta), "--metadata", str(metadata)])
        assert code == 0
        out = capsys.readouterr().out
        assert "1 seqs" not in out
        assert "1 haplotypes" not in out
        assert "1 haplotype " in out

    def test_json_output_is_valid_json_with_the_expected_shape(self, fasta, metadata, capsys):
        code = main(["analyse", "--fasta", str(fasta), "--metadata", str(metadata), "--json"])
        assert code == 0
        payload = json.loads(capsys.readouterr().out)
        assert payload["sequences_read"] == 5
        assert "sensitivity" in payload
        assert "all_sequences" in payload["sensitivity"]

    def test_missing_fasta_file_is_a_clean_error_not_a_traceback(self, tmp_path, metadata, capsys):
        code = main(
            ["analyse", "--fasta", str(tmp_path / "nope.fasta"), "--metadata", str(metadata)]
        )
        assert code == 1
        err = capsys.readouterr().err
        assert "no such file" in err.lower()

    def test_empty_fasta_is_a_clean_error(self, tmp_path, metadata, capsys):
        empty = tmp_path / "empty.fasta"
        empty.write_text("", encoding="utf-8")
        code = main(["analyse", "--fasta", str(empty), "--metadata", str(metadata)])
        assert code == 1
        assert "FASTA" in capsys.readouterr().err

    def test_metadata_is_required(self, fasta, capsys):
        with pytest.raises(SystemExit):
            main(["analyse", "--fasta", str(fasta)])

    def test_metadata_with_no_matching_names_fails_clearly(self, fasta, tmp_path):
        unrelated = tmp_path / "meta.csv"
        unrelated.write_text("name,date,location\nZZZ,2019,Pakistan\n", encoding="utf-8")
        with pytest.raises(SystemExit, match="not enough to compare"):
            main(["analyse", "--fasta", str(fasta), "--metadata", str(unrelated)])

    def test_duplicate_metadata_name_is_a_clean_error(self, fasta, tmp_path, capsys):
        dup = tmp_path / "meta.csv"
        dup.write_text("name,date,location\nA1,2019,Pakistan\nA1,2020,Pakistan\n", encoding="utf-8")
        code = main(["analyse", "--fasta", str(fasta), "--metadata", str(dup)])
        assert code == 1
        assert "twice" in capsys.readouterr().err

    def test_min_samples_sensitivity_table_is_present_in_json(self, fasta, metadata, capsys):
        code = main(
            [
                "analyse",
                "--fasta",
                str(fasta),
                "--metadata",
                str(metadata),
                "--min-samples",
                "2",
                "--json",
            ]
        )
        assert code == 0
        payload = json.loads(capsys.readouterr().out)
        rows = payload["sensitivity"]["all_sequences"]
        assert len(rows) > 1
        assert {row["min_samples"] for row in rows} == {r["min_samples"] for r in rows}


class TestExport:
    def test_writes_a_fasta_file(self, fasta, metadata, tmp_path, capsys):
        out = tmp_path / "out.fasta"
        code = main(
            ["export", "--fasta", str(fasta), "--metadata", str(metadata), "--out", str(out)]
        )
        assert code == 0
        assert out.exists()
        text = out.read_text(encoding="utf-8")
        assert text.startswith(";")
        assert "wrote" in capsys.readouterr().out

    def test_exported_file_is_valid_fasta_that_reads_back_in(self, fasta, metadata, tmp_path):
        out = tmp_path / "out.fasta"
        main(["export", "--fasta", str(fasta), "--metadata", str(metadata), "--out", str(out)])
        from clcuv.io import read_fasta

        records = read_fasta(out)
        assert len(records) >= 1
