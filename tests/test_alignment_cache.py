"""The real-data script's alignment cache: a cache that cannot change the answer.

`scripts/real_data.py` is the only place that touches the network, which is why it is a
script and not a module. The caching it does is pure, though, and a cache that can
serve a wrong alignment would quietly invalidate every number the README quotes - so
the pure part is tested here: the key, the verification, and the misses.
"""

from __future__ import annotations

import gzip
import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


@pytest.fixture
def real_data(tmp_path, monkeypatch):
    """The script, loaded as a module, with its data directory moved into tmp_path."""
    spec = importlib.util.spec_from_file_location(
        "real_data_under_test", ROOT / "scripts" / "real_data.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    monkeypatch.setattr(module, "DATA", tmp_path)
    monkeypatch.setattr(module, "ALIGNMENT_CACHE", tmp_path / "clcuv_aligned.fasta.gz")
    yield module
    del sys.modules[spec.name]


SEQUENCES = ["ACGTACGTAA", "ACGTACGTAA", "ACGTCGTAA"]
ALIGNED = ["ACGTACGTAA", "ACGTACGTAA", "ACGT-CGTAA"]


def test_a_missing_cache_is_a_miss(real_data):
    assert real_data._read_cached_alignment(SEQUENCES) is None


def test_a_written_cache_round_trips(real_data):
    real_data._write_cached_alignment(SEQUENCES, ALIGNED)
    assert real_data._read_cached_alignment(SEQUENCES) == ALIGNED


def test_a_different_corpus_does_not_hit_the_cache(real_data):
    real_data._write_cached_alignment(SEQUENCES, ALIGNED)
    assert real_data._read_cached_alignment([*SEQUENCES, "ACGTACGTAA"]) is None
    assert real_data._read_cached_alignment(list(reversed(SEQUENCES))) is None
    assert real_data._read_cached_alignment(["ACGTACGTAA", "ACGTACGTAA", "ACGTCGTAT"]) is None


def test_rows_that_are_not_the_sequences_they_claim_are_rejected(real_data):
    """The digest alone would accept this; the row check is what refuses it."""
    real_data._write_cached_alignment(SEQUENCES, ALIGNED)
    with gzip.open(real_data.ALIGNMENT_CACHE, "rt", encoding="utf-8") as handle:
        header = handle.readline()
    tampered = ["ACGTACGTAA", "ACGTACGTAA", "ACGTCGTAAA"]  # third row no longer matches
    with gzip.open(real_data.ALIGNMENT_CACHE, "wt", encoding="utf-8", newline="\n") as handle:
        handle.write(header)
        handle.writelines(f"{row}\n" for row in tampered)
    assert real_data._read_cached_alignment(SEQUENCES) is None


def test_ragged_rows_are_rejected(real_data):
    """Rows of unequal length are not an alignment, whatever the key says."""
    real_data._write_cached_alignment(SEQUENCES, ["ACGTACGTAA", "ACGTACGTAA", "ACGTCGTAA"])
    assert real_data._read_cached_alignment(SEQUENCES) is None


def test_a_corrupt_file_is_a_miss_not_a_crash(real_data):
    real_data.ALIGNMENT_CACHE.write_bytes(b"not gzip at all")
    assert real_data._read_cached_alignment(SEQUENCES) is None


def test_the_cache_is_used_instead_of_realigning(real_data, monkeypatch):
    real_data._write_cached_alignment(SEQUENCES, ALIGNED)
    monkeypatch.setattr(
        real_data, "align", lambda *a, **k: pytest.fail("aligned despite a valid cache")
    )
    assert real_data.aligned_corpus(SEQUENCES) == ALIGNED


def test_no_cache_realigns_and_does_not_read_the_file(real_data, monkeypatch):
    real_data._write_cached_alignment(SEQUENCES, ALIGNED)
    monkeypatch.setattr(real_data, "align", lambda *a, **k: ["A", "B", "C"])
    assert real_data.aligned_corpus(SEQUENCES, use_cache=False) == ["A", "B", "C"]
    # and the bogus result was not written over the good cache
    assert real_data._read_cached_alignment(SEQUENCES) == ALIGNED


def test_a_cold_run_writes_the_cache(real_data, monkeypatch):
    calls = []

    def fake_align(sequences, **kwargs):
        calls.append(sequences)
        return ALIGNED

    monkeypatch.setattr(real_data, "align", fake_align)
    assert real_data.aligned_corpus(SEQUENCES) == ALIGNED
    assert len(calls) == 1
    monkeypatch.setattr(
        real_data, "align", lambda *a, **k: pytest.fail("second run realigned anyway")
    )
    assert real_data.aligned_corpus(SEQUENCES) == ALIGNED
