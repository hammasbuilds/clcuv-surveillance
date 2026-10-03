"""Tests for the banded centre-star aligner."""

import importlib
import random

import pytest

from clcuv.align import (
    GAP,
    AlignmentImpossible,
    Scoring,
    align,
    align_pair,
    alignment_report,
    check_comparable,
    choose_centre,
    low_identity_rows,
    odd_prefixes,
    reverse_complement,
    stderr_progress,
)

BASES = "ACGT"


def _random_sequence(n: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice(BASES) for _ in range(n))


def _mutate(sequence: str, *, sites: int, seed: int) -> str:
    rng = random.Random(seed)
    chars = list(sequence)
    for position in rng.sample(range(len(chars)), sites):
        chars[position] = rng.choice([b for b in BASES if b != chars[position]])
    return "".join(chars)


# --- pairwise -------------------------------------------------------------


def test_identical_sequences_align_without_gaps():
    sequence = _random_sequence(300, seed=1)
    top, bottom = align_pair(sequence, sequence)
    assert top == bottom == sequence
    assert GAP not in top


def test_substitutions_do_not_introduce_gaps():
    a = _random_sequence(400, seed=2)
    b = _mutate(a, sites=12, seed=3)
    top, bottom = align_pair(a, b)
    assert len(top) == len(bottom) == 400
    assert GAP not in top and GAP not in bottom
    assert sum(1 for x, y in zip(top, bottom, strict=True) if x != y) == 12


def test_a_deletion_is_placed_as_a_gap():
    a = _random_sequence(300, seed=4)
    b = a[:150] + a[156:]  # six bases removed from the middle
    top, bottom = align_pair(a, b)
    assert len(top) == len(bottom)
    assert bottom.count(GAP) == 6
    assert top.replace(GAP, "") == a
    assert bottom.replace(GAP, "") == b


def test_an_insertion_gaps_the_other_row():
    a = _random_sequence(300, seed=5)
    b = a[:100] + "ACGTACGT" + a[100:]
    top, bottom = align_pair(a, b)
    assert top.count(GAP) == 8
    assert top.replace(GAP, "") == a
    assert bottom.replace(GAP, "") == b


def test_original_sequences_are_always_recoverable():
    """The invariant that makes an alignment an alignment: no base is lost."""
    a = _random_sequence(250, seed=6)
    b = _mutate(a[:120] + a[125:], sites=8, seed=7)
    top, bottom = align_pair(a, b)
    assert top.replace(GAP, "") == a
    assert bottom.replace(GAP, "") == b


def test_ambiguity_codes_are_neither_rewarded_nor_punished():
    scoring = Scoring()
    assert scoring.pair("A", "A") == scoring.match
    assert scoring.pair("A", "C") == scoring.mismatch
    assert scoring.pair("N", "A") == 0
    assert scoring.pair("A", "N") == 0


def test_empty_input_degrades_rather_than_crashing():
    top, bottom = align_pair("", "ACGT")
    assert top == GAP * 4
    assert bottom == "ACGT"


# --- the guard ------------------------------------------------------------


def test_rotated_circular_genomes_are_refused():
    """The failure this guard exists for: two rotations are biologically identical
    and share no aligned column, and nothing else in the pipeline would notice."""
    sequence = _random_sequence(600, seed=8)
    rotated = sequence[300:] + sequence[:300]
    with pytest.raises(AlignmentImpossible):
        check_comparable([sequence, sequence, rotated, rotated, rotated])


def test_co_oriented_genomes_pass_the_guard():
    a = _random_sequence(600, seed=9)
    check_comparable([a, _mutate(a, sites=30, seed=10), _mutate(a, sites=45, seed=11)])


def test_wildly_different_lengths_are_refused():
    a = _random_sequence(600, seed=12)
    with pytest.raises(AlignmentImpossible):
        check_comparable([a, a, a[:400]])


def test_a_single_sequence_is_not_an_alignment():
    with pytest.raises(AlignmentImpossible):
        check_comparable(["ACGT"])


def test_one_odd_sequence_among_many_does_not_trip_the_batch_guard():
    """80% agreement is the bar: one bad submission should not block the analysis.

    The batch guard is about the batch. Tolerating the odd sequence here is only
    defensible because `align` refuses to hand it back silently - see the tests below,
    which is the half that was missing.
    """
    a = _random_sequence(600, seed=13)
    others = [_mutate(a, sites=20, seed=s) for s in range(14, 24)]
    rotated = a[300:] + a[:300]
    check_comparable([a, *others, rotated])
    assert odd_prefixes([a, *others, rotated]) == [11]


def test_the_odd_sequence_is_reported_individually():
    a = _random_sequence(600, seed=13)
    others = [_mutate(a, sites=20, seed=s) for s in range(14, 24)]
    assert odd_prefixes([a, *others]) == []


def _batch_with_one_reverse_complement():
    a = _random_sequence(600, seed=40)
    others = [_mutate(a, sites=20, seed=s) for s in range(41, 49)]
    return [a, *others, reverse_complement(a)]


def test_a_reverse_complemented_submission_is_refused_not_aligned():
    """It passes the prefix guard as a minority of one, then aligns to noise. Every
    variant called from that row would be an artefact."""
    sequences = _batch_with_one_reverse_complement()
    check_comparable(sequences)  # the batch guard lets it through

    with pytest.raises(AlignmentImpossible) as caught:
        align(sequences)
    assert "identity to the centre" in str(caught.value)


def test_the_refusal_names_the_sequence_and_suggests_the_strand():
    sequences = _batch_with_one_reverse_complement()
    names = [f"ACC{i}" for i in range(len(sequences))]
    with pytest.raises(AlignmentImpossible) as caught:
        align(sequences, names=names)
    message = str(caught.value)
    assert "ACC9" in message
    assert "reverse complement" in message


def test_a_rotated_submission_is_refused_too():
    a = _random_sequence(600, seed=50)
    others = [_mutate(a, sites=20, seed=s) for s in range(51, 59)]
    rotated = a[300:] + a[:300]
    with pytest.raises(AlignmentImpossible, match="identity to the centre"):
        align([a, *others, rotated])


def test_the_junk_row_can_still_be_aligned_deliberately():
    """Refusing by default is right; refusing with no way through is not. The caller
    has now been told which row is junk."""
    sequences = _batch_with_one_reverse_complement()
    aligned = align(sequences, min_identity=0.0)
    assert len({len(row) for row in aligned}) == 1
    assert low_identity_rows(aligned, choose_centre(sequences))


def test_good_sequences_are_unaffected_by_the_identity_floor():
    a = _random_sequence(600, seed=60)
    sequences = [a, *[_mutate(a, sites=25, seed=s) for s in range(61, 66)]]
    assert len(align(sequences)) == 6


def test_names_must_match_the_sequences():
    a = _random_sequence(300, seed=70)
    with pytest.raises(AlignmentImpossible, match="names for"):
        align([a, a], names=["only-one"])


# --- multiple alignment ---------------------------------------------------


def test_centre_is_the_sequence_nearest_the_median_length():
    assert choose_centre(["A" * 100, "A" * 300, "A" * 305]) == 1


def test_multiple_alignment_returns_equal_lengths():
    a = _random_sequence(400, seed=30)
    sequences = [a, _mutate(a, sites=10, seed=31), a[:200] + a[205:], a + "ACGTA"]
    aligned = align(sequences)
    assert len({len(row) for row in aligned}) == 1


def test_multiple_alignment_preserves_every_base():
    a = _random_sequence(400, seed=32)
    sequences = [a, _mutate(a, sites=10, seed=33), a[:200] + a[205:], a + "ACGTA"]
    for original, row in zip(sequences, align(sequences), strict=True):
        assert row.replace(GAP, "") == original


def test_insertions_in_different_sequences_get_separate_columns():
    """The merge step's actual job: two sequences with insertions at the same place
    must not be forced to share one column, or a deletion appears from nowhere."""
    a = _random_sequence(300, seed=34)
    sequences = [a, a[:150] + "TTTT" + a[150:], a[:150] + "GGGGGG" + a[150:]]
    aligned = align(sequences)
    assert all(row.replace(GAP, "") == s for row, s in zip(aligned, sequences, strict=True))
    assert len(aligned[0]) >= 306


def test_identical_sequences_align_to_themselves():
    a = _random_sequence(200, seed=35)
    aligned = align([a, a, a])
    assert aligned == [a, a, a]


def test_report_counts_invariant_columns():
    a = _random_sequence(200, seed=36)
    report = alignment_report(align([a, a, a]))
    assert report["sequences"] == 3
    assert report["width"] == 200
    assert report["invariant_columns"] == 200
    assert report["invariant_fraction"] == 1.0
    assert report["gap_fraction"] == 0.0


def test_report_notices_variation():
    a = _random_sequence(300, seed=37)
    report = alignment_report(align([a, _mutate(a, sites=30, seed=38)]))
    assert report["invariant_columns"] == 270
    assert report["invariant_fraction"] < 1.0


def test_report_on_nothing():
    assert alignment_report([]) == {"sequences": 0}


# --- the optimisations do not move the answer -----------------------------
#
# `align_pair` is a hand-tuned banded DP with two shortcuts in it (identical sequences
# return immediately; duplicate sequences are aligned once). Each is justified by an
# argument about the scoring function, and an argument is not a test. So the banded
# implementation is checked against a plain, unbanded Needleman-Wunsch written for
# readability, with the same scores and the same tie-breaking.


def _reference_align_pair(a: str, b: str, scoring: Scoring | None = None) -> tuple[str, str]:
    """Textbook full-matrix Needleman-Wunsch. Slow, and obviously correct."""
    scoring = scoring or Scoring()
    n, m = len(a), len(b)
    score = [[0] * (m + 1) for _ in range(n + 1)]
    move = [[""] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        score[i][0] = score[i - 1][0] + scoring.gap
        move[i][0] = "U"
    for j in range(1, m + 1):
        score[0][j] = score[0][j - 1] + scoring.gap
        move[0][j] = "L"

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # Diagonal first, then up, then left - the same order of preference the
            # banded implementation uses, so ties land the same way.
            best = score[i - 1][j - 1] + scoring.pair(a[i - 1], b[j - 1])
            chosen = "D"
            if score[i - 1][j] + scoring.gap > best:
                best, chosen = score[i - 1][j] + scoring.gap, "U"
            if score[i][j - 1] + scoring.gap > best:
                best, chosen = score[i][j - 1] + scoring.gap, "L"
            score[i][j], move[i][j] = best, chosen

    i, j = n, m
    top: list[str] = []
    bottom: list[str] = []
    while i > 0 or j > 0:
        step = move[i][j]
        if step == "D":
            top.append(a[i - 1])
            bottom.append(b[j - 1])
            i -= 1
            j -= 1
        elif step == "U":
            top.append(a[i - 1])
            bottom.append(GAP)
            i -= 1
        else:
            top.append(GAP)
            bottom.append(b[j - 1])
            j -= 1
    return "".join(reversed(top)), "".join(reversed(bottom))


@pytest.mark.parametrize("seed", range(12))
def test_the_banded_aligner_agrees_with_full_needleman_wunsch(seed):
    rng = random.Random(seed)
    a = _random_sequence(rng.randrange(40, 160), seed=seed + 500)
    chars = list(a)
    for _ in range(rng.randrange(0, 8)):
        position = rng.randrange(len(chars))
        action = rng.choice("sid")
        if action == "s":
            chars[position] = rng.choice("ACGTN")
        elif action == "i":
            chars.insert(position, rng.choice(BASES))
        elif len(chars) > 2:
            del chars[position]
    b = "".join(chars)

    assert align_pair(a, b) == _reference_align_pair(a, b)


def test_the_identical_sequence_shortcut_matches_the_full_matrix():
    a = _random_sequence(120, seed=600)
    assert align_pair(a, a) == _reference_align_pair(a, a)


def test_duplicate_sequences_are_aligned_once(monkeypatch):
    # `clcuv.align` the function is re-exported from the package and shadows
    # `clcuv.align` the module, so the module has to be fetched by name.
    module = importlib.import_module("clcuv.align")

    calls = []
    real = module.align_pair

    def counted(x, y, **kwargs):
        calls.append(y)
        return real(x, y, **kwargs)

    monkeypatch.setattr(module, "align_pair", counted)

    a = _random_sequence(200, seed=601)
    b = _mutate(a, sites=5, seed=602)
    aligned = align([a, b, b, b, a])

    assert aligned[1] == aligned[2] == aligned[3]
    assert aligned[0] == aligned[4]
    # Four non-centre rows, two distinct sequences among them.
    assert len(set(calls)) == len(calls) <= 2


def test_progress_is_reported_once_per_sequence():
    a = _random_sequence(150, seed=603)
    sequences = [a, _mutate(a, sites=3, seed=604), _mutate(a, sites=4, seed=605)]
    seen = []
    align(sequences, progress=lambda done, total: seen.append((done, total)))
    assert seen == [(1, 3), (2, 3), (3, 3)]


# --- progress reporting ---------------------------------------------------


class _Stream:
    def __init__(self, tty: bool):
        self.tty = tty
        self.text = ""

    def isatty(self):
        return self.tty

    def write(self, text):
        self.text += text

    def flush(self):
        pass


def test_progress_on_a_captured_stream_is_one_line_per_tenth():
    stream = _Stream(tty=False)
    report = stderr_progress(50, stream=stream)
    for done in range(1, 51):
        report(done, 50)
    lines = [line for line in stream.text.splitlines() if line.strip()]
    assert len(lines) == 10
    assert "\r" not in stream.text
    assert lines[0].strip().startswith("aligning 5/50 (10%)")
    assert lines[-1].strip().startswith("aligning 50/50 (100%)")


def test_progress_on_a_terminal_rewrites_one_line():
    stream = _Stream(tty=True)
    report = stderr_progress(3, stream=stream)
    for done in range(1, 4):
        report(done, 3)
    assert stream.text.count("\r") == 3
    assert "aligning 1/3" in stream.text
    assert stream.text.endswith("\n")


def test_progress_always_reports_completion():
    stream = _Stream(tty=False)
    report = stderr_progress(7, stream=stream)
    for done in range(1, 8):
        report(done, 7)
    assert "aligning 7/7 (100%)" in stream.text
