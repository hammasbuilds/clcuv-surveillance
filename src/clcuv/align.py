"""Multiple alignment for closely related, co-oriented viral genomes.

Every other module in this package requires aligned sequences and says so. That is the
right boundary — MAFFT and MUSCLE exist and are better than anything written here — but
it left the package unable to touch a raw FASTA from NCBI without another tool
installed, which in practice means the analysis does not get run.

So this is a **narrow** aligner, correct for one situation and honest about the rest:

  closely related    strains of one virus species, typically 85-99% identical
  co-oriented        all sequences starting at the same genomic position
  similar length     differing by tens of bases, not hundreds

Begomovirus genomes satisfy all three when taken from GenBank, because submissions are
conventionally rotated to begin at the same point in the intergenic region. **This is
checked rather than assumed** — `check_comparable()` refuses sequences that do not share
a common prefix, because a rotated circular genome aligned naively produces a result that
looks fine and is meaningless.

Method: **centre-star with banded Needleman-Wunsch.** Every sequence is aligned
pairwise to one centre, then the pairwise alignments are merged. The band is what makes
it tractable in pure Python — a full DP over 2,700 × 2,700 is seven million cells per
pair, while a band of ±120 is six hundred thousand, and for sequences this similar the
optimal path never leaves the band.

It will do badly on divergent sequences, large rearrangements, or genomes that are not
co-oriented. Use a real aligner for those; that is what `check_comparable` is telling you.
"""

from __future__ import annotations

import statistics
from collections.abc import Sequence
from dataclasses import dataclass

GAP = "-"


class AlignmentImpossible(ValueError):
    """The inputs violate an assumption this aligner depends on."""


@dataclass(frozen=True)
class Scoring:
    match: int = 2
    mismatch: int = -3
    gap: int = -5

    def pair(self, a: str, b: str) -> int:
        if a == b and a in "ACGT":
            return self.match
        if a not in "ACGT" or b not in "ACGT":
            return 0  # ambiguity codes neither reward nor punish
        return self.mismatch


COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def odd_prefixes(sequences: Sequence[str], *, prefix: int = 24) -> list[int]:
    """Indices whose start disagrees with the majority start.

    Reported per sequence rather than as one batch verdict. The batch verdict alone was
    the bug: a single rotated or reverse-complemented submission among fifty good ones
    left 98% agreement, passed, and was aligned to noise.
    """
    if not sequences:
        return []
    prefixes = [s[:prefix].upper() for s in sequences]
    reference = max(set(prefixes), key=prefixes.count)
    return [i for i, p in enumerate(prefixes) if _identity(p, reference) < 0.7]


def check_comparable(
    sequences: Sequence[str],
    *,
    prefix: int = 24,
    min_agreement: float = 0.8,
    names: Sequence[str] | None = None,
) -> None:
    """Refuse input this aligner cannot handle.

    The failure being prevented is specific. A circular genome can be deposited starting
    at any point, and two rotations of the same genome are 100% identical biologically
    while sharing almost no aligned column. A naive alignment of them produces a dense
    field of apparent mutations, every downstream number is wrong, and nothing errors.

    This is the *batch* guard: it refuses a set that is mostly not co-oriented. It
    deliberately tolerates a minority of odd sequences, so that one bad submission does
    not block an analysis of fifty good ones - but tolerating is not the same as
    accepting, and `align()` will not let an odd sequence through silently. See
    `low_identity_rows`.
    """
    if len(sequences) < 2:
        raise AlignmentImpossible("need at least two sequences")

    odd = odd_prefixes(sequences, prefix=prefix)
    agreeing = len(sequences) - len(odd)

    if agreeing / len(sequences) < min_agreement:
        raise AlignmentImpossible(
            f"only {agreeing}/{len(sequences)} sequences share a common start. These may "
            "be differently rotated circular genomes or reverse-complemented; align them "
            "with MAFFT or MUSCLE instead"
        )

    lengths = [len(s) for s in sequences]
    spread = (max(lengths) - min(lengths)) / max(statistics.median(lengths), 1)
    if spread > 0.15:
        raise AlignmentImpossible(
            f"lengths differ by {spread:.0%} of the median, which exceeds what a banded "
            "aligner can be trusted with; use MAFFT or MUSCLE"
        )


def low_identity_rows(
    aligned: Sequence[str], centre: int, *, min_identity: float = 0.7
) -> list[tuple[int, float]]:
    """Aligned rows too unlike the centre to be a real alignment, as (index, identity).

    This is the check that catches what the prefix guard cannot. Two rotations of one
    genome, or a genome and its reverse complement, align at roughly the rate two random
    sequences do - about 25% of columns matching outright, and near 50% once the aligner
    has spent gaps buying agreement. Anything under 70% against the centre is not a
    poorly-conserved isolate, it is a sequence that should not have been in this
    alignment, and every variant its row contributes downstream is an artefact.
    """
    reference = aligned[centre]
    out = []
    for index, row in enumerate(aligned):
        if index == centre:
            continue
        score = _column_identity(row, reference)
        if score < min_identity:
            out.append((index, round(score, 4)))
    return out


def _column_identity(a: str, b: str) -> float:
    """Identity over columns where both rows have a called base."""
    compared = same = 0
    for x, y in zip(a, b, strict=False):
        if x not in "ACGT" or y not in "ACGT":
            continue
        compared += 1
        if x == y:
            same += 1
    return same / compared if compared else 0.0


def _identity(a: str, b: str) -> float:
    if not a or not b:
        return 0.0
    n = min(len(a), len(b))
    return sum(1 for i in range(n) if a[i] == b[i]) / n


def align_pair(
    a: str, b: str, *, band: int = 120, scoring: Scoring | None = None
) -> tuple[str, str]:
    """Banded global alignment of two sequences.

    The band is measured around the main diagonal, offset for any length difference, so
    a systematic insertion in one sequence does not push the optimal path out of the
    band on every subsequent row.
    """
    scoring = scoring or Scoring()
    a, b = a.upper(), b.upper()
    n, m = len(a), len(b)

    if not n or not m:
        return (a + GAP * m, GAP * n + b)

    # Widen the band to cover the length difference plus the requested slack.
    band = max(band, abs(n - m) + 16)
    drift = (m - n) / n if n else 0.0
    negative_infinity = float("-inf")

    # One row at a time, keeping the traceback for the whole matrix.
    previous = {0: 0.0}
    traceback: list[dict[int, str]] = [{}]

    for i in range(1, n + 1):
        centre = int(i * (1 + drift))
        low = max(0, centre - band)
        high = min(m, centre + band)

        current: dict[int, float] = {}
        row_back: dict[int, str] = {}

        for j in range(low, high + 1):
            if i == 0 and j == 0:
                current[j] = 0.0
                continue

            best = negative_infinity
            move = ""

            diagonal = previous.get(j - 1)
            if diagonal is not None and j >= 1:
                score = diagonal + scoring.pair(a[i - 1], b[j - 1])
                if score > best:
                    best, move = score, "D"

            up = previous.get(j)
            if up is not None:
                score = up + scoring.gap
                if score > best:
                    best, move = score, "U"

            left = current.get(j - 1)
            if left is not None:
                score = left + scoring.gap
                if score > best:
                    best, move = score, "L"

            if move:
                current[j] = best
                row_back[j] = move

        if not current:
            raise AlignmentImpossible(
                f"band of {band} was too narrow at row {i}; widen it or use a real aligner"
            )

        traceback.append(row_back)
        previous = current

    # Walk back from the corner. If the band excluded the true corner, start from the
    # best reachable cell on the last row rather than failing.
    i, j = n, m
    if j not in previous:
        j = max(previous, key=lambda k: previous[k])

    top: list[str] = []
    bottom: list[str] = []

    while i > 0 or j > 0:
        move = traceback[i].get(j) if i < len(traceback) else None
        if i == 0:
            top.append(GAP)
            bottom.append(b[j - 1])
            j -= 1
        elif j == 0 or move == "U":
            top.append(a[i - 1])
            bottom.append(GAP)
            i -= 1
        elif move == "L":
            top.append(GAP)
            bottom.append(b[j - 1])
            j -= 1
        elif move == "D":
            top.append(a[i - 1])
            bottom.append(b[j - 1])
            i -= 1
            j -= 1
        else:
            # Outside the band. Consume the longer side and continue.
            if i >= j:
                top.append(a[i - 1])
                bottom.append(GAP)
                i -= 1
            else:
                top.append(GAP)
                bottom.append(b[j - 1])
                j -= 1

    return "".join(reversed(top)), "".join(reversed(bottom))


def choose_centre(sequences: Sequence[str]) -> int:
    """Index of the sequence to align everything against.

    The one closest to the median length. A centre that is unusually short or long forces
    gaps into every other sequence and inflates the apparent indel count across the whole
    alignment.
    """
    median = statistics.median(len(s) for s in sequences)
    return min(range(len(sequences)), key=lambda i: abs(len(sequences[i]) - median))


def align(
    sequences: Sequence[str],
    *,
    band: int = 120,
    centre: int | None = None,
    check: bool = True,
    names: Sequence[str] | None = None,
    min_identity: float = 0.7,
) -> list[str]:
    """Centre-star multiple alignment. Returns sequences of equal length.

    Each sequence is aligned to the centre, then every gap any pairwise alignment
    introduced *into the centre* is propagated to all the others — which is what turns a
    set of independent pairwise alignments into one consistent multiple alignment.

    After merging, every row is compared with the centre and anything below
    `min_identity` raises, naming the sequence. A rotated or reverse-complemented
    submission passes the prefix guard whenever it is a small enough minority, and then
    aligns to noise: this is where it is caught. Pass `min_identity=0.0` to align
    anyway, having been told which rows are junk.
    """
    if len(sequences) < 2:
        raise AlignmentImpossible("need at least two sequences")
    if names is not None and len(names) != len(sequences):
        raise AlignmentImpossible(f"got {len(names)} names for {len(sequences)} sequences")
    if check:
        check_comparable(sequences, names=names)

    sequences = [s.upper() for s in sequences]
    centre = choose_centre(sequences) if centre is None else centre
    reference = sequences[centre]

    # For each sequence: where gaps must be inserted into the centre, and the partner row.
    pairwise: list[tuple[str, str]] = []
    for index, sequence in enumerate(sequences):
        if index == centre:
            pairwise.append((reference, reference))
        else:
            pairwise.append(align_pair(reference, sequence, band=band))

    # The merged centre needs, at every position, the maximum number of gaps any pairwise
    # alignment inserted there.
    gaps_after: list[int] = [0] * (len(reference) + 1)
    for aligned_reference, _ in pairwise:
        position = 0
        run = 0
        for character in aligned_reference:
            if character == GAP:
                run += 1
            else:
                gaps_after[position] = max(gaps_after[position], run)
                run = 0
                position += 1
        gaps_after[position] = max(gaps_after[position], run)

    merged: list[str] = []
    for aligned_reference, aligned_other in pairwise:
        row: list[str] = []
        position = 0  # index into the ungapped centre
        run: list[str] = []  # characters opposite a gap in the centre

        for reference_char, other_char in zip(aligned_reference, aligned_other, strict=True):
            if reference_char == GAP:
                run.append(other_char)
                continue
            # Flush any insertion that preceded this centre column, padded to the width
            # the widest insertion needs.
            row.append("".join(run).ljust(gaps_after[position], GAP))
            run = []
            row.append(other_char)
            position += 1

        row.append("".join(run).ljust(gaps_after[position], GAP))
        merged.append("".join(row))

    width = max(len(row) for row in merged)
    merged = [row.ljust(width, GAP) for row in merged]

    if min_identity > 0:
        bad = low_identity_rows(merged, centre, min_identity=min_identity)
        if bad:
            raise AlignmentImpossible(
                _low_identity_message(bad, sequences, centre, names, min_identity)
            )

    return merged


def _low_identity_message(
    bad: list[tuple[int, float]],
    sequences: Sequence[str],
    centre: int,
    names: Sequence[str] | None,
    min_identity: float,
) -> str:
    """Name the offending sequences, and say whether reverse-complementing would fix it.

    Checking the reverse complement costs one pairwise alignment per bad row and turns
    "this failed" into "this is backwards", which is the difference between a user
    giving up and a user fixing their input in one step.
    """

    def label(index: int) -> str:
        return names[index] if names else f"sequence {index}"

    lines = []
    for index, score in bad[:5]:
        flipped = align_pair(sequences[centre], reverse_complement(sequences[index]))
        flipped_score = _column_identity(flipped[1], flipped[0])
        hint = (
            f"; its reverse complement aligns at {flipped_score:.0%}, so it is probably "
            "on the opposite strand"
            if flipped_score > score + 0.15
            else "; it may be a differently rotated circular genome, or a different species"
        )
        lines.append(f"  {label(index)}: {score:.0%} identity to the centre{hint}")

    more = f"\n  ... and {len(bad) - 5} more" if len(bad) > 5 else ""
    return (
        f"{len(bad)} of {len(sequences)} sequences aligned at less than "
        f"{min_identity:.0%} identity to the centre ({label(centre)}), which is close "
        "to what two unrelated sequences score. "
        "Every variant called from those rows would be an artefact:\n"
        + "\n".join(lines)
        + more
        + "\nFix or drop them, or pass min_identity=0.0 to align them anyway."
    )


def alignment_report(aligned: Sequence[str]) -> dict:
    """Quality signals for an alignment, so it can be judged rather than trusted."""
    if not aligned:
        return {"sequences": 0}

    width = len(aligned[0])
    gap_columns = 0
    identical_columns = 0
    gaps = 0

    for i in range(width):
        column = [s[i] for s in aligned]
        called = [c for c in column if c in "ACGT"]
        gaps += len(column) - len(called)
        if not called:
            gap_columns += 1
        elif len(set(called)) == 1 and len(called) == len(column):
            identical_columns += 1

    return {
        "sequences": len(aligned),
        "width": width,
        "all_gap_columns": gap_columns,
        "invariant_columns": identical_columns,
        "invariant_fraction": round(identical_columns / width, 4) if width else 0.0,
        "gap_fraction": round(gaps / (width * len(aligned)), 4) if width else 0.0,
    }
