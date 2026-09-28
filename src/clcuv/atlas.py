"""Mutation atlas: what is changing, where, and whether it is spreading.

A mutation atlas is not a list of mutations. Every genome differs from the reference at
dozens of positions and almost none of them matter. The question surveillance exists to
answer is narrower:

    **which variant is rising in frequency, and how fast?**

A mutation at 40% frequency that has sat at 40% for three seasons is background. A
mutation that went 2% → 8% → 25% in three seasons is a lineage winning, and it is worth
knowing about while it is still at 25%.

So the atlas is organised around **trajectories**, not snapshots, and the alerting is on
the derivative rather than the level.
"""

from __future__ import annotations

import math
import re
from collections import Counter, defaultdict
from collections.abc import Collection, Sequence
from dataclasses import dataclass, field

# Periods are ordered by sorting the string, so only forms that sort chronologically
# are accepted. `May-2025` and `Jan-2026` do not: sorted() puts January 2026 first,
# `emerging_variants` then reads the trend backwards, and nothing errors. Rejecting the
# format is the only version of this that cannot fail silently.
PERIOD_PATTERN = re.compile(r"^\d{4}(-(0[1-9]|1[0-2])(-\d{2})?|-Q[1-4])?$")
PERIOD_FORMS = "YYYY, YYYY-MM, YYYY-MM-DD or YYYY-Qn"


def check_period(period: str) -> None:
    """Raise unless `period` sorts chronologically as a string."""
    if period and not PERIOD_PATTERN.match(period):
        raise ValueError(
            f"period {period!r} is not a sortable form. Periods are ordered by sorting "
            f"the string, so use {PERIOD_FORMS} - a form like 'May-2025' sorts before "
            "'Jan-2026' and would reverse the trend silently. Pass period='' for unknown."
        )


@dataclass(frozen=True)
class Isolate:
    """One sequenced sample, with the metadata that makes it surveillance."""

    name: str
    sequence: str
    # Collection period - a season, month or year. Ordered lexically, so the format is
    # checked: see `check_period`.
    period: str = ""
    location: str = ""
    host: str = ""

    def __post_init__(self) -> None:
        check_period(self.period)


@dataclass
class Variant:
    position: int  # 0-based alignment column
    reference: str
    alternate: str
    counts_by_period: dict[str, int] = field(default_factory=dict)
    totals_by_period: dict[str, int] = field(default_factory=dict)
    locations: Counter = field(default_factory=Counter)
    # Cross-tabulated by (period, location). Needed because a frequency that rises
    # between two periods may only reflect a change in *where* the samples came from -
    # see `emerging_variants(stratify=True)`.
    counts_by_stratum: dict[tuple[str, str], int] = field(default_factory=dict)
    totals_by_stratum: dict[tuple[str, str], int] = field(default_factory=dict)

    @property
    def label(self) -> str:
        return f"{self.reference}{self.position + 1}{self.alternate}"

    def frequency(self, period: str) -> float:
        total = self.totals_by_period.get(period, 0)
        return round(self.counts_by_period.get(period, 0) / total, 6) if total else 0.0

    def strata(self) -> list[str]:
        """Locations with data for this position."""
        return sorted({location for _, location in self.totals_by_stratum})

    def trajectory_within(self, location: str) -> list[tuple[str, int, int]]:
        """(period, carrying, sampled) for one location, in period order.

        Counts rather than frequencies, because a stratified test needs the
        denominators and a frequency has already thrown them away.
        """
        return [
            (period, self.counts_by_stratum.get((period, location), 0), total)
            for (period, loc), total in sorted(self.totals_by_stratum.items())
            if loc == location
        ]

    def trajectory(self) -> list[tuple[str, float, int]]:
        """(period, frequency, sample size), in period order."""
        return [
            (p, self.frequency(p), self.totals_by_period.get(p, 0))
            for p in sorted(self.totals_by_period)
        ]

    @property
    def overall_frequency(self) -> float:
        total = sum(self.totals_by_period.values())
        return round(sum(self.counts_by_period.values()) / total, 6) if total else 0.0


def build_atlas(
    isolates: Sequence[Isolate], reference: str, *, min_frequency: float = 0.01
) -> list[Variant]:
    """Every variant against the reference, with its trajectory.

    `min_frequency` drops singletons. In a set of a few hundred genomes, a variant seen
    once is more likely a sequencing error than a lineage, and an atlas full of
    sequencing errors is one nobody reads.
    """
    if not isolates:
        return []

    length = len(reference)
    for isolate in isolates:
        if len(isolate.sequence) != length:
            raise ValueError(
                f"{isolate.name}: length {len(isolate.sequence)} does not match the "
                f"reference ({length}) - sequences must be aligned"
            )

    variants: dict[tuple[int, str], Variant] = {}
    period_totals_by_position: dict[int, Counter] = defaultdict(Counter)
    stratum_totals_by_position: dict[int, Counter] = defaultdict(Counter)

    # Denominators are counted per position, not per isolate: a genome with an N at a
    # position contributes no information there and must not inflate the denominator.
    for isolate in isolates:
        for i, base in enumerate(isolate.sequence.upper()):
            if base in "ACGT":
                period_totals_by_position[i][isolate.period] += 1
                stratum_totals_by_position[i][(isolate.period, isolate.location)] += 1

    for isolate in isolates:
        for i, base in enumerate(isolate.sequence.upper()):
            ref = reference[i].upper()
            if base not in "ACGT" or ref not in "ACGT" or base == ref:
                continue
            key = (i, base)
            variant = variants.setdefault(key, Variant(position=i, reference=ref, alternate=base))
            variant.counts_by_period[isolate.period] = (
                variant.counts_by_period.get(isolate.period, 0) + 1
            )
            stratum = (isolate.period, isolate.location)
            variant.counts_by_stratum[stratum] = variant.counts_by_stratum.get(stratum, 0) + 1
            variant.locations[isolate.location] += 1

    for (position, _), variant in variants.items():
        variant.totals_by_period = dict(period_totals_by_position[position])
        variant.totals_by_stratum = dict(stratum_totals_by_position[position])

    return sorted(
        (v for v in variants.values() if v.overall_frequency >= min_frequency),
        key=lambda v: (-v.overall_frequency, v.position),
    )


def two_proportion_z(successes_a: int, total_a: int, successes_b: int, total_b: int) -> float:
    """Z-statistic for a difference between two observed proportions.

    Necessary because frequency estimates are noisy and the noise scales with sample
    size. Two seasons of 100 genomes each can easily differ by ten points with nothing
    happening at all, and a detector that fires on that is one nobody trusts by the
    third false alarm.
    """
    if total_a <= 0 or total_b <= 0:
        return 0.0
    pooled = (successes_a + successes_b) / (total_a + total_b)
    if pooled in (0.0, 1.0):
        return 0.0
    standard_error = math.sqrt(pooled * (1 - pooled) * (1 / total_a + 1 / total_b))
    if standard_error == 0:
        return 0.0
    return (successes_b / total_b - successes_a / total_a) / standard_error


@dataclass
class Emerging:
    variant: Variant
    first_period: str
    last_period: str
    first_frequency: float
    last_frequency: float
    change: float
    fold: float | None
    z: float = 0.0
    # Locations where the rise holds within that location alone. Empty when the trend
    # only exists in the pooled data, which usually means it is a sampling artefact.
    confirmed_in: list[str] = field(default_factory=list)
    # How many periods had enough sequences to be usable. Two means the comparison rests
    # on the only two periods available, and one genome moving across `min_samples`
    # would have changed which two those are.
    periods_used: int = 2

    @property
    def stratified(self) -> bool:
        return bool(self.confirmed_in)

    def summary(self) -> dict:
        return {
            "variant": self.variant.label,
            "from": f"{self.first_frequency:.1%} ({self.first_period})",
            "to": f"{self.last_frequency:.1%} ({self.last_period})",
            "change": round(self.change, 4),
            "fold": self.fold,
            "z": round(self.z, 3),
            "periods_used": self.periods_used,
            "confirmed_in": self.confirmed_in,
            "locations": dict(self.variant.locations.most_common(5)),
        }


COMPARISONS = ("extremes", "adjacent")


def _candidate_pairs(usable: list, compare: str) -> list[tuple]:
    """Which (earlier, later) period pairs a rise may be claimed between.

    `extremes` compares the first and last usable period. It is the natural reading of
    "is this rising", and it has a failure mode worth naming: the endpoints are chosen
    by `min_samples`, so one genome crossing that threshold swaps an endpoint and can
    change the answer. `adjacent` compares consecutive usable periods instead, which
    does not depend on which period happens to be the earliest, and fires if any single
    step qualifies. Neither is the true one. Running both is how you find out whether a
    result is a finding or an artefact of the setting - see `emergence_sensitivity`.
    """
    if compare == "extremes":
        return [(usable[0], usable[-1])] if len(usable) >= 2 else []
    if compare == "adjacent":
        return list(zip(usable, usable[1:], strict=False))
    raise ValueError(f"compare must be one of {COMPARISONS}, not {compare!r}")


def rises_within_locations(
    variant: Variant,
    *,
    min_change: float = 0.10,
    min_samples: int = 8,
    min_z: float = 1.96,
    compare: str = "extremes",
    eligible_locations: Collection[str] | None = None,
) -> list[str]:
    """Locations where this variant rises *within that location's own samples*.

    The confounder this exists for, found on real GenBank data rather than imagined:
    the set of places sampled changes between periods, so a variant common in one
    province and absent from another appears to "emerge" when the sampling programme
    moves, with a perfectly valid z-statistic, having done nothing at all. The frequency
    genuinely rose; the population being sampled is simply not the same population.

    `eligible_locations` restricts which strata may confirm. This matters because a
    location string is not automatically a location: a record filed as `Pakistan` with
    no province mixes provinces together, so confirming "within Pakistan" reintroduces
    exactly the confounder being controlled for. `clcuv.geo.resolved_strata` supplies
    the set that names a first-level division. `None` means every stratum may confirm,
    which is only safe when the locations are already known to be comparable.
    """
    confirmed: list[str] = []

    for location in variant.strata():
        if eligible_locations is not None and location not in eligible_locations:
            continue

        usable = [
            (period, carrying, sampled)
            for period, carrying, sampled in variant.trajectory_within(location)
            if sampled >= min_samples
        ]

        for (_, first_count, first_total), (_, last_count, last_total) in _candidate_pairs(
            usable, compare
        ):
            change = last_count / last_total - first_count / first_total
            if change < min_change:
                continue
            if two_proportion_z(first_count, first_total, last_count, last_total) < min_z:
                continue
            confirmed.append(location)
            break

    return confirmed


def emerging_variants(
    variants: Sequence[Variant],
    *,
    min_change: float = 0.10,
    min_samples: int = 10,
    min_z: float = 1.96,
    stratify: bool = False,
    compare: str = "extremes",
    eligible_locations: Collection[str] | None = None,
) -> list[Emerging]:
    """Variants whose frequency is rising, and rising by more than chance.

    Two filters, both necessary and neither sufficient.

    `min_samples` excludes periods with too few sequences. A variant at "100%" in a
    period where three genomes were sequenced is not at 100% — it is unmeasured.

    `min_z` requires the rise to be statistically distinguishable from sampling noise.
    Without it the detector fires on the ordinary wobble of two finite samples: this
    was a real false positive here, where a variant sitting at a constant 40% was
    reported as rising because two seasons of 100 genomes happened to land at 33% and
    45%. Effect size alone is not evidence.

    Every one of those thresholds is a judgement call, and on a corpus this size the
    answer moves when they move. Do not read one call of this function as the result:
    run `emergence_sensitivity` and read the shape of the sweep.
    """
    out: list[Emerging] = []

    for variant in variants:
        usable = [(period, freq, n) for period, freq, n in variant.trajectory() if n >= min_samples]

        best: Emerging | None = None
        for (first_period, first_freq, first_n), (
            last_period,
            last_freq,
            last_n,
        ) in _candidate_pairs(usable, compare):
            change = last_freq - first_freq
            if change < min_change:
                continue

            z = two_proportion_z(
                variant.counts_by_period.get(first_period, 0),
                first_n,
                variant.counts_by_period.get(last_period, 0),
                last_n,
            )
            if z < min_z:
                continue

            if best is not None and change <= best.change:
                continue

            best = Emerging(
                variant=variant,
                first_period=first_period,
                last_period=last_period,
                first_frequency=first_freq,
                last_frequency=last_freq,
                change=change,
                # Fold change is undefined from zero. Reporting it as infinite, or as a
                # large number, turns "newly detected" into "exploding" — a different claim.
                fold=round(last_freq / first_freq, 3) if first_freq > 0 else None,
                z=round(z, 4),
                periods_used=len(usable),
            )

        if best is None:
            continue

        best.confirmed_in = rises_within_locations(
            variant,
            min_change=min_change,
            min_samples=min_samples,
            min_z=min_z,
            compare=compare,
            eligible_locations=eligible_locations,
        )
        # With `stratify`, a rise that only exists in the pooled data is discarded. The
        # pooled test is not wrong - the frequency really did rise - but it cannot tell a
        # lineage spreading from a sampling programme moving to a different province.
        if stratify and not best.confirmed_in:
            continue

        out.append(best)

    # Variants confirmed within a location first: those are the ones that survived the
    # question "compared with what?"
    return sorted(out, key=lambda e: (not e.stratified, -e.change))


def group_linked(emerging: Sequence[Emerging]) -> list[list[Emerging]]:
    """Group variants that rise and fall together, because they are one event.

    Mutations on the same lineage travel together. If a lineage carrying seventy-six
    substitutions spreads, all seventy-six rise in lockstep, and a report saying
    "seventy-six variants are emerging" has described one event seventy-six times -
    the same error as counting one clonal batch as eight observations, one level up.

    Variants are grouped when their per-period carrier counts are identical, which is
    exact linkage over the periods observed rather than a correlation threshold. The
    number worth reporting is `len(group_linked(found))` beside `len(found)`: on this
    corpus those are 3 and 76.
    """
    groups: dict[tuple, list[Emerging]] = {}
    for item in emerging:
        signature = tuple(sorted(item.variant.counts_by_period.items()))
        groups.setdefault(signature, []).append(item)
    return sorted(groups.values(), key=lambda g: (-len(g), g[0].variant.position))


def emergence_sensitivity(
    variants: Sequence[Variant],
    *,
    min_samples_values: Sequence[int] = (5, 6, 7, 8, 9, 10),
    min_change: float = 0.10,
    min_z: float = 1.96,
    compare: str = "extremes",
    eligible_locations: Collection[str] | None = None,
) -> list[dict]:
    """How many variants survive each control, across a range of `min_samples`.

    `min_samples` is the least defensible number in this module. It is not estimated
    from anything; it is a line drawn to keep periods with too few genomes out of the
    test, and where it is drawn decides which periods become the endpoints of the
    comparison. On a corpus where whole years hold eight or nine genomes, moving it by
    one moves the answer.

    A single count reported without this sweep is not a finding, it is a setting. The
    sweep is cheap, it is the first thing a reviewer asks for, and if the count is
    stable across it that is worth far more than the count itself.
    """
    rows = []
    for min_samples in min_samples_values:
        pooled = emerging_variants(
            variants,
            min_change=min_change,
            min_samples=min_samples,
            min_z=min_z,
            compare=compare,
        )
        stratified = [
            e
            for e in pooled
            if rises_within_locations(
                e.variant,
                min_change=min_change,
                min_samples=min_samples,
                min_z=min_z,
                compare=compare,
                eligible_locations=eligible_locations,
            )
        ]
        rows.append(
            {
                "min_samples": min_samples,
                "compare": compare,
                "pooled": len(pooled),
                "stratified": len(stratified),
            }
        )
    return rows


def geographic_spread(variant: Variant) -> dict:
    """Where a variant has been seen, and how concentrated it is.

    A variant confined to one district is a local lineage; the same variant across five
    districts is spreading, and the difference changes what anyone should do about it.
    """
    total = sum(variant.locations.values())
    if not total:
        return {"locations": 0, "concentration": 0.0}

    shares = [n / total for n in variant.locations.values()]
    return {
        "variant": variant.label,
        "locations": len(variant.locations),
        "top": dict(variant.locations.most_common(5)),
        # Herfindahl index: 1.0 means a single location, low means widely spread.
        "concentration": round(sum(s * s for s in shares), 6),
        "spreading": len(variant.locations) >= 3 and sum(s * s for s in shares) < 0.5,
    }
