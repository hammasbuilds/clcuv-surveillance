"""Tests for location normalisation and for the thresholds the headline rests on.

Both were found by asking the repository's own question back at it: the stratified
control claims to compare like with like, and the headline claims a number. Neither
survives contact with the raw `/country` strings or with a sweep of `min_samples`
unless the code below does its job.
"""

from __future__ import annotations

import pytest

from clcuv.atlas import (
    Isolate,
    build_atlas,
    check_period,
    emergence_sensitivity,
    emerging_variants,
    group_linked,
    rises_within_locations,
)
from clcuv.geo import normalise_isolates, normalise_location, resolved_strata

REFERENCE = "ACGT" * 25  # 100 bp


def _with_variant(position: int, base: str) -> str:
    chars = list(REFERENCE)
    chars[position] = base
    return "".join(chars)


# --- normalisation --------------------------------------------------------


@pytest.mark.parametrize(
    ("raw", "stratum"),
    [
        # Six spellings of one province, all seen in the committed corpus.
        ("Pakistan: Punjab", "Pakistan: Punjab"),
        ("Pakistan: Punjab province", "Pakistan: Punjab"),
        ("Pakistan: Punjab,Bahawalpur", "Pakistan: Punjab"),
        ("Pakistan: Vehari, Punjab", "Pakistan: Punjab"),
        ("Pakistan: Vehari-Punjab", "Pakistan: Punjab"),
        ("Pakistan:Faisalabad", "Pakistan: Punjab"),
        # Districts resolve to the division that contains them.
        ("Pakistan: Sakrand", "Pakistan: Sindh"),
        ("Pakistan: Sindh,Tando Jam", "Pakistan: Sindh"),
        ("China: GD, Enping", "China: Guangdong"),
        ("China: Xinjiang-Shihezi", "China: Xinjiang"),
        ("India: Hissar, Haryana", "India: Haryana"),
        ("India: Sri-Ganganagar", "India: Rajasthan"),
        ("India:Pusa,New Dehli", "India: Delhi"),
    ],
)
def test_spelling_variants_collapse_to_one_stratum(raw, stratum):
    assert normalise_location(raw).stratum == stratum


def test_the_same_province_name_in_two_countries_stays_two_strata():
    """Indian Punjab and Pakistani Punjab are not one place, and merging them on the
    bare name would be a worse error than not normalising at all."""
    assert normalise_location("India: Ludhiana, Punjob").stratum == "India: Punjab"
    assert normalise_location("Pakistan: Punjab").stratum == "Pakistan: Punjab"


def test_a_country_only_record_is_kept_but_unresolved():
    place = normalise_location("Pakistan")
    assert place.stratum == "Pakistan"
    assert place.resolved is False
    assert "country-only" in place.flags


def test_a_named_place_we_have_no_entry_for_is_unresolved_not_guessed():
    place = normalise_location("Pakistan: Somewhere Nobody Mapped")
    assert place.resolved is False
    assert "unmapped-detail" in place.flags


def test_a_division_in_the_wrong_country_is_flagged_not_moved():
    """Eleven records are filed as Pakistan: Rajasthan, and Rajasthan is in India.
    Reassigning them would be this module inventing metadata."""
    place = normalise_location("Pakistan: Rajesthan")
    assert place.country == "Pakistan"
    assert place.admin1 == "Rajasthan"
    assert place.flags == ("admin1-not-in-country",)
    assert normalise_location("Pakistan: Rajasthan").stratum == place.stratum


def test_empty_and_missing_locations_do_not_raise():
    assert normalise_location("").stratum == "unknown"
    assert normalise_location("   ").flags == ("empty",)


def test_normalise_isolates_rewrites_locations_and_reports_the_mapping():
    isolates = [
        Isolate("a", REFERENCE, period="2019", location="Pakistan: Vehari-Punjab"),
        Isolate("b", REFERENCE, period="2019", location="Pakistan: Punjab province"),
        Isolate("c", REFERENCE, period="2019", location="Pakistan"),
    ]
    rewritten, mapping = normalise_isolates(isolates)
    assert [i.location for i in rewritten] == [
        "Pakistan: Punjab",
        "Pakistan: Punjab",
        "Pakistan",
    ]
    assert resolved_strata(mapping) == {"Pakistan: Punjab"}


# --- what normalisation changes about the stratified test -----------------


def _split_province_confounder():
    """One province spelled two ways, sampled in two years.

    The variant rises within Punjab. Unnormalised, each spelling holds half the
    genomes, no stratum reaches min_samples, and the control silently confirms
    nothing - which reads as "not confirmed" and is really "not tested".
    """
    isolates = []
    for i in range(10):
        spelling = "Pakistan: Punjab" if i % 2 else "Pakistan: Vehari-Punjab"
        isolates.append(Isolate(f"a{i}", REFERENCE, period="2019", location=spelling))
    for i in range(10):
        spelling = "Pakistan: Punjab" if i % 2 else "Pakistan: Punjab province"
        isolates.append(Isolate(f"b{i}", _with_variant(50, "A"), period="2021", location=spelling))
    return isolates


def test_unnormalised_spellings_split_a_province_below_the_testable_size():
    atlas = build_atlas(_split_province_confounder(), REFERENCE)
    variant = next(v for v in atlas if v.label == "G51A")
    assert rises_within_locations(variant, min_samples=8) == []


def test_normalising_restores_the_comparison():
    isolates, _ = normalise_isolates(_split_province_confounder())
    atlas = build_atlas(isolates, REFERENCE)
    variant = next(v for v in atlas if v.label == "G51A")
    assert rises_within_locations(variant, min_samples=8) == ["Pakistan: Punjab"]


def test_a_country_only_stratum_cannot_confirm_a_within_place_rise():
    """`Pakistan` with no province mixes provinces, so confirming "within Pakistan"
    reintroduces the confounder stratification exists to remove."""
    isolates = [Isolate(f"a{i}", REFERENCE, period="2019", location="Pakistan") for i in range(10)]
    isolates += [
        Isolate(f"b{i}", _with_variant(50, "A"), period="2021", location="Pakistan")
        for i in range(10)
    ]
    rewritten, mapping = normalise_isolates(isolates)
    variant = next(v for v in build_atlas(rewritten, REFERENCE) if v.label == "G51A")

    assert rises_within_locations(variant, min_samples=8) == ["Pakistan"]
    assert (
        rises_within_locations(variant, min_samples=8, eligible_locations=resolved_strata(mapping))
        == []
    )


def test_eligible_locations_flows_through_emerging_variants():
    isolates = [Isolate(f"a{i}", REFERENCE, period="2019", location="Pakistan") for i in range(10)]
    isolates += [
        Isolate(f"b{i}", _with_variant(50, "A"), period="2021", location="Pakistan")
        for i in range(10)
    ]
    atlas = build_atlas(isolates, REFERENCE)
    assert emerging_variants(atlas, min_samples=8, stratify=True)
    assert not emerging_variants(
        atlas, min_samples=8, stratify=True, eligible_locations={"Pakistan: Punjab"}
    )


# --- periods --------------------------------------------------------------


@pytest.mark.parametrize("period", ["2019", "2019-05", "2019-05-01", "2019-Q3", ""])
def test_sortable_periods_are_accepted(period):
    check_period(period)
    assert Isolate("a", REFERENCE, period=period).period == period


@pytest.mark.parametrize("period", ["May-2025", "Jan-2026", "p", "19", "2019-13", "2019-Q5"])
def test_unsortable_periods_are_refused(period):
    with pytest.raises(ValueError, match="sortable form"):
        Isolate("a", REFERENCE, period=period)


def test_the_trap_the_period_check_exists_for():
    """Lexical order on month names reverses the trend, and nothing errors."""
    assert sorted(["May-2025", "Jan-2026"]) == ["Jan-2026", "May-2025"]
    assert sorted(["2025-05", "2026-01"]) == ["2025-05", "2026-01"]


# --- threshold sensitivity ------------------------------------------------


def _one_genome_from_the_threshold():
    """2019 has exactly 8 genomes, 2020 has 12, 2021 has 12.

    The variant is flat between 2020 and 2021 and high in 2019, so whether anything
    looks like it is rising depends entirely on whether 2019 is an endpoint - which
    `min_samples=8` decides, and `min_samples=9` decides the other way.
    """
    isolates = [Isolate(f"a{i}", REFERENCE, period="2019", location="P") for i in range(8)]
    isolates += [
        Isolate(
            f"b{i}", _with_variant(50, "A") if i < 8 else REFERENCE, period="2020", location="P"
        )
        for i in range(12)
    ]
    isolates += [
        Isolate(
            f"c{i}", _with_variant(50, "A") if i < 9 else REFERENCE, period="2021", location="P"
        )
        for i in range(12)
    ]
    return build_atlas(isolates, REFERENCE)


def test_the_headline_moves_with_min_samples():
    atlas = _one_genome_from_the_threshold()
    assert [e.variant.label for e in emerging_variants(atlas, min_samples=8)] == ["G51A"]
    assert emerging_variants(atlas, min_samples=9) == []


def test_the_sweep_shows_that_movement_rather_than_hiding_it():
    rows = emergence_sensitivity(_one_genome_from_the_threshold(), min_samples_values=(8, 9))
    assert [(r["min_samples"], r["pooled"]) for r in rows] == [(8, 1), (9, 0)]


def test_periods_used_says_how_thin_the_comparison_was():
    found = emerging_variants(_one_genome_from_the_threshold(), min_samples=8)
    assert found[0].periods_used == 3


def _a_spike_in_the_middle():
    """The variant sweeps in 2020 and is gone again by 2021.

    Comparing the first and last period reports nothing happened, which is true of the
    endpoints and false of the data.
    """
    isolates = [Isolate(f"a{i}", REFERENCE, period="2019", location="P") for i in range(12)]
    isolates += [
        Isolate(
            f"b{i}", _with_variant(50, "A") if i < 9 else REFERENCE, period="2020", location="P"
        )
        for i in range(12)
    ]
    isolates += [Isolate(f"c{i}", REFERENCE, period="2021", location="P") for i in range(12)]
    return build_atlas(isolates, REFERENCE)


def test_comparing_only_the_endpoints_misses_a_rise_between_them():
    atlas = _a_spike_in_the_middle()
    assert emerging_variants(atlas, min_samples=8) == []
    assert [
        e.variant.label for e in emerging_variants(atlas, min_samples=8, compare="adjacent")
    ] == ["G51A"]


def test_linked_variants_are_one_event_not_many():
    """Two mutations carried by the same genomes rise together and are one finding."""
    isolates = [Isolate(f"a{i}", REFERENCE, period="2019", location="P") for i in range(12)]
    both = list(REFERENCE)
    both[50], both[60] = "A", "T"
    isolates += [Isolate(f"b{i}", "".join(both), period="2021", location="P") for i in range(12)]

    found = emerging_variants(build_atlas(isolates, REFERENCE), min_samples=8)
    assert sorted(e.variant.label for e in found) == ["A61T", "G51A"]
    assert len(group_linked(found)) == 1


def test_unlinked_variants_stay_separate():
    isolates = [Isolate(f"a{i}", REFERENCE, period="2019", location="P") for i in range(12)]
    for i in range(12):
        chars = list(REFERENCE)
        chars[50] = "A"
        if i < 6:  # carried by half the 2021 genomes, so a different trajectory
            chars[60] = "T"
        isolates.append(Isolate(f"b{i}", "".join(chars), period="2021", location="P"))

    found = emerging_variants(build_atlas(isolates, REFERENCE), min_samples=8)
    assert len(found) == 2
    assert len(group_linked(found)) == 2


def test_an_unknown_comparison_mode_is_refused():
    with pytest.raises(ValueError, match="compare must be one of"):
        emerging_variants(_one_genome_from_the_threshold(), compare="whatever")
