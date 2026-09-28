"""CLCuV surveillance tests.

Phylogenetics and selection analysis are exact: an additive distance matrix has one
correct tree, and dN/dS on synonymous-only changes is zero. Those are asserted. The
surveillance layer is tested on constructed trajectories, so a rising variant is rising
by construction rather than by luck.
"""

from __future__ import annotations

import math
import random

import pytest

from clcuv.atlas import (
    Isolate,
    build_atlas,
    emerging_variants,
    geographic_spread,
    two_proportion_z,
)
from clcuv.classify import (
    StrainClassifier,
    cosine_distance,
    kmer_profile,
    recombination_signal,
)
from clcuv.codon import (
    CODON_TABLE,
    CodonError,
    amino_acid_changes,
    is_synonymous,
    selection_pressure,
    synonymous_sites,
    translate,
)
from clcuv.phylo import (
    PhyloError,
    clades,
    distance_matrix,
    jukes_cantor,
    neighbour_joining,
    p_distance,
    quote_newick,
    upgma,
)

REF = "ATGGCTAAGCGTCCAGGATCCAAGTTCGATCCGTTAACGGCTAAGGCATCGTAGCTAGCATCG"


class TestGeneticCode:
    @pytest.mark.parametrize(
        ("codon", "residue"),
        [
            ("ATG", "M"),
            ("TGG", "W"),
            ("TAA", "*"),
            ("TAG", "*"),
            ("TGA", "*"),
            ("GGG", "G"),
            ("TTT", "F"),
            ("CTA", "L"),
        ],
    )
    def test_codon_table(self, codon, residue):
        assert CODON_TABLE[codon] == residue

    def test_table_is_complete(self):
        assert len(CODON_TABLE) == 64

    def test_translation(self):
        assert translate("ATGGCTTGGTAA") == "MAW*"

    def test_translation_can_stop_at_a_stop_codon(self):
        assert translate("ATGGCTTGGTAA", to_stop=True) == "MAW"

    def test_a_partial_codon_is_dropped_not_padded(self):
        """Padding invents a residue that is not in the data."""
        assert translate("ATGGCTTG") == "MA"

    def test_an_ambiguous_codon_becomes_x(self):
        assert translate("ATGNNN") == "MX"


class TestSelection:
    def test_methionine_and_tryptophan_have_no_synonymous_sites(self):
        """They have exactly one codon each, so every substitution changes the
        protein."""
        assert synonymous_sites("ATG") == 0.0
        assert synonymous_sites("TGG") == 0.0

    def test_a_fourfold_degenerate_third_position_gives_one_site(self):
        assert synonymous_sites("GGG") == pytest.approx(1.0)

    def test_leucine_has_more_than_one_synonymous_site(self):
        """CTA is Leu; so is TTA, so the first position is partly synonymous too."""
        assert synonymous_sites("CTA") > 1.0

    def test_synonymous_detection(self):
        assert is_synonymous("CTA", "CTG")  # both Leu
        assert not is_synonymous("ATG", "ATA")  # Met -> Ile

    def test_synonymous_only_changes_give_a_ratio_of_zero(self):
        a = "ATGGCTTTAGGGCCCAAA"
        b = "ATGGCCTTGGGACCCAAG"
        result = selection_pressure(a, b)
        assert result.nonsynonymous_differences == 0
        assert result.dn_ds == 0.0
        assert "purifying" in result.interpretation()

    def test_no_synonymous_differences_is_undetermined_not_infinite(self):
        """Reporting infinity turns 'we cannot tell' into 'strong positive selection',
        which is the wrong direction to be wrong in for an alert."""
        result = selection_pressure("ATGGCTTTA", "ATGTCTTTA")
        assert result.dn_ds is None
        assert "undetermined" in result.interpretation()

    def test_identical_sequences_have_no_differences(self):
        result = selection_pressure(REF[:60], REF[:60])
        assert result.synonymous_differences == 0
        assert result.nonsynonymous_differences == 0

    def test_unaligned_sequences_are_refused(self):
        with pytest.raises(CodonError):
            selection_pressure("ATGGCT", "ATGGC")

    def test_amino_acid_changes_use_residue_coordinates(self):
        """That is how the literature names resistance-breaking mutations, and a
        report nobody can cross-reference is one nobody uses."""
        assert amino_acid_changes("ATGGCTTTA", "ATGTCTTTA") == [(2, "A", "S")]

    def test_multiple_substitutions_in_one_codon_are_averaged(self):
        """There is no phylogeny here to pick a mutational order, so averaging is the
        honest treatment rather than assuming one."""
        result = selection_pressure("ATGCTA", "ATGTTG")
        assert result.synonymous_differences + result.nonsynonymous_differences > 0


class TestDistance:
    def test_identical_sequences_are_zero_apart(self):
        assert p_distance(REF, REF) == 0.0

    def test_differences_are_counted(self):
        assert p_distance("ACGT", "ACGA") == 0.25

    def test_gaps_are_skipped_not_counted_as_differences(self):
        """A gap is missing data; counting it makes poorly sequenced strains look
        divergent."""
        assert p_distance("ACGT", "AC-T") == 0.0

    def test_unaligned_sequences_are_refused(self):
        with pytest.raises(PhyloError):
            p_distance("ACGT", "ACG")

    def test_jukes_cantor_exceeds_the_raw_proportion(self):
        """The same site can mutate twice, and the second change hides the first."""
        assert jukes_cantor(0.3) > 0.3

    def test_the_correction_is_negligible_at_low_divergence(self):
        assert jukes_cantor(0.01) == pytest.approx(0.01, abs=0.001)

    def test_saturation_is_infinite_not_a_large_number(self):
        """At 75% divergence the sequences are indistinguishable from random, and the
        honest statement is that the distance cannot be estimated."""
        assert math.isinf(jukes_cantor(0.75))

    def test_matrix_is_symmetric_with_a_zero_diagonal(self):
        names, matrix = distance_matrix([REF, REF[::-1], REF[10:] + REF[:10]])
        assert all(matrix[i][i] == 0.0 for i in range(3))
        assert all(matrix[i][j] == matrix[j][i] for i in range(3) for j in range(3))


class TestTreeBuilding:
    # True tree ((A:5,B:2):1,(C:2,D:3)) — exact additive distances, and A is a
    # fast-evolving lineage, which is what breaks the molecular clock.
    NAMES = ["A", "B", "C", "D"]
    MATRIX = [[0, 7, 8, 9], [7, 0, 5, 6], [8, 5, 0, 5], [9, 6, 5, 0]]
    TRUTH = {frozenset(["A", "B"]), frozenset(["C", "D"])}

    def test_neighbour_joining_recovers_the_true_topology(self):
        assert clades(neighbour_joining(self.NAMES, self.MATRIX)) >= self.TRUTH

    def test_neighbour_joining_recovers_the_true_branch_lengths(self):
        tree = neighbour_joining(self.NAMES, self.MATRIX)
        lengths = {}

        def walk(node):
            for child, length in node.children:
                if child.is_leaf:
                    lengths[child.name] = length
                walk(child)

        walk(tree)
        assert lengths == pytest.approx({"A": 5.0, "B": 2.0, "C": 2.0, "D": 3.0})

    def test_upgma_gets_it_wrong_when_rates_differ(self):
        """UPGMA assumes a molecular clock. A strain under resistance pressure
        evolves faster — which is exactly the strain surveillance cares about — and
        UPGMA places it wrongly, confidently."""
        assert not clades(upgma(self.NAMES, self.MATRIX)) >= self.TRUTH

    def test_upgma_is_correct_when_the_clock_holds(self):
        """It is not a broken method, it is a method with an assumption."""
        names = ["A", "B", "C", "D"]
        # Ultrametric: every leaf equidistant from the root.
        matrix = [[0, 2, 4, 4], [2, 0, 4, 4], [4, 4, 0, 2], [4, 4, 2, 0]]
        truth = {frozenset(["A", "B"]), frozenset(["C", "D"])}
        assert truth <= clades(upgma(names, matrix))

    def test_two_sequences_make_a_tree(self):
        assert len(neighbour_joining(["A", "B"], [[0, 1], [1, 0]]).leaves()) == 2

    def test_one_sequence_is_refused(self):
        with pytest.raises(PhyloError):
            neighbour_joining(["A"], [[0]])

    def test_newick_output_is_produced(self):
        newick = neighbour_joining(self.NAMES, self.MATRIX).newick()
        assert newick.startswith("(") and "A" in newick

    def test_saturated_distances_are_refused_with_the_pair_named(self):
        """`jukes_cantor` returns inf at p >= 0.75 on purpose, and the tree builders
        used to turn that into `TypeError: cannot unpack non-iterable NoneType`."""
        names = ["ingroup", "outgroup"]
        with pytest.raises(PhyloError, match="saturated"):
            neighbour_joining(
                [*names, "third"], [[0, math.inf, 0.1], [math.inf, 0, 0.2], [0.1, 0.2, 0]]
            )
        with pytest.raises(PhyloError, match="outgroup"):
            upgma(names, [[0, math.inf], [math.inf, 0]])

    def test_the_documented_saturation_case_end_to_end(self):
        """The README advertises `inf` at p >= 0.75, so four random sequences - which
        saturate against each other - must give an explanation, not a crash."""
        rng = random.Random(99)
        sequences = ["".join(rng.choice("ACGT") for _ in range(300)) for _ in range(4)]
        names, matrix = distance_matrix(sequences)
        assert any(math.isinf(v) for row in matrix for v in row)
        with pytest.raises(PhyloError, match="too divergent"):
            neighbour_joining(names, matrix)

    def test_uncorrected_distances_still_build_a_tree(self):
        """The error tells the user to pass correct=False, so that has to work."""
        rng = random.Random(99)
        sequences = ["".join(rng.choice("ACGT") for _ in range(300)) for _ in range(4)]
        names, matrix = distance_matrix(sequences, correct=False)
        assert len(neighbour_joining(names, matrix).leaves()) == 4

    @pytest.mark.parametrize(
        ("label", "expected"),
        [
            ("MW183409.1", "MW183409.1"),
            ("MW183409.1 Punjab (2019)", "'MW183409.1 Punjab (2019)'"),
            ("ON312781.1:Sindh", "'ON312781.1:Sindh'"),
            ("MW183409.1 | Punjab | 2019", "'MW183409.1 | Punjab | 2019'"),
            ("O'Brien isolate", "'O''Brien isolate'"),
        ],
    )
    def test_newick_labels_are_quoted_when_they_have_to_be(self, label, expected):
        assert quote_newick(label) == expected

    def test_real_export_labels_produce_parseable_newick(self):
        """`export_alignment.py` writes exactly these labels, and unquoted they turn
        one taxon into several that any tree viewer misreads."""
        names = ["MW183409.1 | Punjab | 2019", "ON312781.1 | Sindh | 2021"]
        newick = neighbour_joining(names, [[0, 0.04], [0.04, 0]]).newick()
        assert newick.count("'") == 4
        # Outside the quoted labels, the only separator left is the one comma Newick
        # uses between the two children.
        assert (
            newick.replace("'MW183409.1 | Punjab | 2019'", "")
            .replace("'ON312781.1 | Sindh | 2021'", "")
            .count(",")
            == 1
        )


class TestAtlas:
    @staticmethod
    def isolates(seed: int = 5):
        rng = random.Random(seed)
        out = []
        for period, rate in (("2024-Q1", 0.02), ("2025-Q1", 0.10), ("2026-Q1", 0.35)):
            for i in range(100):
                s = list(REF)
                if rng.random() < rate:
                    s[20] = "G"
                out.append(
                    Isolate(
                        f"{period}-{i}",
                        "".join(s),
                        period=period,
                        location=rng.choice(["Multan", "Vehari", "Bahawalpur"]),
                    )
                )
        return out

    def test_a_variant_is_found_and_named_conventionally(self):
        atlas = build_atlas(self.isolates(), REF)
        assert atlas[0].label == "C21G"

    def test_unaligned_isolates_are_refused(self):
        with pytest.raises(ValueError):
            build_atlas([Isolate("x", "ACGT")], REF)

    def test_singletons_are_dropped(self):
        """A variant seen once in a few hundred genomes is more likely a sequencing
        error than a lineage."""
        isolates = [Isolate(f"i{i}", REF, period="2024") for i in range(200)]
        odd = list(REF)
        odd[5] = "T" if REF[5] != "T" else "A"
        isolates.append(Isolate("odd", "".join(odd), period="2024"))
        assert build_atlas(isolates, REF, min_frequency=0.01) == []

    def test_a_rising_variant_is_reported_as_emerging(self):
        emerging = emerging_variants(build_atlas(self.isolates(), REF))
        assert emerging
        assert emerging[0].last_frequency > emerging[0].first_frequency
        assert emerging[0].change >= 0.10

    def test_a_stable_variant_is_not_emerging(self):
        """Background variation at a constant 40% is not news."""
        rng = random.Random(1)
        isolates = []
        for period in ("2024-Q1", "2025-Q1", "2026-Q1"):
            for i in range(100):
                s = list(REF)
                if rng.random() < 0.4:
                    s[20] = "G"
                isolates.append(Isolate(f"{period}-{i}", "".join(s), period=period))
        assert emerging_variants(build_atlas(isolates, REF)) == []

    def test_sampling_noise_alone_does_not_trigger_an_alert(self):
        """A real false positive this caught: a variant sitting at a constant 40% was
        reported as rising, because two seasons of 100 genomes happened to land at 33%
        and 45%. Effect size alone is not evidence."""
        variants = build_atlas(
            [
                Isolate(f"a{i}", REF if i >= 33 else REF[:20] + "G" + REF[21:], period="2024")
                for i in range(100)
            ]
            + [
                Isolate(f"b{i}", REF if i >= 45 else REF[:20] + "G" + REF[21:], period="2025")
                for i in range(100)
            ],
            REF,
        )
        assert emerging_variants(variants) == []  # z below 1.96
        assert emerging_variants(variants, min_z=0.0)  # effect size alone would fire

    def test_a_large_rise_clears_significance_easily(self):
        emerging = emerging_variants(build_atlas(self.isolates(), REF))
        assert emerging[0].z > 1.96

    @pytest.mark.parametrize(
        ("a", "na", "b", "nb", "significant"),
        [
            (33, 100, 45, 100, False),  # noise
            (2, 100, 42, 100, True),  # a lineage winning
            (5, 10, 8, 10, False),
        ],  # a big jump on ten samples is not evidence
    )
    def test_two_proportion_z(self, a, na, b, nb, significant):
        assert (two_proportion_z(a, na, b, nb) >= 1.96) is significant

    def test_small_periods_are_excluded_from_the_trend(self):
        """A variant at '100%' in a period with three sequences is not at 100% - it
        is unmeasured."""
        isolates = [Isolate(f"a{i}", REF, period="2024") for i in range(50)]
        mutant = list(REF)
        mutant[20] = "G"
        isolates += [Isolate(f"b{i}", "".join(mutant), period="2025") for i in range(3)]
        assert emerging_variants(build_atlas(isolates, REF), min_samples=10) == []

    def test_fold_change_from_zero_is_undefined_not_infinite(self):
        """'Newly detected' and 'exploding' are different claims."""
        isolates = [Isolate(f"a{i}", REF, period="2024") for i in range(50)]
        mutant = list(REF)
        mutant[20] = "G"
        isolates += [Isolate(f"b{i}", REF, period="2025") for i in range(25)]
        isolates += [Isolate(f"c{i}", "".join(mutant), period="2025") for i in range(25)]
        emerging = emerging_variants(build_atlas(isolates, REF))
        assert emerging and emerging[0].fold is None

    def test_geographic_spread_distinguishes_local_from_widespread(self):
        atlas = build_atlas(self.isolates(), REF)
        spread = geographic_spread(atlas[0])
        assert spread["locations"] == 3
        assert spread["spreading"] is True

    def test_a_single_location_is_concentrated(self):
        isolates = [
            Isolate(
                f"i{i}",
                REF if i % 2 else REF[:20] + "G" + REF[21:],
                period="2026",
                location="Multan",
            )
            for i in range(100)
        ]
        spread = geographic_spread(build_atlas(isolates, REF)[0])
        assert spread["concentration"] == 1.0
        assert spread["spreading"] is False


class TestClassification:
    @staticmethod
    def fitted():
        return StrainClassifier(k=4).fit([("Burewala", REF)] * 5 + [("Rajasthan", REF[::-1])] * 5)

    def test_a_known_strain_is_assigned(self):
        assert self.fitted().classify(REF).strain == "Burewala"

    def test_a_novel_genome_is_refused_rather_than_forced(self):
        """The novel strain is the one the system exists to find, and forced-choice
        classification is exactly what hides it."""
        rng = random.Random(9)
        novel = "".join(rng.choice("ACGT") for _ in range(len(REF)))
        result = self.fitted().classify(novel)
        assert result.strain is None
        assert result.novel
        assert "exceeds the threshold" in result.reason

    def test_an_unfitted_classifier_refuses_to_guess(self):
        with pytest.raises(ValueError):
            StrainClassifier().classify(REF)

    def test_kmer_profiles_skip_ambiguous_windows(self):
        assert kmer_profile("ACGTNACGT", k=4)["ACGT"] == 2

    def test_cosine_distance_is_length_invariant(self):
        """A partial genome should be classified by composition, not pushed away from
        every centroid for being short."""
        full = kmer_profile(REF * 3, k=4)
        partial = kmer_profile(REF, k=4)
        assert cosine_distance(full, partial) < 0.01

    def test_an_empty_profile_is_maximally_distant(self):
        assert cosine_distance(kmer_profile("", k=4), kmer_profile(REF, k=4)) == 1.0


class TestRecombination:
    def test_a_block_swap_is_detected(self):
        """Recombination moves a block, so similarity to one parent rises while
        similarity to the other falls at a specific coordinate."""
        rng = random.Random(4)
        a = "".join(rng.choice("ACGT") for _ in range(1000))
        b = "".join(rng.choice("ACGT") for _ in range(1000))
        recombinant = a[:500] + b[500:]
        result = recombination_signal(recombinant, a, b, window=100, step=50)
        assert result["recombinant"]
        assert result["switches"] >= 1
        assert any(400 <= bp <= 600 for bp in result["breakpoints"])

    def test_a_non_recombinant_shows_no_switch(self):
        rng = random.Random(4)
        a = "".join(rng.choice("ACGT") for _ in range(1000))
        b = "".join(rng.choice("ACGT") for _ in range(1000))
        assert recombination_signal(a, a, b, window=100, step=50)["switches"] == 0

    def test_unaligned_sequences_are_refused(self):
        with pytest.raises(ValueError):
            recombination_signal("ACGT", "ACGT", "ACG")
