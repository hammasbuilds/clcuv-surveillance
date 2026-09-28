"""Distances and tree building.

Two distance corrections and two clustering methods, and in both cases the second
exists because the first is wrong in a way that matters for surveillance.

**p-distance vs Jukes-Cantor.** Counting differences undercounts divergence, because
the same site can mutate twice and the second change hides the first. At 5% divergence
this hardly matters; at 30% it matters a great deal, and begomovirus genomes routinely
reach that between species. Jukes-Cantor corrects for multiple hits.

**UPGMA vs neighbour-joining.** UPGMA assumes a molecular clock — every lineage
evolving at the same rate. Real viral lineages do not: a strain under host-resistance
pressure accumulates substitutions faster than one in a susceptible cultivar, which is
exactly the strain surveillance cares about. UPGMA puts fast-evolving lineages in the
wrong place, and it does so confidently. Neighbour-joining makes no clock assumption.

So UPGMA is here to be *compared against*, and the comparison is a test.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass, field


class PhyloError(ValueError):
    pass


# Characters that Newick gives a meaning to. A label containing one of these must be
# quoted or the tree reparses as a different tree - and `export_alignment.py` writes
# labels like "MW183409.1 | Punjab | 2019", which contains two of them.
NEWICK_SPECIALS = set(" \t\n()[]':;,")


def quote_newick(label: str) -> str:
    """A label safe to write into a Newick string.

    Unquoted Newick treats a space as a name terminator and `:` as the start of a
    branch length, so `MW183409.1 Punjab (2019)` silently becomes a malformed tree that
    viewers read as several taxa. Single quotes are the spec's escape, with an internal
    quote written twice.
    """
    if label and not (NEWICK_SPECIALS & set(label)):
        return label
    return "'" + label.replace("'", "''") + "'"


def _check_finite(names: Sequence[str], matrix: list[list[float]]) -> None:
    """Refuse a matrix holding `inf`, and say which pair caused it.

    `jukes_cantor` returns infinity at p >= 0.75 on purpose: the distance genuinely
    cannot be estimated there. Tree building cannot use that value - in
    neighbour-joining `inf - inf` becomes NaN, every Q comparison is then false, and
    the algorithm falls off the end with a TypeError that says nothing about the data.
    Failing here, naming the pair, is the difference between a bug and an answer.
    """
    for i, row in enumerate(matrix):
        for j, value in enumerate(row):
            if not math.isfinite(value):
                raise PhyloError(
                    f"distance between {names[i]} and {names[j]} is saturated "
                    "(p >= 0.75, Jukes-Cantor undefined), so no tree can be built from "
                    "this matrix. These sequences are too divergent to place with a "
                    "distance method - drop the outgroup, or pass correct=False to "
                    "distance_matrix to use uncorrected p-distances."
                )


def p_distance(a: str, b: str) -> float:
    """Proportion of compared sites that differ.

    Positions where either sequence is a gap or an ambiguity are skipped rather than
    counted as differences — a gap is missing data, and treating it as a substitution
    makes poorly sequenced strains look divergent.
    """
    if len(a) != len(b):
        raise PhyloError("sequences must be aligned and the same length")

    compared = differences = 0
    for x, y in zip(a.upper(), b.upper(), strict=False):
        if x not in "ACGT" or y not in "ACGT":
            continue
        compared += 1
        if x != y:
            differences += 1

    return round(differences / compared, 6) if compared else 0.0


def jukes_cantor(p: float) -> float:
    """Correct a p-distance for multiple substitutions at the same site.

        d = −(3/4) ln(1 − (4/3)p)

    Undefined at p ≥ 0.75, where the sequences are statistically indistinguishable
    from random. Returning infinity there is correct and returning a large finite
    number is not: the honest statement is that the distance cannot be estimated.
    """
    if p <= 0:
        return 0.0
    if p >= 0.75:
        return math.inf
    return round(-0.75 * math.log(1 - (4 / 3) * p), 6)


def distance_matrix(
    sequences: Sequence[str], *, names: Sequence[str] | None = None, correct: bool = True
) -> tuple[list[str], list[list[float]]]:
    names = list(names) if names else [f"seq-{i}" for i in range(len(sequences))]
    if len(names) != len(sequences):
        raise PhyloError("names and sequences must be the same length")

    n = len(sequences)
    matrix = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            p = p_distance(sequences[i], sequences[j])
            d = jukes_cantor(p) if correct else p
            matrix[i][j] = matrix[j][i] = d
    return names, matrix


@dataclass
class Node:
    name: str = ""
    children: list[tuple[Node, float]] = field(default_factory=list)

    @property
    def is_leaf(self) -> bool:
        return not self.children

    def leaves(self) -> list[str]:
        if self.is_leaf:
            return [self.name]
        out: list[str] = []
        for child, _ in self.children:
            out.extend(child.leaves())
        return out

    def newick(self) -> str:
        if self.is_leaf:
            return quote_newick(self.name)
        inner = ",".join(f"{child.newick()}:{length:.6f}" for child, length in self.children)
        return f"({inner})"


def _closest_pair(matrix: list[list[float]], active: list[int]) -> tuple[int, int]:
    best = (active[0], active[1])
    best_distance = math.inf
    for a_index, i in enumerate(active):
        for j in active[a_index + 1 :]:
            if matrix[i][j] < best_distance:
                best_distance = matrix[i][j]
                best = (i, j)
    return best


def upgma(names: Sequence[str], matrix: list[list[float]]) -> Node:
    """Average-linkage clustering. Assumes a molecular clock — see the module note."""
    if len(names) < 2:
        raise PhyloError("need at least two sequences")
    _check_finite(names, matrix)

    working = [row[:] for row in matrix]
    nodes: dict[int, Node] = {i: Node(name=n) for i, n in enumerate(names)}
    sizes: dict[int, int] = dict.fromkeys(nodes, 1)
    heights: dict[int, float] = dict.fromkeys(nodes, 0.0)
    active = list(nodes)

    while len(active) > 1:
        i, j = _closest_pair(working, active)
        height = working[i][j] / 2

        parent = Node(
            children=[
                (nodes[i], round(height - heights[i], 6)),
                (nodes[j], round(height - heights[j], 6)),
            ]
        )

        new_index = max(nodes) + 1
        nodes[new_index] = parent
        sizes[new_index] = sizes[i] + sizes[j]
        heights[new_index] = height

        for row in working:
            row.append(0.0)
        working.append([0.0] * (len(working) + 1))

        for k in active:
            if k in (i, j):
                continue
            # Average linkage, weighted by cluster size.
            merged = (working[i][k] * sizes[i] + working[j][k] * sizes[j]) / (sizes[i] + sizes[j])
            working[new_index][k] = working[k][new_index] = merged

        active = [k for k in active if k not in (i, j)] + [new_index]

    return nodes[active[0]]


def neighbour_joining(names: Sequence[str], matrix: list[list[float]]) -> Node:
    """Saitou-Nei neighbour-joining. Makes no molecular-clock assumption.

    The Q-criterion is what makes it clock-free: rather than joining the closest pair,
    it joins the pair that is closest *after* correcting for how far each is from
    everything else. A fast-evolving lineage is far from everything, and that offset
    cancels instead of dragging it to the wrong place in the tree.
    """
    n = len(names)
    if n < 2:
        raise PhyloError("need at least two sequences")
    _check_finite(names, matrix)
    if n == 2:
        root = Node()
        root.children = [
            (Node(name=names[0]), round(matrix[0][1] / 2, 6)),
            (Node(name=names[1]), round(matrix[0][1] / 2, 6)),
        ]
        return root

    working = {i: {j: matrix[i][j] for j in range(n) if j != i} for i in range(n)}
    nodes: dict[int, Node] = {i: Node(name=names[i]) for i in range(n)}
    next_index = n

    while len(working) > 2:
        active = list(working)
        size = len(active)
        totals = {i: sum(working[i][j] for j in active if j != i) for i in active}

        best_pair, best_q = None, math.inf
        for a_index, i in enumerate(active):
            for j in active[a_index + 1 :]:
                q = (size - 2) * working[i][j] - totals[i] - totals[j]
                if q < best_q:
                    best_q, best_pair = q, (i, j)

        if best_pair is None:
            # Unreachable while `_check_finite` runs first, and kept because the failure
            # it replaces was a TypeError from unpacking None several lines later.
            raise PhyloError("no pair could be joined; the distance matrix is not usable")

        i, j = best_pair
        d_ij = working[i][j]
        limb_i = 0.5 * d_ij + (totals[i] - totals[j]) / (2 * (size - 2))
        limb_j = d_ij - limb_i

        parent = Node(
            children=[
                (nodes[i], round(max(limb_i, 0.0), 6)),
                (nodes[j], round(max(limb_j, 0.0), 6)),
            ]
        )
        nodes[next_index] = parent

        working[next_index] = {}
        for k in active:
            if k in (i, j):
                continue
            d = 0.5 * (working[i][k] + working[j][k] - d_ij)
            working[next_index][k] = d
            working[k][next_index] = d
            del working[k][i]
            del working[k][j]

        del working[i]
        del working[j]
        next_index += 1

    left, right = list(working)
    root = Node(
        children=[
            (nodes[left], round(working[left][right] / 2, 6)),
            (nodes[right], round(working[left][right] / 2, 6)),
        ]
    )
    return root


def clades(tree: Node) -> set[frozenset[str]]:
    """Every set of leaves that forms a clade. The topology, stripped of branch
    lengths and of the arbitrary left/right ordering — which is what makes two trees
    comparable."""
    found: set[frozenset[str]] = set()

    def walk(node: Node) -> list[str]:
        if node.is_leaf:
            return [node.name]
        leaves: list[str] = []
        for child, _ in node.children:
            leaves.extend(walk(child))
        found.add(frozenset(leaves))
        return leaves

    walk(tree)
    return found
