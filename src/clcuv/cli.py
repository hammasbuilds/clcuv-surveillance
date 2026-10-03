"""`clcuv` on the command line: a FASTA and a metadata sheet in, an answer out.

    clcuv analyse --fasta isolates.fa --metadata isolates.csv
    clcuv analyse --fasta isolates.fa --metadata isolates.csv --json
    clcuv export  --fasta isolates.fa --metadata isolates.csv --out deduped.fa

Everything here is assembly. The reason it exists is that the analysis was only
reachable by importing five modules in the right order and knowing which controls to
apply - which meant the controls were optional in practice, and the whole argument of
this package is that they are not.

`analyse` prints the emergence result at one setting and the sweep around it, together,
on purpose. A single count is a setting rather than a finding, and printing it alone
invites exactly the sentence this repository exists to argue against.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

from .align import align, stderr_progress
from .atlas import (
    build_atlas,
    emergence_sensitivity,
    emerging_variants,
    group_linked,
)
from .geo import normalise_isolates, resolved_strata
from .haplotype import collapse_clonal, effective_sample_sizes
from .io import FastaError, isolates_from, read_fasta, read_metadata, write_fasta

# A CLCuMuV isolate cannot predate the virus's description by decades. One record in the
# public corpus is dated Dec-1918; keeping it creates a surveillance period seventy
# years before the first report, holding one genome.
MIN_YEAR = 1980
MAX_YEAR = 2100


def consensus(aligned: list[str]) -> str:
    """Majority base per column, ambiguous columns marked N so they cannot be variants.

    Using one isolate as the reference is the alternative, and it is worse: every
    position where that isolate happens to be unusual becomes a variant present in
    almost every other genome.
    """
    out = []
    for column in zip(*aligned, strict=True):
        called = Counter(base for base in column if base in "ACGT")
        if not called:
            out.append("N")
            continue
        ((base, count),) = called.most_common(1)
        out.append(base if count > len(aligned) / 2 else "N")
    return "".join(out)


def _load(args) -> tuple[list, dict, list[str]]:
    """FASTA + metadata -> isolates, the location mapping, and what was dropped."""
    records = read_fasta(args.fasta)
    metadata = read_metadata(args.metadata) if args.metadata else {}
    if not metadata:
        raise SystemExit(
            "--metadata is required: an emergence test needs a period and a place for "
            "every sequence, and there is nothing in a FASTA header this can trust to "
            "supply them."
        )

    isolates, skipped = isolates_from(records, metadata)

    implausible = [i.name for i in isolates if not (MIN_YEAR <= int(i.period[:4]) <= MAX_YEAR)]
    isolates = [i for i in isolates if i.name not in set(implausible)]

    notes = []
    if skipped:
        notes.append(f"{len(skipped)} sequences have no usable period or location: {skipped[:5]}")
    if implausible:
        notes.append(f"{len(implausible)} dropped for an implausible year: {implausible[:5]}")
    return isolates, records, notes


def _aligned_sequences(isolates, *, band: int) -> list[str]:
    widths = {len(i.sequence) for i in isolates}
    if len(widths) == 1 and any("-" in i.sequence for i in isolates):
        return [i.sequence for i in isolates]  # already aligned
    # Unaligned input means a full centre-star alignment, which is minutes for a few
    # hundred genomes. Progress goes to stderr so stdout (including --json) stays clean.
    return align(
        [i.sequence for i in isolates],
        band=band,
        names=[i.name for i in isolates],
        progress=stderr_progress(len(isolates)),
    )


def analyse(args) -> int:
    isolates, records, notes = _load(args)
    if len(isolates) < 2:
        raise SystemExit(
            f"only {len(isolates)} of {len(records)} sequences have both a period and a "
            "location, which is not enough to compare anything."
        )

    aligned = _aligned_sequences(isolates, band=args.band)
    isolates = [
        type(i)(name=i.name, sequence=s, period=i.period, location=i.location, host=i.host)
        for i, s in zip(isolates, aligned, strict=True)
    ]
    reference = consensus(aligned)

    normalised, mapping = normalise_isolates(isolates)
    eligible = resolved_strata(mapping)
    collapsed, collapse_report = collapse_clonal(normalised)

    atlas = build_atlas(normalised, reference)
    collapsed_atlas = build_atlas(collapsed, reference)

    def counts(variants):
        pooled = emerging_variants(variants, min_samples=args.min_samples, compare=args.compare)
        stratified = emerging_variants(
            variants,
            min_samples=args.min_samples,
            compare=args.compare,
            stratify=True,
            eligible_locations=eligible,
        )
        return pooled, stratified

    pooled, stratified = counts(atlas)
    pooled_c, stratified_c = counts(collapsed_atlas)

    sweep = {
        "all_sequences": emergence_sensitivity(
            atlas, compare=args.compare, eligible_locations=eligible
        ),
        "one_per_haplotype": emergence_sensitivity(
            collapsed_atlas, compare=args.compare, eligible_locations=eligible
        ),
    }

    result = {
        "sequences_read": len(records),
        "surveillable": len(normalised),
        "notes": notes,
        "alignment": {"width": len(reference), "uncalled_columns": reference.count("N")},
        "locations": {
            "raw_strings": len(mapping),
            "strata": len({p.stratum for p in mapping.values()}),
            "resolved_strata": sorted(eligible),
            "unresolved_strata": sorted({p.stratum for p in mapping.values() if not p.resolved}),
            "flagged": sorted(
                {p.raw for p in mapping.values() if p.flags and p.flags[0] != "country-only"}
            ),
        },
        "clonality": collapse_report.summary(),
        "settings": {"min_samples": args.min_samples, "compare": args.compare},
        "emerging": {
            "all_sequences": {"pooled": len(pooled), "stratified": len(stratified)},
            "one_per_haplotype": {"pooled": len(pooled_c), "stratified": len(stratified_c)},
        },
        "distinct_events": {
            "all_sequences": len(group_linked(stratified)),
            "one_per_haplotype": len(group_linked(stratified_c)),
        },
        "sensitivity": sweep,
        "variants": [e.summary() for e in stratified_c[: args.top]],
    }

    if args.json:
        json.dump(result, sys.stdout, indent=2)
        sys.stdout.write("\n")
        return 0

    _print_report(result, normalised, mapping)
    return 0


def _print_report(result: dict, isolates, mapping) -> None:
    print(f"{result['sequences_read']} sequences read, {result['surveillable']} surveillable")
    for note in result["notes"]:
        print(f"  note: {note}")

    alignment = result["alignment"]
    print(f"\nalignment: {alignment['width']} columns, {alignment['uncalled_columns']} uncalled")

    locations = result["locations"]
    print(
        f"\nlocations: {locations['raw_strings']} raw strings -> {locations['strata']} strata "
        f"({len(locations['resolved_strata'])} name a province and may confirm a rise)"
    )
    if locations["unresolved_strata"]:
        print(
            "  unresolved (counted in the pooled test, cannot confirm): "
            + ", ".join(locations["unresolved_strata"])
        )
    for flagged in locations["flagged"]:
        print(f"  flagged: {flagged!r} names a division that is not in that country")

    print("\nindependent genomes per stratum:")
    sizes = effective_sample_sizes(isolates)
    for (period, location), size in sorted(sizes.items()):
        if size["sequences"] < 5:
            note = "   <- too few to test"
        elif size["haplotypes"] < 5:
            note = "   <- too clonal to test"
        else:
            note = ""
        sequences, haplotypes = size["sequences"], size["haplotypes"]
        print(
            f"  {period}  {location:<26} {sequences:>3} "
            f"{'sequence ' if sequences == 1 else 'sequences'} -> "
            f"{haplotypes:>3} {'haplotype ' if haplotypes == 1 else 'haplotypes'}"
            f"  (x{size['inflation']}){note}"
        )

    print(f"\n  {json.dumps(result['clonality'])}")

    settings = result["settings"]
    print(
        f"\n--- what survives each control (min_samples={settings['min_samples']}, "
        f"compare={settings['compare']}) ---"
    )
    labels = (("all_sequences", "all sequences"), ("one_per_haplotype", "one per haplotype"))
    for key, label in labels:
        counts = result["emerging"][key]
        print(
            f"  {label:<20} pooled={counts['pooled']:<4} stratified={counts['stratified']:<4} "
            f"distinct events={result['distinct_events'][key]}"
        )

    print("\n--- and how much of that is the threshold talking ---")
    print("  min_samples  all sequences (pooled/strat)   one per haplotype (pooled/strat)")
    sensitivity = result["sensitivity"]
    rows = zip(sensitivity["all_sequences"], sensitivity["one_per_haplotype"], strict=True)
    for all_row, hap_row in rows:
        print(
            f"  {all_row['min_samples']:>11}  {all_row['pooled']:>13}/{all_row['stratified']:<15}"
            f"{hap_row['pooled']:>13}/{hap_row['stratified']}"
        )
    print(
        "\n  A count that moves across this table is a setting, not a finding. Read the\n"
        "  rows where the endpoint periods are large enough to mean anything."
    )


def export(args) -> int:
    """Write one sequence per haplotype, aligned, for a downstream tool."""
    isolates, records, notes = _load(args)
    for note in notes:
        print(f"note: {note}")

    aligned = _aligned_sequences(isolates, band=args.band)
    isolates = [
        type(i)(name=i.name, sequence=s, period=i.period, location=i.location, host=i.host)
        for i, s in zip(isolates, aligned, strict=True)
    ]
    normalised, _ = normalise_isolates(isolates)
    collapsed, report = collapse_clonal(normalised)

    written = write_fasta(
        args.out,
        ((f"{i.name} | {i.location or '?'} | {i.period or '?'}", i.sequence) for i in collapsed),
        header=[
            "one sequence per haplotype per (period, location)",
            f"{len(collapsed)} of {len(normalised)} records, "
            f"aligned to {len(collapsed[0].sequence)} columns by clcuv-surveillance",
        ],
    )
    print(f"{report.summary()}")
    print(f"wrote {written} sequences to {args.out}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="clcuv",
        description="Genomic surveillance for Cotton Leaf Curl Virus.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    def common(sub):
        sub.add_argument("--fasta", required=True, type=Path, help="sequences, aligned or not")
        sub.add_argument(
            "--metadata",
            type=Path,
            help="CSV with a name column plus period/date and location/country",
        )
        sub.add_argument(
            "--band",
            type=int,
            default=120,
            help=(
                "aligner band width in bases (default 120). Runtime is roughly linear "
                "in it: on 70 unaligned 2.7 kb genomes --band 40 took 57s against "
                "1m58s at 120. It is not free - a band narrower than the largest indel "
                "between a sequence and the centre cannot place that indel, and those "
                "same 70 genomes aligned to 2,968 columns instead of 3,162. Ignored "
                "when the input is already aligned."
            ),
        )

    run = subparsers.add_parser("analyse", help="which variant is rising, and is that real")
    common(run)
    run.add_argument("--min-samples", type=int, default=8, dest="min_samples")
    run.add_argument("--compare", choices=("extremes", "adjacent"), default="extremes")
    run.add_argument("--top", type=int, default=10, help="variants to list in --json")
    run.add_argument("--json", action="store_true", help="machine-readable output")
    run.set_defaults(func=analyse)

    out = subparsers.add_parser("export", help="write the deduplicated alignment as FASTA")
    common(out)
    out.add_argument("--out", required=True, type=Path)
    out.set_defaults(func=export)

    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except FileNotFoundError as error:
        print(f"error: {error.filename or error}: no such file", file=sys.stderr)
        return 1
    except (FastaError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
