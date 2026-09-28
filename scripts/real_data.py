"""Run the whole pipeline on real GenBank genomes, end to end.

    uv run python scripts/real_data.py fetch     # ~60 genomes from NCBI, cached
    uv run python scripts/real_data.py analyse   # align, atlas, emergence test

Nothing here is in the library. The library takes aligned sequences and metadata; this
is the glue that gets them out of NCBI, and it is a script rather than a module because
it depends on a network service whose behaviour is not ours to test.

Only the standard library is used, so it runs in a bare checkout.

Run `analyse` to see the finding this project's last two controls exist for: a double
digit number of variants look like they are emerging, most of them survive stratification
by region once mis-spelled and mis-cased location strings are merged into one, and
**none of them survive being asked how many independent genomes are behind them**. The
exact counts are printed, not written down here, because they move as the corpus grows -
an earlier version of this paragraph said "nine" after the corpus had already made it "52".
"""

from __future__ import annotations

import json
import re
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from clcuv.align import align, alignment_report  # noqa: E402
from clcuv.atlas import Isolate, build_atlas, emerging_variants  # noqa: E402
from clcuv.cli import MAX_YEAR, MIN_YEAR  # noqa: E402
from clcuv.geo import normalise_isolates, resolved_strata  # noqa: E402
from clcuv.haplotype import collapse_clonal, effective_sample_sizes  # noqa: E402

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
QUERY = '"Cotton leaf curl Multan virus"[Organism] AND 2500:3000[SLEN]'
# Every record the query matches. The first version asked for 60, which was a round
# number and not a reason - and the headline finding is that there are not enough
# INDEPENDENT genomes to support an emergence claim, which is a statement about sample
# size. Making it 60 by choice and then reporting a sample-size limit would be circular.
RETMAX = 300
# Committed to the repo: 532 KB, and it makes the headline result reproducible with
# no network at all. Delete it and `fetch` downloads it again.
DATA = Path(__file__).resolve().parent.parent / "data"

SPECIES = "Cotton leaf curl Multan virus"


# --- fetch ----------------------------------------------------------------


def _get(endpoint: str, **params) -> str:
    params.setdefault("db", "nucleotide")
    params.setdefault("tool", "clcuv-surveillance")
    url = f"{EUTILS}/{endpoint}.fcgi?{urllib.parse.urlencode(params)}"
    with urllib.request.urlopen(url, timeout=120) as response:  # noqa: S310
        return response.read().decode("utf-8", "replace")


def fetch(force: bool = False) -> Path:
    """Download GenBank records, once. NCBI asks for three requests a second at most."""
    DATA.mkdir(exist_ok=True)
    target = DATA / "clcuv.gb"
    if target.exists() and not force:
        print(f"using {target} ({target.stat().st_size:,} bytes)")
        print("pass --force to re-download, e.g. after raising RETMAX")
        return target

    search = _get("esearch", term=QUERY, retmax=RETMAX)
    ids = re.findall(r"<Id>(\d+)</Id>", search)
    total = re.search(r"<Count>(\d+)</Count>", search)
    print(f"{total.group(1) if total else '?'} records match; fetching {len(ids)}")

    time.sleep(0.4)
    target.write_text(_get("efetch", id=",".join(ids), rettype="gb", retmode="text"))
    print(f"wrote {target} ({target.stat().st_size:,} bytes)")
    return target


# --- parsing --------------------------------------------------------------


def parse_genbank(text: str) -> list[dict]:
    """Accession, organism, sequence, and the two qualifiers that make it surveillance.

    A record missing `/country` or `/collection_date` is kept and reported rather than
    dropped, because how much of GenBank lacks usable metadata is itself a finding.
    """
    records = []
    for block in text.split("\n//\n"):
        if "ORIGIN" not in block:
            continue
        accession = re.search(r"^VERSION\s+(\S+)", block, re.M)
        organism = re.search(r"^\s+ORGANISM\s+(.+)$", block, re.M)
        country = re.search(r'/(?:country|geo_loc_name)="([^"]+)"', block)
        date = re.search(r'/collection_date="([^"]+)"', block)
        host = re.search(r'/host="([^"]+)"', block)

        sequence = "".join(
            re.sub(r"[^acgtnACGTN]", "", line)
            for line in block.split("ORIGIN", 1)[1].splitlines()[1:]
        ).upper()
        if not accession or not sequence:
            continue

        records.append(
            {
                "accession": accession.group(1),
                "organism": organism.group(1).strip() if organism else "",
                "country": country.group(1) if country else "",
                "date": date.group(1) if date else "",
                "host": host.group(1) if host else "",
                "sequence": sequence,
            }
        )
    return records


def year_of(date: str) -> str:
    """The year out of a GenBank collection_date, which has no single format.

    Seen in this dataset alone: `2019`, `May-2019`, `01-May-2019`, `2015-01`. Taking the
    last four characters produces `5-01` for the fourth, which then becomes its own
    surveillance period and quietly splits the data.
    """
    match = re.search(r"(19|20)\d{2}", date)
    return match.group(0) if match else ""


# --- analysis -------------------------------------------------------------


def _consensus(aligned: list[str]) -> str:
    """Majority base per column, ambiguous columns marked N so they cannot be variants.

    Using one isolate as the reference is the alternative, and it is worse: every
    position where that isolate happens to be unusual becomes a variant present in
    almost every other genome. An early run of this made exactly that mistake and
    reported 620 variants at 98% frequency.
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


def analyse() -> None:
    records = parse_genbank(fetch().read_text())
    print(f"\n{len(records)} records parsed")
    print(f"  with country : {sum(1 for r in records if r['country'])}")
    print(f"  with date    : {sum(1 for r in records if r['date'])}")

    # One species only. Mixing Multan, Kokhran and Burewala virus into one alignment
    # measures the distance between species and calls it within-species variation.
    muv = [r for r in records if r["organism"].startswith(SPECIES)]
    print(f"\n{len(muv)} are {SPECIES}; the rest are other CLCuV species and excluded")

    print("\naligning ...", flush=True)
    started = time.time()
    aligned = align([r["sequence"] for r in muv])
    print(f"  {time.time() - started:.1f}s")
    print(" ", json.dumps(alignment_report(aligned)))

    reference = _consensus(aligned)
    print(f"  consensus: {len(reference)} bp, {reference.count('N')} uncalled")

    isolates = [
        Isolate(
            name=r["accession"],
            sequence=sequence,
            period=year_of(r["date"]),
            location=r["country"],
            host=r["host"],
        )
        for r, sequence in zip(muv, aligned, strict=True)
        if year_of(r["date"]) and r["country"] and MIN_YEAR <= int(year_of(r["date"])) <= MAX_YEAR
    ]
    print(f"\n{len(isolates)} have both a year and a place, and can be surveilled")

    # Raw `/country` strings ("Pakistan: Punjab", "Pakistan: Punjab province", "Pakistan:
    # Punjab,Bahawalpur") are spellings of the same stratum, not different ones - see
    # clcuv.geo. Stratifying on the raw string splits one province's evidence across
    # several rows and undercounts every one of them.
    isolates, mapping = normalise_isolates(isolates)
    eligible = resolved_strata(mapping)

    print("\n--- how many independent genomes are actually here? ---")
    for (period, location), size in effective_sample_sizes(isolates).items():
        note = "" if size["usable_for_statistics"] else "   <- too clonal to test"
        print(
            f"  {period}  {location:<30} {size['sequences']:>2} seqs "
            f"-> {size['haplotypes']:>2} haplotypes  (x{size['inflation']}){note}"
        )

    collapsed, report = collapse_clonal(isolates)
    print("\n ", json.dumps(report.summary()))

    atlas = build_atlas(isolates, reference)
    print(f"\n{len(atlas)} variants above 1% against the consensus")

    print("\n--- what survives each control ---")
    counts: dict[str, int] = {}
    for name, pool in (("all sequences", isolates), ("one per haplotype", collapsed)):
        variants = build_atlas(pool, reference)
        pooled = emerging_variants(variants, min_samples=8)
        stratified = emerging_variants(
            variants, min_samples=8, stratify=True, eligible_locations=eligible
        )
        counts[name] = len(stratified)
        print(f"  {name:<20} n={len(pool):<3} pooled={len(pooled):<3} stratified={len(stratified)}")

    # Computed, not written down. An earlier version said "the nine that survive" as prose.
    # The corpus then went from 60 genomes to 254, the number became 52, and the sentence
    # was still confidently saying nine - a narrative that does not read the data it
    # describes is a claim waiting to go stale.
    survived = counts["all sequences"]
    independent = counts["one per haplotype"]
    print(
        f"\nThe {survived} that survive stratification do so on pooled counts that still\n"
        "contain clonal duplicates. Collapse each clonal group to one haplotype, ask again,\n"
        f"and {independent} are left. Asked for independent evidence they have none, which is\n"
        "the honest answer this dataset supports."
    )


if __name__ == "__main__":
    command = sys.argv[1] if len(sys.argv) > 1 else "analyse"
    if command == "fetch":
        fetch(force="--force" in sys.argv)
    elif command == "analyse":
        analyse()
    else:
        sys.exit(f"usage: {sys.argv[0]} [fetch|analyse]")
