"""Run the whole pipeline on real GenBank genomes, end to end.

    uv run python scripts/real_data.py fetch                # genomes from NCBI, cached
    uv run python scripts/real_data.py analyse              # align, atlas, emergence test
    uv run python scripts/real_data.py analyse --no-cache   # realign from scratch

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

import gzip
import hashlib
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from clcuv.align import align, alignment_report, stderr_progress  # noqa: E402
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


# --- alignment: progress, and not doing it twice --------------------------

ALIGNMENT_CACHE = DATA / "clcuv_aligned.fasta.gz"
BAND = 120


def _digest(sequences: list[str]) -> str:
    """Identity of this exact input to the aligner, order included."""
    payload = hashlib.sha256()
    payload.update(f"v1 band={BAND} n={len(sequences)}\n".encode())
    for sequence in sequences:
        payload.update(sequence.encode())
        payload.update(b"\n")
    return payload.hexdigest()


def _read_cached_alignment(sequences: list[str]) -> list[str] | None:
    """The cached alignment, if it is an alignment of exactly these sequences.

    The cache is keyed by a digest of the input, and then *verified*: every row with
    its gaps removed must be the sequence it claims to align, and all rows must be the
    same width. A stale or hand-edited cache is therefore not a wrong result, it is a
    cache miss.
    """
    if not ALIGNMENT_CACHE.exists():
        return None
    try:
        with gzip.open(ALIGNMENT_CACHE, "rt", encoding="utf-8") as handle:
            header = handle.readline().strip()
            rows = [line.strip() for line in handle if line.strip()]
    except (OSError, EOFError, UnicodeDecodeError):
        return None

    if not header.startswith(">") or header[1:] != _digest(sequences):
        return None
    if len(rows) != len(sequences) or len({len(row) for row in rows}) != 1:
        return None
    if any(row.replace("-", "") != sequence for row, sequence in zip(rows, sequences, strict=True)):
        return None
    return rows


def _write_cached_alignment(sequences: list[str], aligned: list[str]) -> None:
    DATA.mkdir(exist_ok=True)
    with gzip.open(ALIGNMENT_CACHE, "wt", encoding="utf-8", newline="\n") as handle:
        handle.write(f">{_digest(sequences)}\n")
        for row in aligned:
            handle.write(f"{row}\n")


def aligned_corpus(sequences: list[str], *, use_cache: bool = True) -> list[str]:
    """Align the corpus, reusing a cached alignment of the same input if there is one.

    The committed `data/clcuv.gb` never changes between runs, so neither does its
    alignment. Recomputing it on every run spends minutes to arrive at a byte-identical
    answer; the cache is checked against the input rather than trusted, so it cannot
    change the result. `--no-cache` forces the full computation.
    """
    if use_cache:
        cached = _read_cached_alignment(sequences)
        if cached is not None:
            print(f"\nalignment cached in {ALIGNMENT_CACHE.name}, reusing it", flush=True)
            return cached

    print(
        f"\naligning {len(sequences)} genomes (pure Python; minutes, not seconds) ...",
        flush=True,
    )
    aligned = align(sequences, band=BAND, progress=stderr_progress(len(sequences)))
    if use_cache:
        _write_cached_alignment(sequences, aligned)
    return aligned


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


def analyse(*, use_cache: bool = True) -> None:
    records = parse_genbank(fetch().read_text())
    print(f"\n{len(records)} records parsed")
    print(f"  with country : {sum(1 for r in records if r['country'])}")
    print(f"  with date    : {sum(1 for r in records if r['date'])}")

    # One species only. Mixing Multan, Kokhran and Burewala virus into one alignment
    # measures the distance between species and calls it within-species variation.
    muv = [r for r in records if r["organism"].startswith(SPECIES)]
    print(f"\n{len(muv)} are {SPECIES}; the rest are other CLCuV species and excluded")

    started = time.time()
    aligned = aligned_corpus([r["sequence"] for r in muv], use_cache=use_cache)
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
        sequences = size["sequences"]
        haplotypes = size["haplotypes"]
        print(
            f"  {period}  {location:<30} {sequences:>2} "
            f"{'sequence ' if sequences == 1 else 'sequences'} "
            f"-> {haplotypes:>2} {'haplotype ' if haplotypes == 1 else 'haplotypes'}"
            f"  (x{size['inflation']}){note}"
        )

    collapsed, report = collapse_clonal(isolates)
    print("\n ", json.dumps(report.summary()))

    atlas = build_atlas(isolates, reference)
    print(f"\n{len(atlas)} variants above 1% against the consensus")

    print("\n--- what survives each control ---")
    counts: dict[str, int] = {}
    for name, pool in (("all sequences", isolates), ("one per haplotype", collapsed)):
        # The pooled atlas is the one already built above; rebuilding it was a scan of
        # 250 x 3,162 columns for a result that cannot differ.
        variants = atlas if pool is isolates else build_atlas(pool, reference)
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
        analyse(use_cache="--no-cache" not in sys.argv)
    else:
        sys.exit(f"usage: {sys.argv[0]} [fetch|analyse] [--no-cache]")
