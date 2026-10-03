"""Export the aligned, clonally-deduplicated CLCuV set as FASTA.

Written so another project can use real sequences without importing this one. The output is
committed there as input data with its provenance in the header, rather than the two repos
sharing a runtime dependency.

**One sequence per haplotype, not one per record.** Conservation measured across clonal
duplicates is inflated by exactly the factor this repository exists to point out: a 2021
Punjab submission of eight genomes that collapses to one haplotype would otherwise vote
eight times on how conserved a site is. A primer designed against that number has been told
a position is safer than the evidence supports.

    uv run python scripts/export_alignment.py ../primer-designer/data/clcuv_aligned.fasta
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from real_data import SPECIES, aligned_corpus, fetch, parse_genbank, year_of  # noqa: E402

from clcuv.atlas import Isolate  # noqa: E402
from clcuv.haplotype import collapse_clonal  # noqa: E402


def main() -> int:
    dest = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("clcuv_aligned.fasta")

    records = parse_genbank(fetch().read_text(encoding="utf-8", errors="replace"))
    muv = [r for r in records if r["organism"].startswith(SPECIES)]
    print(f"{len(records)} records, {len(muv)} are {SPECIES}")

    started = time.time()
    # Same cached alignment the analysis uses, progress reported on stderr: the export
    # recomputing it was six minutes spent on a byte-identical result.
    aligned = aligned_corpus([r["sequence"] for r in muv])
    print(f"  {time.time() - started:.1f}s, {len(aligned[0])} columns")

    isolates = [
        Isolate(
            name=r["accession"],
            sequence=sequence,
            period=year_of(r["date"]),
            location=r["country"],
            host=r["host"],
        )
        for r, sequence in zip(muv, aligned, strict=True)
    ]

    collapsed, report = collapse_clonal(isolates)
    print(f"{len(collapsed)} of {len(isolates)} survive collapsing clonal groups")
    print(f"  {report.summary()}")

    width = len(collapsed[0].sequence)
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", encoding="utf-8", newline="\n") as fh:
        fh.write(
            "; Cotton leaf curl Multan virus, NCBI GenBank, one sequence per haplotype\n"
            f"; {len(collapsed)} of {len(isolates)} records, clonal groups collapsed\n"
            f"; aligned by clcuv-surveillance/src/clcuv/align.py to {width} columns\n"
        )
        for isolate in collapsed:
            loc = isolate.location or "?"
            fh.write(f">{isolate.name} | {loc} | {isolate.period or '?'}\n")
            seq = isolate.sequence
            for i in range(0, len(seq), 70):
                fh.write(seq[i : i + 70] + "\n")

    print(f"wrote {dest} ({dest.stat().st_size:,} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
