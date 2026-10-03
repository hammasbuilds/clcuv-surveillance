"""Reading sequences and metadata off disk.

The package could do the analysis and could not open a file, which meant every user had
to write this themselves before they could use any of it - and write it correctly,
because a FASTA reader that silently drops a record or keeps a `*` in a sequence moves
every downstream coordinate.

Deliberately small: FASTA in, a metadata table joined onto it, nothing else. Biopython
does all of this better and this package has no dependencies, which is the whole trade.
"""

from __future__ import annotations

import csv
import re
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

from .atlas import Isolate, check_period


class FastaError(ValueError):
    """The file is not FASTA, or not FASTA this can use."""


@dataclass(frozen=True)
class Record:
    name: str
    sequence: str
    description: str = ""


# Characters a sequence may contain. Gaps are kept, because an already-aligned FASTA is
# the common input; everything else is rejected rather than stripped, since silently
# deleting a character shifts every position after it.
ALLOWED = set("ACGTURYKMSWBDHVNacgturykmswbdhvn-.")


def parse_fasta(text: str, *, allow_gaps: bool = True) -> list[Record]:
    """Parse FASTA text into records.

    `;` comment lines are skipped: `scripts/export_alignment.py` writes its provenance
    header that way, so this package's own output has to round-trip through here.
    """
    records: list[Record] = []
    name = ""
    description = ""
    chunks: list[str] = []
    seen: set[str] = set()

    def flush() -> None:
        if not name:
            return
        sequence = "".join(chunks)
        if not sequence:
            raise FastaError(f"{name!r} has a header but no sequence")
        bad = sorted(set(sequence) - ALLOWED)
        if bad:
            raise FastaError(
                f"{name!r} contains {bad[:5]}, which is not nucleotide sequence. "
                "If this is protein, this package works on nucleotides only."
            )
        if not allow_gaps and ("-" in sequence or "." in sequence):
            raise FastaError(f"{name!r} contains gaps, but unaligned sequence was expected")
        if name in seen:
            raise FastaError(
                f"{name!r} appears twice. Duplicate names silently overwrite each other "
                "in every join downstream, so they are refused here."
            )
        seen.add(name)
        records.append(Record(name=name, sequence=sequence.upper(), description=description))

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith(">"):
            flush()
            header = line[1:].strip()
            name, _, description = header.partition(" ")
            description = description.strip()
            chunks = []
            if not name:
                raise FastaError("a record has an empty name")
        elif name:
            chunks.append(line.replace(" ", ""))
        else:
            raise FastaError("sequence data appears before the first '>' header")

    flush()
    if not records:
        raise FastaError("no FASTA records found - is this a FASTA file?")
    return records


def read_fasta(path: str | Path, *, allow_gaps: bool = True) -> list[Record]:
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    return parse_fasta(text, allow_gaps=allow_gaps)


def write_fasta(
    path: str | Path,
    records: Iterable[tuple[str, str]],
    *,
    header: Sequence[str] = (),
    width: int = 70,
) -> int:
    """Write FASTA with LF endings, whatever the platform.

    CRLF here is not cosmetic: a file written on Windows and committed shows as changed
    on every refetch, and some downstream parsers keep the `\\r` as sequence.
    """
    path = Path(path)
    if path.parent != Path(""):
        path.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for line in header:
            handle.write(f"; {line}\n")
        for name, sequence in records:
            handle.write(f">{name}\n")
            for i in range(0, len(sequence), width):
                handle.write(sequence[i : i + width] + "\n")
            written += 1
    return written


# Column names accepted for each field, in order of preference. GenBank exports, lab
# spreadsheets and this package's own output all name these differently, and rejecting a
# sheet because it says "district" rather than "location" helps nobody.
COLUMNS: dict[str, tuple[str, ...]] = {
    "name": ("name", "accession", "id", "isolate", "sequence_id", "seqid"),
    "period": ("period", "year", "date", "collection_date", "season"),
    "location": ("location", "country", "region", "district", "province", "geo_loc_name"),
    "host": ("host", "cultivar", "variety"),
}


def _column(fieldnames: Sequence[str], field: str) -> str | None:
    lowered = {n.strip().lower(): n for n in fieldnames if n}
    for candidate in COLUMNS[field]:
        if candidate in lowered:
            return lowered[candidate]
    return None


def read_metadata(path: str | Path) -> dict[str, dict[str, str]]:
    """Read a CSV of per-sequence metadata, keyed by sequence name.

    Only four fields are used, and the name column is the only one that must be there:
    without it there is nothing to join on, and a join on row order is the kind of thing
    that works on the sample and silently mismatches on the real file.
    """
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        name_column = _column(fieldnames, "name")
        if not name_column:
            raise ValueError(
                f"no name column in {path}: expected one of {COLUMNS['name']}, found {fieldnames}"
            )
        columns = {field: _column(fieldnames, field) for field in COLUMNS}

        out: dict[str, dict[str, str]] = {}
        for row in reader:
            key = (row.get(name_column) or "").strip()
            if not key:
                continue
            if key in out:
                raise ValueError(f"{key!r} appears twice in {path}")
            out[key] = {
                field: (row.get(column) or "").strip() if column else ""
                for field, column in columns.items()
            }
    return out


def period_from(value: str) -> str:
    """Coerce a date column into a sortable period, or "" if it holds no year.

    GenBank `collection_date` has no single format - `2019`, `May-2019`, `01-May-2019`
    and `2015-01` all appear in one corpus - and taking the last four characters yields
    `5-01` for the fourth, which becomes its own surveillance period and splits the data
    in half without saying so. A value that is already sortable is left alone.
    """
    text = (value or "").strip()
    if not text:
        return ""
    try:
        check_period(text)
    except ValueError:
        match = re.search(r"(19|20)\d{2}", text)
        return match.group(0) if match else ""
    return text


def isolates_from(
    records: Sequence[Record],
    metadata: dict[str, dict[str, str]] | None = None,
    *,
    period_of=period_from,
) -> tuple[list[Isolate], list[str]]:
    """Join records to metadata and build isolates.

    Returns the isolates and the names that had no usable metadata, because a surveillance
    run that quietly analyses 40 of 60 sequences has answered a different question from
    the one it was asked.
    """
    metadata = metadata or {}
    isolates: list[Isolate] = []
    skipped: list[str] = []

    for record in records:
        row = metadata.get(record.name)
        if row is None:
            # `>MW183409.1 | Punjab | 2019` - the name is the first token, so a sheet
            # keyed on the accession still joins.
            row = metadata.get(record.name.split()[0]) if record.name else None
        if not row:
            skipped.append(record.name)
            continue

        period = row.get("period", "")
        if period_of is not None:
            period = period_of(period)
        if not period or not row.get("location"):
            skipped.append(record.name)
            continue

        try:
            check_period(period)
        except ValueError as error:
            raise ValueError(f"{record.name}: {error}") from error

        isolates.append(
            Isolate(
                name=record.name,
                sequence=record.sequence,
                period=period,
                location=row.get("location", ""),
                host=row.get("host", ""),
            )
        )

    return isolates, skipped
