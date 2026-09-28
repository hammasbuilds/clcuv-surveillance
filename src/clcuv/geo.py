"""Turning GenBank `/country` strings into strata that can be compared.

The stratified control in `atlas.py` only works if "the same location" means the same
thing in both periods. GenBank's `/country` qualifier is free text, and in this corpus
alone it arrives as:

    Pakistan
    Pakistan: Punjab
    Pakistan: Punjab province
    Pakistan: Punjab,Bahawalpur
    Pakistan: Vehari, Punjab
    Pakistan: Vehari-Punjab
    Pakistan:Faisalabad
    Pakistan: NIAB, FSD

Those are six spellings of one province plus two of its districts, and as raw strings
they are eight separate strata. Splitting one province eight ways does not make the
comparison conservative — it makes each stratum too small to test, which silently
removes the control rather than tightening it.

Two failures are distinguished here, because they need opposite treatment:

**Spelling.** `Punjab province` and `Vehari-Punjab` are Punjab. These are merged.

**Missing detail.** A bare `Pakistan` names no province at all. It cannot be merged
with `Pakistan: Punjab` (that would reintroduce the confounder the stratification
exists to remove) and it cannot be split (there is nothing to split on). So it is kept
as its own stratum and marked unresolved, and the stratified test is told not to let
an unresolved stratum *confirm* anything. It still counts in the pooled test, which is
where it belongs: it is evidence about the country and not about a province.

One inconsistency is reported rather than repaired. Eleven records are filed as
`Pakistan: Rajasthan` / `Rajesthan`, and Rajasthan is an Indian state. Reassigning them
to India would be this module inventing metadata; they are normalised for spelling,
flagged `admin1-not-in-country`, and left where the submitter put them.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass

from .atlas import Isolate

# Country name as written -> canonical country.
COUNTRY_ALIASES: dict[str, str] = {
    "pakistan": "Pakistan",
    "india": "India",
    "china": "China",
    "philippines": "Philippines",
    "taiwan": "Taiwan",
    "thailand": "Thailand",
}

# First-level administrative divisions, per country, as lowercase alias -> canonical.
# Cities and districts map to the division that contains them. Only places that occur
# in this corpus are listed: a table of every district in South Asia would be a
# liability, because an entry nobody has checked against the data is a guess.
ADMIN1_ALIASES: dict[str, dict[str, str]] = {
    "Pakistan": {
        "punjab": "Punjab",
        "punjab province": "Punjab",
        "bahawalpur": "Punjab",
        "lodhran": "Punjab",
        "vehari": "Punjab",
        "multan": "Punjab",
        "faisalabad": "Punjab",
        "fsd": "Punjab",  # NIAB, Faisalabad
        "niab": "Punjab",
        "sindh": "Sindh",
        "tando jam": "Sindh",
        "sakrand": "Sindh",
        "islamabad": "Islamabad",
        # Not a Pakistani division. Normalised for spelling and flagged - see the
        # module note.
        "rajasthan": "Rajasthan",
        "rajesthan": "Rajasthan",
    },
    "India": {
        "punjab": "Punjab",
        "punjob": "Punjab",
        "ludhiana": "Punjab",
        "haryana": "Haryana",
        "hissar": "Haryana",
        "sirsa": "Haryana",
        "sahanwala": "Haryana",
        "rajasthan": "Rajasthan",
        "sriganganagar": "Rajasthan",
        "sri-ganganagar": "Rajasthan",
        "delhi": "Delhi",
        "new dehli": "Delhi",
        "new delhi": "Delhi",
        "pusa": "Delhi",
        "uttar pradesh": "Uttar Pradesh",
        "gorakhpur": "Uttar Pradesh",
        "mirzapur": "Uttar Pradesh",
        "chhattisgarh": "Chhattisgarh",
        "bastar": "Chhattisgarh",
        "raipur": "Chhattisgarh",
    },
    "China": {
        "guangdong": "Guangdong",
        "gd": "Guangdong",
        "foshan": "Guangdong",
        "caozhou": "Guangdong",
        "dianbai": "Guangdong",
        "enping": "Guangdong",
        "gaoliang": "Guangdong",
        "huangpu": "Guangdong",
        "jiangmen": "Guangdong",
        "kaiping": "Guangdong",
        "lianzhou": "Guangdong",
        "maonan": "Guangdong",
        "taishan": "Guangdong",
        "tianhe": "Guangdong",
        "yangcheng": "Guangdong",
        "yangchun": "Guangdong",
        "yangdong": "Guangdong",
        "yangxi": "Guangdong",
        "zengcheng": "Guangdong",
        "zhanjiang": "Guangdong",
        "jieyang": "Guangdong",
        "maoming": "Guangdong",
        "shanwei": "Guangdong",
        "xuwen": "Guangdong",
        "yangjiang": "Guangdong",
        "fujian": "Fujian",
        "xiamen": "Fujian",
        "guangxi": "Guangxi",
        "nanning": "Guangxi",
        "hainan": "Hainan",
        "jiangsu": "Jiangsu",
        "xinjiang": "Xinjiang",
        "shihezi": "Xinjiang",
        "wulumuqi": "Xinjiang",
        "yunnan": "Yunnan",
    },
    "Philippines": {
        "davao": "Davao",
        "laguna": "Calabarzon",
    },
    "Taiwan": {"taitung": "Taitung"},
    "Thailand": {"surat thani": "Surat Thani"},
}

# Divisions that are real, but not in the country the record claims.
FOREIGN_ADMIN1: dict[str, set[str]] = {"Pakistan": {"Rajasthan"}}


@dataclass(frozen=True)
class Place:
    """A parsed `/country` string."""

    raw: str
    country: str
    admin1: str = ""
    flags: tuple[str, ...] = ()

    @property
    def resolved(self) -> bool:
        """Whether a first-level division is known.

        Unresolved places are real data and are kept. They just cannot answer the
        question "did this rise *within* one place", because they do not name one.
        """
        return bool(self.admin1)

    @property
    def stratum(self) -> str:
        """The label to group by."""
        if not self.country:
            return self.raw or "unknown"
        return f"{self.country}: {self.admin1}" if self.admin1 else self.country


def _tokens(detail: str) -> list[str]:
    """Split the part after the colon into candidate place names, most specific last.

    Separators seen in this corpus: `,` and `-`, sometimes with spaces, sometimes not
    (`Punjab,Bahawalpur`, `Vehari-Punjab`, `Xinjiang-Shihezi`). A hyphen inside a name
    that we know (`Sri-Ganganagar`) is looked up whole before splitting.
    """
    return [part.strip() for part in re.split(r"[,/]", detail) if part.strip()]


def normalise_location(raw: str) -> Place:
    """Parse one GenBank `/country` value into country + first-level division."""
    text = (raw or "").strip()
    if not text:
        return Place(raw=raw or "", country="", flags=("empty",))

    country_part, _, detail = text.partition(":")
    country = COUNTRY_ALIASES.get(country_part.strip().lower(), country_part.strip())
    if not country:
        return Place(raw=text, country="", flags=("unparsed",))

    table = ADMIN1_ALIASES.get(country, {})
    if not detail.strip():
        return Place(raw=text, country=country, flags=("country-only",))

    for token in _tokens(detail):
        key = token.lower()
        # Whole token first: `Sri-Ganganagar` is a name, `Vehari-Punjab` is two.
        candidates = [key, *[p.strip() for p in key.split("-") if p.strip()]]
        for candidate in candidates:
            admin1 = table.get(candidate)
            if admin1:
                flags = (
                    ("admin1-not-in-country",)
                    if admin1 in FOREIGN_ADMIN1.get(country, set())
                    else ()
                )
                return Place(raw=text, country=country, admin1=admin1, flags=flags)

    # A named place we have no entry for. Unresolved rather than guessed: inventing a
    # division for it would be exactly the silent mistake this module exists to stop.
    return Place(raw=text, country=country, flags=("unmapped-detail",))


def normalise_isolates(
    isolates: Sequence[Isolate],
) -> tuple[list[Isolate], dict[str, Place]]:
    """Rewrite every isolate's `location` to its normalised stratum.

    Returns the isolates and the raw-string -> Place mapping, so a run can print what
    it merged. A normalisation nobody can audit is indistinguishable from a bug.
    """
    mapping = {isolate.location: normalise_location(isolate.location) for isolate in isolates}
    rewritten = [
        Isolate(
            name=isolate.name,
            sequence=isolate.sequence,
            period=isolate.period,
            location=mapping[isolate.location].stratum,
            host=isolate.host,
        )
        for isolate in isolates
    ]
    return rewritten, mapping


def resolved_strata(mapping: dict[str, Place]) -> set[str]:
    """Strata that name a first-level division, and so may confirm a within-place rise."""
    return {place.stratum for place in mapping.values() if place.resolved}
