# clcuv-surveillance

Genomic surveillance for Cotton Leaf Curl Virus, the most serious threat to Punjab's cotton crop, on real NCBI
GenBank genomes. We built it for plant virologists, breeders and diagnostic labs who need to know which variant is
rising, whether anything is selecting for it, how the strains are related, and whether a genome is something never
seen before, without mistaking sampling bias for a rising variant. Pure Python, zero dependencies.

## What it does

- **Parses GenBank records** and normalises the fields that make them surveillance data: organism, place
  (`Pakistan: Punjab`, `Pakistan: Punjab province` and `Pakistan: Punjab,Bahawalpur` are one stratum) and collection
  date (`2019`, `May-2019`, `01-May-2019`, `2015-01`).
- **Aligns the genomes** with banded centre-star Needleman-Wunsch; refuses differently rotated circular genomes,
  whose naive alignment looks fine and means nothing. The alignment of the committed corpus is cached and verified
  row by row against the sequences.
- **Builds a mutation atlas** and flags variants whose frequency is rising, with three controls applied in turn:
  1. a two-proportion z-test beside the effect-size threshold, and periods with too few sequences excluded;
  2. `stratify=True`: the rise is re-tested inside each location, so a variant common in a newly sampled region
     does not appear to emerge;
  3. `collapse_clonal`: each (year, location) is reduced to one sequence per haplotype, so one submission batch of
     identical genomes is not counted as many independent observations.
- **Measures selection pressure** with Nei-Gojobori dN/dS; with no synonymous differences the ratio is `None`, not
  infinity.
- **Builds trees** with neighbour-joining (and UPGMA, which assumes a molecular clock) on Jukes-Cantor distances;
  at p >= 0.75 the distance is `inf`.
- **Classifies strains with k-mer profiles** and refuses to assign a genome beyond each strain's calibrated distance
  threshold: `strain=None, novel=True` is a finding. A genome close to two strains is flagged.
- **Detects recombination** as a switch in similarity between two named parents at a coordinate.
- **Exports one sequence per haplotype** as aligned FASTA for
  [primer-designer](https://github.com/hammasbuilds/primer-designer), which checks primer conservation on it.

## Sample input / output

Input: `data/clcuv.gb`, 254 GenBank records as NCBI serves them (committed), for example:

```
LOCUS       PZ547674                2739 bp    DNA     linear   VRL 02-SEP-2026
  ORGANISM  Cotton leaf curl Multan virus
     source          1..2739
                     /host="cotton"
                     /geo_loc_name="India"
                     /collection_date="21-Aug-2022"
```

`python demo.py` (real output, trimmed in the middle):

```
254 records parsed
  with country : 248
  with date    : 229

250 are Cotton leaf curl Multan virus; the rest are other CLCuV species and excluded

alignment cached in clcuv_aligned.fasta.gz, reusing it
  0.1s
  {"sequences": 250, "width": 3162, "all_gap_columns": 0, "invariant_columns": 962, "invariant_fraction": 0.3042, "gap_fraction": 0.1335}
  consensus: 3162 bp, 449 uncalled

228 have both a year and a place, and can be surveilled

--- how many independent genomes are actually here? ---
  ...
  2019  Pakistan: Punjab                8 sequences ->  3 haplotypes  (x2.67)   <- too clonal to test
  2020  Pakistan: Punjab               13 sequences -> 11 haplotypes  (x1.18)
  2021  China: Guangdong               25 sequences -> 21 haplotypes  (x1.19)
  2021  Pakistan: Punjab                8 sequences ->  1 haplotype   (x8.0)   <- too clonal to test
  2021  Pakistan: Sindh                 5 sequences ->  1 haplotype   (x5.0)   <- too clonal to test
  ...

837 variants above 1% against the consensus

--- what survives each control ---
  all sequences        n=228 pooled=92  stratified=76
  one per haplotype    n=184 pooled=0   stratified=0
```

## Results

On the 250 public Cotton leaf curl Multan virus genomes (254 records, 228 with both a year and a place):

| Control | Variants flagged as emerging |
|---|---:|
| Effect size + z-test, pooled over all sequences | 92 |
| ... re-tested inside each location | 76 |
| ... with each clonal group counted once (184 haplotypes) | **0** |

- 228 isolates are 184 haplotypes across 67 (year, location) strata; the largest clonal group is 8 genomes (Punjab
  2021, all 28 pairs identical).
- Merging location spellings removed 16 of the 92 pooled candidates; collapsing clonal batches removed the other 76.
- No variant in the public CLCuMuV genomes can be shown to be emerging: there are not enough independent genomes
  from any single place and pair of years. The CLI prints the count across a sweep of `min_samples` values, and every
  setting ends at zero after the clonal collapse.
- Alignment: 3,162 columns, 30.4% invariant. The export to primer-designer is 250 records, 210 haplotypes.
- Timing on an 8-core laptop under load: `python demo.py` 15-31 s with the committed alignment; a full realignment
  (`--no-cache`) 7-13 minutes, about 1.6 s per genome pair.

## Tech stack

- Python 3.10+, standard library only (no Biopython, BLAST or MAFFT)
- Own implementations of banded Needleman-Wunsch, Jukes-Cantor, neighbour-joining, UPGMA, Nei-Gojobori dN/dS,
  k-mer classification and a recombination scan
- NCBI E-utilities over `urllib` for `scripts/real_data.py`; `argparse` CLI (`clcuv`)
- pytest and ruff for development, `uv` for the scripts

## Run it

```bash
git clone https://github.com/hammasbuilds/clcuv-surveillance
cd clcuv-surveillance
pip install -e ".[dev]"  # zero runtime dependencies; pytest/ruff for development
pytest -q                # 231 tests, no sequence download
python demo.py           # the whole finding on the committed corpus, offline
```

Re-download from NCBI and rerun everything (cached after the first run):

```bash
uv run python scripts/real_data.py analyse              # add --no-cache to realign from scratch
uv run python scripts/export_alignment.py ../primer-designer/data/clcuv_aligned.fasta
```

On your own genomes:

```bash
clcuv analyse --fasta isolates.fasta --metadata isolates.csv [--json]
clcuv export  --fasta isolates.fasta --metadata isolates.csv --out deduped.fasta
```

`isolates.csv` needs a name column (`name`, `accession`, `id`, ...), a period (`period`, `year`, `date`, ...) and a
location (`location`, `country`, `region`, ...) for every sequence. Unaligned FASTA is aligned first: 70 genomes of
2.7 kb took 1 m 58 s at the default `--band 120` and 57 s at `--band 40` (a band narrower than an indel cannot place
it). Already-aligned input (equal length, gaps present) skips the aligner. The aligner is meant for closely related,
co-oriented, similar-length genomes such as begomovirus submissions; use MAFFT for divergent sets.

From Python:

```python
from clcuv import Isolate, build_atlas, emerging_variants, selection_pressure
from clcuv import distance_matrix, neighbour_joining, StrainClassifier

isolates = [Isolate(name, seq, period="2026-Q1", location="Multan") for name, seq in ...]

atlas = build_atlas(isolates, reference)
for e in emerging_variants(atlas):
    print(e.summary())      # e.g. C21G: 1.0% (2024-Q1) -> 42.0% (2026-Q1), z=7.06

selection_pressure(coat_protein_a, coat_protein_b).interpretation()

names, matrix = distance_matrix(sequences, names=labels)
print(neighbour_joining(names, matrix).newick())

clf = StrainClassifier().fit(labelled_genomes)
clf.classify(new_genome)    # strain=None means "I have not seen this before"
```

## Licence

MIT. The genomes are public GenBank records (NCBI).
