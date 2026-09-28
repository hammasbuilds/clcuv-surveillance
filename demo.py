"""Show what clcuv-surveillance does, in one command.

    python demo.py

Cotton leaf curl virus genomes from NCBI, aligned, with a variant atlas built
over them - and then the question that decides whether any of it means
anything: **how many of those genomes are independent observations?**

Sequences deposited from one outbreak are near-identical. Counting them as
separate evidence inflates every variant frequency, and the inflation is
largest exactly where sampling was heaviest, which is where people look for
signal. The pipeline therefore runs three readings of the same data:

  pooled               every sequence counts
  stratified           by year and place, so one heavily-sampled site cannot
                       carry a variant on its own
  one per haplotype    clonal duplicates collapsed, so each distinct genome
                       speaks once

The number that survives all three is the number the dataset actually
supports, and it is much smaller than the first one.

Downloads ~250 genomes from NCBI on the first run and caches them, so a second
run is offline. Standard library only - no BLAST, no MAFFT, no API key.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def main() -> int:
    print("clcuv-surveillance: how much of this dataset is independent evidence?", flush=True)
    print(flush=True)
    print("First run downloads genomes from NCBI and caches them under data/;", flush=True)
    print("later runs are offline. Aligning 250 genomes takes roughly two minutes.", flush=True)
    print(flush=True)

    env = {**os.environ, "PYTHONIOENCODING": "utf-8"}
    result = subprocess.run(
        [sys.executable, "scripts/real_data.py", "analyse"],
        cwd=ROOT,
        env=env,
        check=False,
    )
    if result.returncode != 0:
        print(flush=True)
        print("The analysis needs NCBI on the first run. If it is unreachable, the", flush=True)
        print("cached corpus under data/ is used instead - and if that is absent too,", flush=True)
        print("there is nothing to analyse and saying so is the correct outcome.", flush=True)
        return result.returncode

    print(flush=True)
    print("The last two lines are the finding. Variants that look well-supported", flush=True)
    print("on pooled counts are supported by repeated sampling of the same", flush=True)
    print("outbreak, not by repeated observation of the same mutation.", flush=True)
    print(flush=True)
    print("Full write-up: README.md", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
