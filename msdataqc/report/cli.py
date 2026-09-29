"""msqc manifest.tsv --out report.html [--title NAME] [--tolerance 0.10]"""

import argparse

import pandas as pd

from .core import page, run


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="msqc", description=__doc__)
    ap.add_argument("manifest")
    ap.add_argument("--out", required=True)
    ap.add_argument("--title", default="msqc")
    ap.add_argument("--tolerance", type=float, default=0.10)
    a = ap.parse_args(argv)
    res = run(pd.read_csv(a.manifest, sep="\t"), a.tolerance)
    with open(a.out, "w") as fh:
        fh.write(page(res, a.title))
    print(a.out)
    return 0
