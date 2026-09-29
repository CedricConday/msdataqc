"""lesioncount mask.nii.gz  |  lesioncount --manifest m.tsv (subject, mask)"""

import argparse
import json
import sys

import pandas as pd

from .core import grid, spread


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="lesioncount", description=__doc__)
    ap.add_argument("mask", nargs="?")
    ap.add_argument("--manifest")
    a = ap.parse_args(argv)
    if a.manifest:
        rows = []
        for _, r in pd.read_csv(a.manifest, sep="\t").iterrows():
            t = grid(r["mask"])
            for _, x in t.iterrows():
                rows.append({"subject": r["subject"], **x.to_dict()})
            print(r["subject"], file=sys.stderr, flush=True)
        pd.DataFrame(rows).to_csv(sys.stdout, sep="\t", index=False)
        return 0
    t = grid(a.mask)
    print(t.to_string(index=False))
    print(json.dumps(spread(t)))
    return 0
