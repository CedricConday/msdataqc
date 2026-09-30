"""intensitydrift --manifest m.tsv (subject, visit, t1, flair, brainmask[, mask]) [--tolerance 0.10]"""

import argparse
import sys

import pandas as pd

from .core import drift, visit_contrast


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="intensitydrift", description=__doc__)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--tolerance", type=float, default=0.10)
    ap.add_argument("--visits-out")
    a = ap.parse_args(argv)
    m = pd.read_csv(a.manifest, sep="\t")
    rows = []
    from msdataqc.tissue import TissueError

    for _, r in m.iterrows():
        try:
            c = visit_contrast({"T1w": r["t1"], "FLAIR": r["flair"]}, r["brainmask"], r.get("mask") if "mask" in m else None,
                               r.get("tissue") if "tissue" in m else None)
        except TissueError as e:
            print(r["subject"], r["visit"], "SKIPPED:", e, file=sys.stderr, flush=True)
            continue
        for con, v in c.items():
            rows.append({"subject": r["subject"], "visit": r["visit"], "contrast": con, **v})
        print(r["subject"], r["visit"], file=sys.stderr, flush=True)
    t = pd.DataFrame(rows)
    if a.visits_out:
        t.to_csv(a.visits_out, sep="\t", index=False)
    drift(t, a.tolerance).to_csv(sys.stdout, sep="\t", index=False)
    return 0
