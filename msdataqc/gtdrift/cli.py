"""gtdrift mask_t0.nii.gz mask_t1.nii.gz [mask_t2 ...]  |  gtdrift --manifest m.tsv (subject, time, mask)"""

import argparse
import json
import sys

import pandas as pd

from .core import audit


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="gtdrift", description=__doc__)
    ap.add_argument("masks", nargs="*")
    ap.add_argument("--manifest", help="TSV: subject, time, mask (masks of one subject on one grid)")
    ap.add_argument("--merge-mm", type=float, default=2.0)
    ap.add_argument("--min-voxels", type=int, default=3)
    ap.add_argument("--sites-out")
    a = ap.parse_args(argv)
    if a.manifest:
        m = pd.read_csv(a.manifest, sep="\t")
        rows = []
        for sub, g in m.sort_values(["subject", "time"]).groupby("subject", sort=False):
            _, s = audit(list(g["mask"]), list(g["time"]), a.merge_mm, a.min_voxels)
            rows.append({"subject": sub, **{k: v for k, v in s.items() if k != "times"}})
        pd.DataFrame(rows).to_csv(sys.stdout, sep="\t", index=False)
        return 0
    if len(a.masks) < 2:
        ap.error("give at least two masks, or --manifest")
    df, s = audit(a.masks, None, a.merge_mm, a.min_voxels)
    if a.sites_out:
        df.to_csv(a.sites_out, sep="\t", index=False)
    print(json.dumps(s, indent=2))
    return 0
