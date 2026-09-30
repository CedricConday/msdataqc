"""Temporal consistency of a series of co-registered lesion masks.

Sites are the connected components of the union of all timepoints, with pieces whose gap is at most
``merge_mm`` (default 2 mm, one empty voxel at 1 mm spacing) joined, so a lesion that shifts by a voxel stays
one site; ``merge_mm=0`` joins nothing. For each site the presence pattern over time is
recorded (a timepoint counts as present when at least ``min_voxels`` of its mask fall in the site).

Patterns a chronic T2 lesion should not show:
* **flicker**: present, then absent, then present again (``1..0..1``);
* **transient**: absent at the first visit, present in between, absent at the last (``0..1..0``);
* **vanished**: present at the first visit and absent at the last.
Real lesions can shrink and, rarely, disappear on T2; the counts are an audit of the annotation,
not a claim about biology. The volume-weighted share tells how much of the ground truth is involved.
"""

from __future__ import annotations

import re

import nibabel as nib
import numpy as np
import pandas as pd
from scipy import ndimage as ndi


def audit(masks: list, times: list | None = None, merge_mm: float = 2.0, min_voxels: int = 3) -> tuple[pd.DataFrame, dict]:
    imgs = [nib.load(m) for m in masks]
    ref = imgs[0]
    for im in imgs[1:]:
        if im.shape != ref.shape or not np.allclose(im.affine, ref.affine, atol=1e-3):
            raise ValueError("all masks must share one grid (register them first)")
    zooms = tuple(float(z) for z in ref.header.get_zooms()[:3])
    vox = float(np.prod(zooms))
    stack = np.stack([np.asarray(im.dataobj) > 0.5 for im in imgs])
    union = stack.any(axis=0)
    st = ndi.generate_binary_structure(3, 2)
    if merge_mm > 0:
        # voxels of the union closer than merge_mm (millimetres, any voxel spacing) belong to one site
        grown = ndi.distance_transform_edt(~union, sampling=zooms) <= merge_mm / 2
        sites, n = ndi.label(grown, st)
        sites = np.where(union, sites, 0)
    else:
        sites, n = ndi.label(union, st)
    t = len(masks)
    counts = np.zeros((n + 1, t), np.int64)
    flat = sites.ravel()
    for k in range(t):
        counts[:, k] = np.bincount(flat[stack[k].ravel()], minlength=n + 1)
    rows = []
    for s in range(1, n + 1):
        pres = counts[s] >= min_voxels
        if not pres.any():
            continue
        pat = "".join("1" if x else "0" for x in pres)
        flick = bool(re.search("10+1", pat))
        trans = pat[0] == "0" and pat[-1] == "0"
        van = pat[0] == "1" and pat[-1] == "0"
        rows.append({"site": s, "pattern": pat, "max_volume_mm3": float(counts[s].max() * vox),
                     **{f"volume_t{k}_mm3": float(counts[s, k] * vox) for k in range(t)},
                     "flicker": flick, "transient": trans, "vanished": van})
    df = pd.DataFrame(rows)
    vol = df["max_volume_mm3"].sum() if len(df) else 0.0
    def share(col):
        return float(df.loc[df[col], "max_volume_mm3"].sum() / vol) if vol else 0.0
    summary = {"timepoints": t, "times": times, "sites": len(df),
               "flicker": int(df["flicker"].sum()) if len(df) else 0,
               "transient": int(df["transient"].sum()) if len(df) else 0,
               "vanished": int(df["vanished"].sum()) if len(df) else 0,
               "flicker_volume_share": share("flicker") if len(df) else 0.0,
               "vanished_volume_share": share("vanished") if len(df) else 0.0,
               "implausible_sites": int((df["flicker"] | df["transient"] | df["vanished"]).sum()) if len(df) else 0}
    summary["implausible_share"] = summary["implausible_sites"] / summary["sites"] if summary["sites"] else 0.0
    return df, summary
