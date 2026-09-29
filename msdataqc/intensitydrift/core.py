"""Scale-free contrast per visit, and its drift between visits.

Raw intensities differ between visits for trivial reasons (gain, normalisation), so every measure is a
ratio to the median of normal-appearing white matter (NAWM) of the same image. Tissue comes from
three-class k-means on the T1w inside the brain mask, lesions excluded. Per visit and contrast:
GM/NAWM, CSF/NAWM, lesion/NAWM, and the NAWM coefficient of variation. A lesion segmenter and a
registration both assume these ratios are stable; drift beyond ``tolerance`` (relative) flags the pair.
"""

from __future__ import annotations

import nibabel as nib
import numpy as np
import pandas as pd


def _kmeans3(x, iters=30):
    c = np.percentile(x, [15, 50, 85]).astype(float)
    for _ in range(iters):
        lab = np.argmin(np.abs(x[:, None] - c[None]), axis=1)
        new = np.array([x[lab == k].mean() if (lab == k).any() else c[k] for k in range(3)])
        if np.allclose(new, c):
            break
        c = new
    return lab


def visit_contrast(images: dict, brain, lesion=None, tissue_path=None, gate: bool = True) -> dict:
    """images: {'T1w': path, 'FLAIR': path, ...} on one grid; tissue classes from the gated T1w split
    (msdataqc.tissue) or an external segmentation. Raises TissueError when the split is implausible."""
    from msdataqc.tissue import classes

    bi = nib.load(brain)
    b = np.asarray(bi.dataobj) > 0
    les = np.asarray(nib.load(lesion).dataobj) > 0.5 if lesion else np.zeros_like(b)
    t1 = np.asarray(nib.load(images["T1w"]).dataobj, dtype=np.float32)
    lab = classes(t1, b, les, tuple(float(z) for z in bi.header.get_zooms()[:3]), tissue_path, gate=gate)
    out = {}
    for name, p in images.items():
        d = np.asarray(nib.load(p).dataobj, dtype=np.float32)
        cover = d > 0  # a contrast may cover less of the brain than the T1w; zeros are outside its field
        def med(mask, d=d, cover=cover):
            v = d[mask & cover]
            return float(np.median(v)) if v.size else None
        wm = med(lab == 2)
        if not wm:
            out[name] = {"gm_wm": None, "csf_wm": None, "wm_cv": None, "raw_wm_median": None, "lesion_wm": None,
                         "coverage": float((cover & b).mean())}
            continue
        g, c, lz = med(lab == 1), med(lab == 0), med(les)
        rec = {"gm_wm": g / wm if g else None, "csf_wm": c / wm if c else None,
               "wm_cv": float(np.std(d[(lab == 2) & cover]) / wm), "raw_wm_median": wm,
               "lesion_wm": lz / wm if lz else None, "coverage": float((cover & b).sum() / b.sum())}
        out[name] = rec
    return out


def drift(table: pd.DataFrame, tolerance: float = 0.10) -> pd.DataFrame:
    """table: subject, visit, contrast, gm_wm, csf_wm, lesion_wm, wm_cv (one row per visit and contrast)."""
    rows = []
    for (sub, con), g in table.sort_values("visit").groupby(["subject", "contrast"]):
        g = g.reset_index(drop=True)
        for i in range(1, len(g)):
            a, b = g.loc[0], g.loc[i]
            rec = {"subject": sub, "contrast": con, "from": a["visit"], "to": b["visit"]}
            worst = 0.0
            for k in ("gm_wm", "csf_wm", "lesion_wm", "wm_cv"):
                if pd.notna(a[k]) and pd.notna(b[k]) and a[k]:
                    r = (b[k] - a[k]) / abs(a[k])
                    rec[f"d_{k}"] = float(r)
                    if k != "wm_cv":
                        worst = max(worst, abs(r))
            rec["max_rel_drift"] = worst
            rec["flag"] = worst > tolerance
            rows.append(rec)
    return pd.DataFrame(rows)
