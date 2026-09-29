"""Run the longitudinal dataset checks and write one page.

Manifest: subject, visit, t1, flair, brainmask[, mask]. Visits of one subject must share a grid for the
lesion-flicker audit (MNI-space datasets do; otherwise register first, e.g. with lesiontrack) and the audit
is skipped with a note when they do not.
"""

from __future__ import annotations

import html

import nibabel as nib
import numpy as np
import pandas as pd


def run(manifest: pd.DataFrame, tolerance: float = 0.10, gate: bool = True) -> dict:
    from msdataqc.gtdrift.core import audit
    from msdataqc.intensitydrift.core import drift, visit_contrast
    from msdataqc.tissue import TissueError

    rows, burden, flick, refused = [], [], [], []
    for _, r in manifest.iterrows():
        mask = r["mask"] if "mask" in manifest and pd.notna(r.get("mask")) else None
        try:
            vc = visit_contrast({"T1w": r["t1"], "FLAIR": r["flair"]}, r["brainmask"], mask,
                                r["tissue"] if "tissue" in manifest and pd.notna(r.get("tissue")) else None, gate=gate)
        except TissueError as e:
            refused.append({"subject": r["subject"], "visit": r["visit"], "reason": str(e)})
            vc = {}
        for con, v in vc.items():
            rows.append({"subject": r["subject"], "visit": r["visit"], "contrast": con, **v})
        if mask:
            im = nib.load(mask)
            burden.append({"subject": r["subject"], "visit": r["visit"],
                           "lesion_ml": float((np.asarray(im.dataobj) > 0.5).sum() * np.prod(im.header.get_zooms()[:3]) / 1000)})
    vis = pd.DataFrame(rows)
    dr = drift(vis, tolerance) if len(vis) else pd.DataFrame()
    if "mask" in manifest:
        for sub, g in manifest.sort_values("visit").groupby("subject"):
            if len(g) < 2:
                continue
            try:
                _, s = audit(list(g["mask"]))
                flick.append({"subject": sub, **{k: v for k, v in s.items() if k != "times"}})
            except ValueError as e:
                flick.append({"subject": sub, "note": str(e)})
    return {"visits": vis, "drift": dr, "burden": pd.DataFrame(burden), "flicker": pd.DataFrame(flick),
            "refused": pd.DataFrame(refused, columns=["subject", "visit", "reason"])}


def _table(df: pd.DataFrame, cols) -> str:
    cols = [c for c in cols if c in df]
    head = "".join(f"<th>{html.escape(c)}</th>" for c in cols)
    body = "".join("<tr>" + "".join(f"<td>{html.escape(f'{v:.3g}' if isinstance(v, float) else str(v))}</td>" for v in r) + "</tr>"
                   for r in df[cols].itertuples(index=False))
    return f"<table><tr>{head}</tr>{body}</table>"


def page(res: dict, title: str = "msqc") -> str:
    dr, fl, bu = res["drift"], res["flicker"], res["burden"]
    lines = []
    if len(dr):
        for con in sorted(dr["contrast"].unique()):
            x = dr[dr["contrast"] == con]
            lines.append(f"{con}: {int(x['flag'].sum())} of {len(x)} visit pairs change contrast beyond tolerance; "
                         f"median grey/white change {x['d_gm_wm'].abs().median() * 100:.0f} %.")
    if len(fl) and "sites" in fl:
        lines.append(f"Lesion ground truth: {int(fl['implausible_sites'].sum())} of {int(fl['sites'].sum())} lesion sites flicker, "
                     f"appear transiently or vanish across visits.")
    rf = res.get("refused", pd.DataFrame())
    if len(rf):
        lines.insert(0, f"T1w tissue split refused for {len(rf)} visits (grey and white matter do not separate); "
                        "their contrast is not measured. Give an external tissue segmentation for them.")
    summary = "".join(f"<li>{html.escape(s)}</li>" for s in lines) or "<li>nothing to report</li>"
    return (f"<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>"
            f"<title>{html.escape(title)}</title><style>body{{font-family:system-ui,sans-serif;max-width:1000px;margin:24px auto;padding:0 16px}}"
            f"table{{border-collapse:collapse;font-size:.85rem;margin-bottom:24px}}td,th{{border-bottom:1px solid #ddd;padding:4px 8px;text-align:left}}"
            f"</style></head><body><h1>{html.escape(title)}</h1><ul>{summary}</ul>"
            f"<h2>Contrast drift (intensitydrift)</h2>{_table(dr, ['subject', 'contrast', 'from', 'to', 'd_gm_wm', 'd_lesion_wm', 'flag'])}"
            f"<h2>Lesion ground-truth consistency (gtdrift)</h2>{_table(fl, ['subject', 'timepoints', 'sites', 'flicker', 'transient', 'vanished', 'note'])}"
            f"<h2>Lesion burden per visit</h2>{_table(bu, ['subject', 'visit', 'lesion_ml'])}"
            f"<h2>Visits refused by the tissue gate</h2>{_table(rf, ['subject', 'visit', 'reason'])}</body></html>")
