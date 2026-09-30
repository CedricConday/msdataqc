"""Classify the anatomical images of each BIDS session, test them against the consensus, compare visits."""

from __future__ import annotations

import json
from itertools import pairwise
from pathlib import Path

import nibabel as nib
import pandas as pd

from . import rules as R


def _meta(nii: Path) -> dict:
    js = Path(str(nii).replace(".nii.gz", ".json").replace(".nii", ".json"))
    d = json.loads(js.read_text()) if js.exists() else {}
    img = nib.load(nii)
    z = [float(x) for x in img.header.get_zooms()[:3]]
    d["voxel_mm"] = z
    d["_file"] = str(nii)
    d["_has_sidecar"] = js.exists()
    return d


def classify(nii: Path, d: dict) -> str | None:
    name = nii.name
    desc = (d.get("SeriesDescription", "") + " " + d.get("ProtocolName", "")).lower()
    three_d = d.get("MRAcquisitionType") == "3D" or (max(d["voxel_mm"]) / min(d["voxel_mm"]) <= R.ISO_TOLERANCE
                                                     and max(d["voxel_mm"]) <= 1.3)
    if name.endswith(("_FLAIR.nii.gz", "_FLAIR.nii")) or "flair" in desc:
        return "flair3d" if three_d else "flair2d"
    if "_T1w" in name:
        post = "ce-" in name or "gad" in desc or "post" in desc or d.get("ContrastBolusIngredient")
        return "t1post" if post else "t1pre"
    if "_T2w" in name or "_PDT2" in name or "_PD" in name:
        return "t2"
    if "_dwi" in name:
        return "dwi"
    if "_T2starw" in name or "_swi" in name.lower() or "swi" in desc:
        return "swi"
    return None


def check_image(kind: str, d: dict) -> list[dict]:
    out = []
    def add(rule, ok, detail, level="core"):
        out.append({"image": Path(d["_file"]).name, "kind": kind, "rule": rule, "pass": ok, "level": level, "detail": detail})
    f = d.get("MagneticFieldStrength")
    add("field strength >= 1.5 T", None if f is None else f >= R.FIELD_MIN_T, f"{f} T" if f else "not in sidecar")
    if f is not None and f < R.FIELD_PREFERRED_T:
        add("3 T preferred", False, f"{f} T", "optional")
    z = d["voxel_mm"]
    if kind in ("flair3d", "t1pre", "t1post") and (d.get("MRAcquisitionType") == "3D" or kind == "flair3d"):
        iso = max(z) / min(z) <= R.ISO_TOLERANCE
        add("3D voxel ~1 mm isotropic", iso and max(z) <= R.MAX_3D_VOXEL_MM + 0.05, "x".join(f"{v:.2f}" for v in z))
    else:
        st = d.get("SliceThickness", max(z))
        add("2D slice thickness <= 3 mm", st <= R.MAX_2D_SLICE_MM + 1e-6, f"{st:.2f} mm")
        add("in-plane <= 1 mm", sorted(z)[1] <= R.MAX_INPLANE_MM + 0.05, "x".join(f"{v:.2f}" for v in z), "optional")
    if kind == "flair2d":
        add("FLAIR acquired in 3D", False, "2D FLAIR; the consensus recommends 3D", "core")
    if kind == "t1post":
        delay = d.get("ContrastDelayMinutes")
        add("post-contrast delay >= 5 min", None if delay is None else delay >= R.POST_CONTRAST_DELAY_MIN,
            f"{delay} min" if delay is not None else "delay not recorded")
    if not d["_has_sidecar"]:
        add("JSON sidecar present", False, "without it field strength and sequence cannot be checked", "core")
    return out


def check_session(files: list[Path], purpose: str = "diagnosis") -> tuple[pd.DataFrame, dict]:
    rows, have = [], {}
    for f in files:
        d = _meta(f)
        k = classify(f, d)
        if k is None:
            continue
        have.setdefault(k, []).append(d)
        rows += check_image(k, d)
    present = []
    for key, lab, level in R.PROTOCOL[purpose]:
        ok = key in have
        present.append({"image": "-", "kind": key, "rule": f"{level}: {lab}", "pass": ok, "level": level,
                        "detail": "present" if ok else ("2D FLAIR only" if key == "flair3d" and "flair2d" in have else "missing")})
    df = pd.DataFrame(present + rows)
    core = df[df["level"] == "core"]
    summary = {"purpose": purpose, "images": sum(len(v) for v in have.values()), "kinds": sorted(have),
               "core_failed": int((core["pass"] == False).sum()),
               "core_unknown": int(core["pass"].isna().sum()),
               "core_complete": bool((core["pass"] == True).all())}
    return df, summary


def compare(a: dict, b: dict) -> list[dict]:
    out = []
    for key, how in R.COMPARABLE.items():
        va, vb = a.get(key), b.get(key)
        if va is None or vb is None:
            out.append({"key": key, "a": va, "b": vb, "match": None})
            continue
        if how in ("exact", "warn"):
            m = va == vb
        else:
            if isinstance(va, list):
                m = all(abs(x - y) <= how * max(abs(x), 1e-9) for x, y in zip(va, vb))
            else:
                m = abs(va - vb) <= how * max(abs(va), 1e-9)
        out.append({"key": key, "a": va, "b": vb, "match": m, "severity": "warn" if how == "warn" else "break"})
    return out


def scan_bids(root: Path, purpose: str = "diagnosis") -> dict:
    root = Path(root)
    report = {"sessions": [], "comparisons": []}
    subs = sorted(p for p in root.glob("sub-*") if p.is_dir())
    for s in subs:
        sess = sorted(p for p in s.glob("ses-*") if p.is_dir()) or [s]
        flairs = []
        for ss in sess:
            files = sorted((ss / "anat").glob("*.nii*")) if (ss / "anat").exists() else []
            if not files:
                continue
            df, summ = check_session(files, purpose)
            report["sessions"].append({"subject": s.name, "session": ss.name, **summ,
                                       "failed": df[df["pass"] == False][["image", "rule", "detail"]].to_dict("records")})
            fl = [f for f in files if "FLAIR" in f.name]
            if fl:
                flairs.append((ss.name, _meta(fl[0])))
        for (n0, d0), (n1, d1) in pairwise(flairs):
            c = compare(d0, d1)
            breaks = [x["key"] for x in c if x["match"] is False and x.get("severity") == "break"]
            report["comparisons"].append({"subject": s.name, "from": n0, "to": n1, "comparable": not breaks,
                                          "differs": breaks, "detail": c})
    return report
