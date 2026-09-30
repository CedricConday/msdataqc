"""MSLesSeg layout -> BIDS.

Input: the unzipped ``MSLesSeg Dataset`` folder (train/P*/T*/P*_T*_{T1,T2,FLAIR,MASK,brainmask}.nii.gz and
info_dataset/clinical_data.csv). Output:

    sub-P01/ses-T1/anat/sub-P01_ses-T1_{T1w,T2w,FLAIR}.nii.gz (+ JSON sidecar noting the MNI preprocessing)
    derivatives/mslesseg-masks/sub-P01/ses-T1/anat/sub-P01_ses-T1_space-MNI152NLin2009aSym_desc-{lesion,brain}_mask.nii.gz
    participants.tsv (sex, course), sub-*/sub-*_sessions.tsv (age, EDSS, lesion volume)

Files are hard-linked when possible (no copy), else copied. The images in MSLesSeg are already
co-registered to MNI space and skull-stripped; the sidecars say so, since BIDS raw data are normally
unprocessed.
"""

from __future__ import annotations

import csv
import json
import os
import shutil
from pathlib import Path

import pandas as pd

SUFFIX = {"T1": "T1w", "T2": "T2w", "FLAIR": "FLAIR"}


def _link(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        return
    try:
        os.link(src, dst)
    except OSError:
        shutil.copy2(src, dst)


def _num(x):
    x = (x or "").strip().replace(",", ".")
    try:
        return float(x)
    except ValueError:
        return None


def convert(src, dst, split: str = "train") -> dict:
    src, dst = Path(src), Path(dst)
    root = src / split
    clin = {}
    cf = src / "info_dataset" / "clinical_data.csv"
    if cf.exists():
        with open(cf, encoding="utf-8-sig") as fh:
            for r in list(csv.reader(fh, delimiter=";"))[1:]:
                if len(r) >= 6:
                    clin[(r[0], r[1])] = {"age": _num(r[2]), "sex": r[3], "course": r[4], "edss": _num(r[5]),
                                          "lesion_volume_mm3": _num(r[7]) if len(r) > 7 else None}
    parts, n_img = [], 0
    for pdir in sorted(root.glob("P*"), key=lambda p: int(p.name[1:])):
        sub = f"sub-P{int(pdir.name[1:]):02d}"
        sess = []
        for tdir in sorted(pdir.glob("T*")):
            ses = f"ses-{tdir.name}"
            for mod, suf in SUFFIX.items():
                f = tdir / f"{pdir.name}_{tdir.name}_{mod}.nii.gz"
                if f.exists():
                    out = dst / sub / ses / "anat" / f"{sub}_{ses}_{suf}.nii.gz"
                    _link(f, out)
                    out.with_suffix("").with_suffix(".json").write_text(json.dumps(
                        {"Description": "MSLesSeg preprocessed image: co-registered to MNI space, skull-stripped", "SourceFile": f.name}, indent=2))
                    n_img += 1
            for mod, desc in (("MASK", "lesion"), ("brainmask", "brain")):
                f = tdir / f"{pdir.name}_{tdir.name}_{mod}.nii.gz"
                if f.exists():
                    _link(f, dst / "derivatives" / "mslesseg-masks" / sub / ses / "anat" /
                          f"{sub}_{ses}_space-MNI152NLin2009aSym_desc-{desc}_mask.nii.gz")
            c = clin.get((pdir.name, tdir.name), {})
            sess.append({"session_id": ses, "age": c.get("age"), "edss": c.get("edss"), "lesion_volume_mm3": c.get("lesion_volume_mm3")})
        if sess:
            pd.DataFrame(sess).to_csv(dst / sub / f"{sub}_sessions.tsv", sep="\t", index=False, na_rep="n/a")
            c0 = clin.get((pdir.name, "T1"), {})
            parts.append({"participant_id": sub, "sex": c0.get("sex"), "course": c0.get("course"), "sessions": len(sess)})
    pd.DataFrame(parts).to_csv(dst / "participants.tsv", sep="\t", index=False, na_rep="n/a")
    (dst / "dataset_description.json").write_text(json.dumps({"Name": "MSLesSeg (BIDS conversion)", "BIDSVersion": "1.9.0",
        "DatasetType": "raw", "GeneratedBy": [{"Name": "mslesseg2bids", "Version": "0.1.0"}],
        "HowToAcknowledge": "Cite the MSLesSeg dataset paper (Guarnera et al.)"}, indent=2))
    dd = dst / "derivatives" / "mslesseg-masks"
    dd.mkdir(parents=True, exist_ok=True)
    (dd / "dataset_description.json").write_text(json.dumps({"Name": "MSLesSeg expert lesion and brain masks", "BIDSVersion": "1.9.0",
        "DatasetType": "derivative", "GeneratedBy": [{"Name": "MSLesSeg annotators"}]}, indent=2))
    return {"participants": len(parts), "images": n_img, "longitudinal": sum(p["sessions"] > 1 for p in parts)}
