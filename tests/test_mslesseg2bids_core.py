import pandas as pd

from msdataqc.mslesseg2bids.core import convert


def test_layout(tmp_path):
    s = tmp_path / "src"
    for p, ts in (("P1", ("T1", "T2")), ("P2", ("T1",))):
        for t in ts:
            d = s / "train" / p / t
            d.mkdir(parents=True)
            for m in ("T1", "FLAIR", "MASK", "brainmask"):
                (d / f"{p}_{t}_{m}.nii.gz").write_bytes(b"x")
    (s / "info_dataset").mkdir()
    (s / "info_dataset" / "clinical_data.csv").write_text("﻿Patient;Timepoint;Age;Sex;MS Type;EDSS;Lesion Number;Lesion Volume;\n"
                                                          "P1;T1;28,1;F;SMRR;3,5;;20674.0;18\nP1;T2;28,5;F;SMRR;3,0;;17877.0;18\n")
    r = convert(s, tmp_path / "bids")
    assert r == {"participants": 2, "images": 6, "longitudinal": 1}
    b = tmp_path / "bids"
    assert (b / "sub-P01/ses-T2/anat/sub-P01_ses-T2_FLAIR.nii.gz").exists()
    assert (b / "derivatives/mslesseg-masks/sub-P01/ses-T1/anat/sub-P01_ses-T1_space-MNI152NLin2009aSym_desc-lesion_mask.nii.gz").exists()
    ses = pd.read_csv(b / "sub-P01/sub-P01_sessions.tsv", sep="\t")
    assert ses["edss"].tolist() == [3.5, 3.0]
