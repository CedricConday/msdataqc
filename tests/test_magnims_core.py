import json

import nibabel as nib
import numpy as np

from msdataqc.magnims.core import scan_bids


def _img(p, zooms, meta):
    im = nib.Nifti1Image(np.zeros((8, 8, 8), np.int16), np.diag([*zooms, 1]))
    im.header.set_zooms(zooms)
    nib.save(im, p)
    p.with_suffix("").with_suffix(".json").write_text(json.dumps(meta))


def test_monitoring_and_comparability(tmp_path):
    base = {"MagneticFieldStrength": 3, "Manufacturer": "Siemens", "ManufacturersModelName": "Prisma",
            "MRAcquisitionType": "3D", "RepetitionTime": 5.0, "EchoTime": 0.39, "InversionTime": 1.8, "FlipAngle": 120,
            "ReceiveCoilName": "HN64", "SliceThickness": 1.0}
    for ses, meta in (("ses-1", base), ("ses-2", {**base, "MagneticFieldStrength": 1.5})):
        d = tmp_path / "sub-01" / ses / "anat"
        d.mkdir(parents=True)
        _img(d / f"sub-01_{ses}_FLAIR.nii.gz", (1.0, 1.0, 1.0), meta)
    r = scan_bids(tmp_path, "monitoring")
    assert all(s["core_complete"] for s in r["sessions"])
    c = r["comparisons"][0]
    assert not c["comparable"] and c["differs"] == ["MagneticFieldStrength"]


def test_2d_flair_fails_diagnosis(tmp_path):
    d = tmp_path / "sub-01" / "anat"
    d.mkdir(parents=True)
    _img(d / "sub-01_FLAIR.nii.gz", (0.9, 0.9, 5.0), {"MagneticFieldStrength": 1.5, "MRAcquisitionType": "2D", "SliceThickness": 5.0})
    s = scan_bids(tmp_path, "diagnosis")["sessions"][0]
    rules = {f["rule"] for f in s["failed"]}
    assert not s["core_complete"] and "2D slice thickness <= 3 mm" in rules and "FLAIR acquired in 3D" in rules
