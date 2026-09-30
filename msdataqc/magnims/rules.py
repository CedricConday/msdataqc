"""The rules, as read from Wattjes et al., Lancet Neurol 2021;20:653-670 (2021 MAGNIMS-CMSC-NAIMS
consensus on MRI in MS), brain protocol. Each rule carries its level: ``core`` (recommended for every
scan of that purpose) or ``optional``. They are the author's reading, kept in one file so a reader
can check each against the paper and change it.
"""

FIELD_MIN_T = 1.5          # minimum field strength; 3 T preferred
FIELD_PREFERRED_T = 3.0
MAX_2D_SLICE_MM = 3.0      # 2D sequences: slice thickness <= 3 mm, no gap
MAX_3D_VOXEL_MM = 1.0      # 3D sequences: (near-)isotropic ~1 mm
ISO_TOLERANCE = 1.25       # 3D voxel counts as isotropic when max/min edge <= this
MAX_INPLANE_MM = 1.0
POST_CONTRAST_DELAY_MIN = 5.0

# purpose -> list of (key, label, level)
PROTOCOL = {
    "diagnosis": [
        ("flair3d", "3D T2-FLAIR (sagittal acquisition, ~1 mm isotropic)", "core"),
        ("t2", "axial T2 or PD/T2 (2D <= 3 mm, or 3D)", "core"),
        ("t1post", "post-contrast T1 (at least 5 min after gadolinium)", "core"),
        ("t1pre", "pre-contrast T1 (3D preferred, also for atrophy)", "optional"),
        ("dwi", "diffusion-weighted imaging", "optional"),
        ("swi", "susceptibility-based (SWI/T2*) for the central vein sign or paramagnetic rim", "optional"),
    ],
    "monitoring": [
        ("flair3d", "3D T2-FLAIR, same protocol as the reference scan", "core"),
        ("t2", "axial T2 or PD/T2", "optional"),
        ("t1post", "post-contrast T1 (when inflammatory activity must be excluded)", "optional"),
        ("t1pre", "3D T1 (for atrophy)", "optional"),
    ],
}

# keys that must match between visits for longitudinal measurement
COMPARABLE = {"MagneticFieldStrength": "exact", "Manufacturer": "exact", "ManufacturersModelName": "exact",
              "MRAcquisitionType": "exact", "voxel_mm": 0.10, "SliceThickness": 0.10,
              "RepetitionTime": 0.10, "EchoTime": 0.10, "InversionTime": 0.10, "FlipAngle": 0.10,
              "ReceiveCoilName": "exact", "SoftwareVersions": "warn"}
