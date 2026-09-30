"""lesionmap manifest.tsv --out DIR [--group column] (manifest: subject, mask[, t1, group columns])"""

import argparse
import json

from .core import run


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="lesionmap", description=__doc__)
    ap.add_argument("manifest")
    ap.add_argument("--out", required=True)
    ap.add_argument("--group")
    ap.add_argument("--title", default="lesion probability")
    a = ap.parse_args(argv)
    print(json.dumps(run(a.manifest, a.out, a.group, a.title), indent=2))
    return 0
