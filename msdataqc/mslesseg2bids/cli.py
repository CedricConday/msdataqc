"""mslesseg2bids "MSLesSeg Dataset" bids_out [--split train|test]"""

import argparse
import json

from .core import convert


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="mslesseg2bids", description=__doc__)
    ap.add_argument("src")
    ap.add_argument("dst")
    ap.add_argument("--split", default="train")
    a = ap.parse_args(argv)
    print(json.dumps(convert(a.src, a.dst, a.split)))
    return 0
