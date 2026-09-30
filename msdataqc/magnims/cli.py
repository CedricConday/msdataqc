"""magnims-check /path/to/bids [--purpose diagnosis|monitoring] [--json report.json]"""

import argparse
import json
import sys

from .core import scan_bids


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="magnims-check", description=__doc__)
    ap.add_argument("bids")
    ap.add_argument("--purpose", choices=("diagnosis", "monitoring"), default="diagnosis")
    ap.add_argument("--json")
    a = ap.parse_args(argv)
    r = scan_bids(a.bids, a.purpose)
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(r, fh, indent=2, default=str)
    for s in r["sessions"]:
        state = "complete" if s["core_complete"] else f"{s['core_failed']} core failed, {s['core_unknown']} unknown"
        print(f"{s['subject']}/{s['session']}: {state}; kinds {','.join(s['kinds'])}")
        for f in s["failed"]:
            print(f"    - {f['rule']}: {f['detail']} ({f['image']})")
    for c in r["comparisons"]:
        print(f"{c['subject']} {c['from']} -> {c['to']}: {'comparable' if c['comparable'] else 'NOT comparable: ' + ', '.join(c['differs'])}")
    ok = all(s["core_complete"] for s in r["sessions"]) and all(c["comparable"] for c in r["comparisons"])
    print(f"sessions {len(r['sessions'])}, core complete {sum(s['core_complete'] for s in r['sessions'])}", file=sys.stderr)
    return 0 if ok else 1
