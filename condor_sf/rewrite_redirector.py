"""Rewrite XRootD redirectors in an infiles JSON.

Why this exists: ~97 files in `Muon1_Run2023C-v4` failed with
`OSError: Failed to open file: [ERROR] Operation expired` when accessed via
the FNAL-local redirector `cmsxrootd.fnal.gov`. The CMS global redirector
(`cms-xrd-global.cern.ch`) routes to any storage site that hosts a replica
and tends to recover from FNAL-local outages.

Default behavior is non-destructive: writes to <input>_global.json next to
the original. Pass --inplace to overwrite.

Usage:
    python3 rewrite_redirector.py infiles/2023/2023_MuonData.json
    python3 rewrite_redirector.py infiles/2023/2023_MuonData.json --inplace
    python3 rewrite_redirector.py infiles/2023/2023_MuonData.json \\
        --redirector cms-xrd-global.cern.ch \\
        --output infiles/2023/2023_MuonData_global.json
"""
from __future__ import annotations
import argparse
import json
import os
import re
import sys

REDIRECTOR_RE = re.compile(r"^root://[^/]+//")


def rewrite_value(val, new_prefix: str):
    """Rewrite a single URL or a list of URLs. Returns (new_value, n_changed)."""
    if isinstance(val, str):
        if REDIRECTOR_RE.match(val):
            new_val = REDIRECTOR_RE.sub(new_prefix, val)
            return new_val, 1 if new_val != val else 0
        return val, 0
    if isinstance(val, list):
        n = 0
        out = []
        for x in val:
            nx, dn = rewrite_value(x, new_prefix)
            n += dn
            out.append(nx)
        return out, n
    if isinstance(val, dict):
        n = 0
        out = {}
        for k, v in val.items():
            nv, dn = rewrite_value(v, new_prefix)
            n += dn
            out[k] = nv
        return out, n
    return val, 0


def main():
    p = argparse.ArgumentParser()
    p.add_argument("input", help="Input JSON file (e.g. infiles/2023/2023_MuonData.json)")
    p.add_argument("--redirector", default="cms-xrd-global.cern.ch",
                   help="New XRootD redirector host (default: cms-xrd-global.cern.ch)")
    p.add_argument("--output", default=None,
                   help="Output path; default <input>_<short>.json where <short> is "
                        "extracted from --redirector (e.g. _global.json).")
    p.add_argument("--inplace", action="store_true",
                   help="Overwrite the input file (a .bak copy is kept).")
    p.add_argument("--dry-run", action="store_true",
                   help="Print the change count but do not write.")
    args = p.parse_args()

    if not os.path.exists(args.input):
        sys.exit(f"ERROR: input {args.input} not found")

    new_prefix = f"root://{args.redirector}//"
    short = "global" if "global" in args.redirector else args.redirector.split(".")[0]
    if args.inplace:
        out_path = args.input
    elif args.output:
        out_path = args.output
    else:
        base, ext = os.path.splitext(args.input)
        out_path = f"{base}_{short}{ext}"

    with open(args.input) as f:
        data = json.load(f)

    new_data, n_changed = rewrite_value(data, new_prefix)
    print(f"Rewrote {n_changed} URL(s) to {new_prefix}")

    if args.dry_run:
        print("(--dry-run; not writing)")
        return
    if n_changed == 0:
        print("Nothing to write; skipping.")
        return

    if args.inplace:
        bak = args.input + ".bak"
        if not os.path.exists(bak):
            os.replace(args.input, bak)
            print(f"  backed up original to {bak}")
        else:
            print(f"  {bak} already exists; not overwriting it")

    with open(out_path, "w") as f:
        json.dump(new_data, f, indent=2)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
