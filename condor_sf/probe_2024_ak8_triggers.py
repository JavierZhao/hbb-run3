"""One-off probe: list HLT branches in 2024 Muon0/MC files that LOOK LIKE
HLT_AK8PFJet425_SoftDropMass40 (which is missing in 2024 NanoAODv15 per
the quick-test scan on 2026-05-29).

Open one file per dataset, scan Events.keys() for any HLT branch matching
several candidate patterns (varying jet pT thresholds, mass-window cuts,
PNet tags, etc.), then summarize unique matches across all probed files.

Run via shell/condor with the coffea image so uproot is available.
"""
from __future__ import annotations
import json
import os
import re
import sys
import time

import uproot

print = lambda *a, _p=print, **kw: _p(*a, **{**kw, "flush": True})

INFILES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "infiles", "2024")

# One representative file per dataset is enough — HLT menu is constant within an era.
PROBE_KEYS = {
    "2024_MuonData.json": ["Muon0_Run2024C", "Muon0_Run2024F", "Muon0_Run2024H", "Muon0_Run2024I"],
    "2024_ttbar.json":    ["TTtoLNu2Q"],
    "2024_QCD_MuEnriched.json": ["QCD_MuEnriched_PT_300to470"],
}

# Patterns we want to surface as POTENTIAL stand-ins for AK8PFJet425_SoftDropMass40.
# The pT threshold may have moved (400/420/425/450), the SoftDropMass cut may be
# different (30/40/50), and the branch may carry a tag suffix (_PNetBB*, _PNetTau*).
SEARCH_REGEXES = [
    re.compile(r"^HLT_AK8PFJet\d+_SoftDropMass\d+(?:_.*)?$"),
    re.compile(r"^HLT_AK8PFJet\d+_PNetBB.*$"),
    re.compile(r"^HLT_AK8PFJet\d+_Mass.*$"),
    re.compile(r"^HLT_AK8DiPFJet.*SoftDropMass.*$"),
    # Catch-all: any AK8 trigger with a softdrop or mass qualifier
    re.compile(r"^HLT_AK8.*[Ss]oftDrop.*$"),
]

# Exact name we KNOW is missing — sanity check that it is in fact absent
KNOWN_MISSING = "HLT_AK8PFJet425_SoftDropMass40"


def probe(path, label):
    try:
        f = uproot.open(path, timeout=60)
    except Exception as exc:
        print(f"  [{label}] open failed: {exc}")
        return set()
    try:
        events = f["Events"]
    except Exception as exc:
        print(f"  [{label}] no Events tree: {exc}")
        f.close()
        return set()
    try:
        keys = set(events.keys())
    except Exception as exc:
        print(f"  [{label}] keys() failed: {exc}")
        f.close()
        return set()
    finally:
        try: f.close()
        except Exception: pass

    matches = set()
    for r in SEARCH_REGEXES:
        matches.update(k for k in keys if r.match(k))
    print(f"  [{label}] {len(keys)} total branches; {len(matches)} AK8/SoftDrop matches; "
          f"AK8PFJet425_SoftDropMass40 present? {KNOWN_MISSING in keys}")
    return matches


def main():
    all_matches = {}
    for json_name, keys in PROBE_KEYS.items():
        ds_path = os.path.join(INFILES_DIR, json_name)
        if not os.path.exists(ds_path):
            print(f"missing {ds_path}, skipping"); continue
        with open(ds_path) as f:
            d = json.load(f)
        print(f"\n=== {json_name} ===")
        for key in keys:
            files = d.get(key, [])
            if not files:
                print(f"  [{key}] no files in JSON"); continue
            sample = files[0]
            print(f"  [{key}] probing 1 of {len(files)} files: {sample}")
            t0 = time.time()
            m = probe(sample, key)
            print(f"    ({time.time()-t0:.1f}s)")
            for x in m:
                all_matches.setdefault(x, []).append(key)

    print("\n" + "=" * 60)
    print("UNIQUE AK8 / SoftDrop-style HLT BRANCHES SEEN IN 2024:")
    print("=" * 60)
    for name in sorted(all_matches):
        seen_in = sorted(set(all_matches[name]))
        marker = "  <-- candidate" if "AK8PFJet" in name and ("SoftDrop" in name or "Mass" in name or "PNet" in name) else ""
        print(f"  {name}  [seen in: {','.join(seen_in)}]{marker}")
    if not all_matches:
        print("  (none found)")


if __name__ == "__main__":
    main()
