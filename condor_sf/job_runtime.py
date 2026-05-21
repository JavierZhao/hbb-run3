"""Parse condor .log files and report wall + run times per cluster.

Condor .log files record event lines like:
  000 (12345.000.000) 2026-05-20 15:03:47 Job submitted from host: ...
  001 (12345.000.000) 2026-05-20 15:03:50 Job executing on host: ...
  006 (12345.000.000) 2026-05-20 15:04:00 Image size of job updated: ...
  005 (12345.000.000) 2026-05-21 13:17:02 Job terminated.
  009 (12345.000.000) 2026-05-21 14:00:00 Job was aborted by the user
  012 (12345.000.000) 2026-05-21 14:00:00 Job was held.

This script picks the first 000 (submit), first 001 (start), and the last
terminal event (005/009/012) for each cluster, computes wall and run times,
and prints a clean table.

Wall time = terminate - submit  (includes idle queue time)
Run time  = terminate - execute (actual time on a worker)

Usage:
    # All known job directories under condor_sf/:
    python3 job_runtime.py

    # Specific .log files:
    python3 job_runtime.py sf_2d_2022/59052944_0.log extra_trig_2022/*.log

    # Glob across all years:
    python3 job_runtime.py 'sf_2d_*/{*.log,*.log.archived}'
"""
from __future__ import annotations
import argparse
import datetime as dt
import glob
import os
import re
import sys

LOG_LINE_RE = re.compile(
    r"^(?P<code>\d{3}) \((?P<cluster>\d+)\.\d+\.\d+\)\s+(?P<ts>\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})"
)

# 005=terminated, 009=aborted, 012=held (terminal-ish for our purposes too)
TERMINAL_CODES = {"005", "009"}
HOLD_CODE = "012"


def parse_log(path: str):
    """Return dict cluster -> (submit_dt, exec_dt, term_dt, term_code, held_dt or None)."""
    events = {}
    try:
        with open(path) as f:
            for line in f:
                m = LOG_LINE_RE.match(line)
                if not m:
                    continue
                code = m.group("code")
                cl = m.group("cluster")
                ts = dt.datetime.strptime(m.group("ts"), "%Y-%m-%d %H:%M:%S")
                ev = events.setdefault(cl, {"submit": None, "exec": None,
                                            "term": None, "term_code": None,
                                            "held": None})
                if code == "000" and ev["submit"] is None:
                    ev["submit"] = ts
                elif code == "001" and ev["exec"] is None:
                    ev["exec"] = ts
                elif code in TERMINAL_CODES:
                    ev["term"] = ts
                    ev["term_code"] = code
                elif code == HOLD_CODE and ev["held"] is None:
                    ev["held"] = ts
    except (OSError, UnicodeDecodeError) as exc:
        print(f"  WARN: cannot read {path}: {exc}", file=sys.stderr)
    return events


def fmt_delta(d: dt.timedelta) -> str:
    s = int(d.total_seconds())
    h, rem = divmod(s, 3600)
    m, s = divmod(rem, 60)
    if h:
        return f"{h}h{m:02d}m"
    if m:
        return f"{m}m{s:02d}s"
    return f"{s}s"


def _merge_events(dst, src, log_path):
    """Merge src events into dst, keeping earliest submit/exec and latest term.

    Tracks the .log file(s) the cluster's events came from.
    """
    for cl, ev in src.items():
        d = dst.setdefault(cl, {"submit": None, "exec": None,
                                "term": None, "term_code": None,
                                "held": None, "logs": []})
        d["logs"].append(log_path)
        if ev["submit"] and (d["submit"] is None or ev["submit"] < d["submit"]):
            d["submit"] = ev["submit"]
        if ev["exec"] and (d["exec"] is None or ev["exec"] < d["exec"]):
            d["exec"] = ev["exec"]
        # Later terminal event wins (e.g. an abort after a retry)
        if ev["term"] and (d["term"] is None or ev["term"] > d["term"]):
            d["term"] = ev["term"]
            d["term_code"] = ev["term_code"]
        if ev["held"] and (d["held"] is None or ev["held"] < d["held"]):
            d["held"] = ev["held"]


def report(paths, since_iso=None):
    since_dt = None
    if since_iso:
        since_dt = dt.datetime.strptime(since_iso, "%Y-%m-%d %H:%M:%S")

    # Collect events from every log first, deduped by cluster
    by_cluster = {}
    for p in paths:
        if not os.path.isfile(p):
            continue
        _merge_events(by_cluster, parse_log(p), os.path.relpath(p))

    rows = []
    for cl, ev in by_cluster.items():
        sub, exe, term, code, held = (ev["submit"], ev["exec"], ev["term"],
                                      ev["term_code"], ev["held"])
        if since_dt and sub and sub < since_dt:
            continue
        status = (
            "TERM"      if code == "005" else
            "ABORT"     if code == "009" else
            "HELD"      if held and not term else
            "RUN/IDLE"  if exe else
            "PENDING"
        )
        wall = (term - sub) if (term and sub) else None
        run  = (term - exe) if (term and exe) else None
        # Pick a representative log path (shortest = least-suffixed = canonical)
        log_repr = min(ev["logs"], key=len) if ev["logs"] else "-"
        rows.append({
            "cluster": cl,
            "log": log_repr,
            "status": status,
            "submit": sub.strftime("%m-%d %H:%M") if sub else "-",
            "start":  exe.strftime("%m-%d %H:%M") if exe else "-",
            "end":    term.strftime("%m-%d %H:%M") if term else "-",
            "wall":   fmt_delta(wall) if wall else "-",
            "run":    fmt_delta(run)  if run  else "-",
        })

    if not rows:
        print("(no clusters parsed)")
        return

    rows.sort(key=lambda r: (r["submit"], r["cluster"]))
    hdr = f"{'cluster':>10s}  {'status':<8s}  {'submit':<11s}  {'start':<11s}  {'end':<11s}  {'wall':>8s}  {'run':>8s}  log"
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        print(f"{r['cluster']:>10s}  {r['status']:<8s}  {r['submit']:<11s}  {r['start']:<11s}  {r['end']:<11s}  {r['wall']:>8s}  {r['run']:>8s}  {r['log']}")


def discover_logs(roots):
    """Find .log files under given roots (default current dir)."""
    out = []
    for root in roots:
        if os.path.isfile(root):
            out.append(root)
            continue
        for sub in ("sf_2d_*", "extra_trig_*", "scan_*"):
            for d in glob.glob(os.path.join(root, sub)):
                for ext in (".log", ".log.archived", ".log.stuck",
                            ".log.killed_pre_round_fix",
                            ".log.biased_pre_v3plus_fix",
                            ".log.data_xrootd_failed",
                            ".log.failed_no_max"):
                    out.extend(glob.glob(os.path.join(d, f"*{ext}")))
    return sorted(set(out))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("paths", nargs="*", default=[],
                   help="Explicit .log files. If none, auto-discovers under cwd "
                        "(sf_2d_*/*.log, extra_trig_*/*.log, scan_*/*.log + archived).")
    p.add_argument("--since", default=None,
                   help='Only show clusters submitted at/after this time '
                        '(format "YYYY-MM-DD HH:MM:SS").')
    args = p.parse_args()

    paths = args.paths if args.paths else discover_logs(["."])
    report(paths, since_iso=args.since)


if __name__ == "__main__":
    main()
