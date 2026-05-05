"""
summarize_outputs.py
====================
Walk the trigger-SF pipeline output and produce a single Markdown index of every
figure and table, optionally compiled to PDF.

Inputs (auto-detected under --base-outdir):

  Scale-factor pipeline (scale_factor_2d.py):
    figures_sf_2d/<year>/<prod>/<region>/<mc_sample>/*.png
    output/<year>/<prod>/<region>/<mc_sample>/per_trigger_efficiency_*.csv
    output/<year>/<prod>/<region>/<mc_sample>/per_era_efficiency_*.csv

  VBF parking pipeline (parking_vbf_eff.py):
    figures_parking_vbf/<year>/<method>/*.png
    output/<year>/<method>/*.csv

Outputs:
    <out>           markdown file (default: summary.md under --base-outdir)
    <out>.pdf       (when --pdf, if pandoc / weasyprint / wkhtmltopdf available)

Usage:
    python3 summarize_outputs.py
    python3 summarize_outputs.py --base-outdir my_test_run_both --pdf
    python3 summarize_outputs.py --out report.md --pdf

The Markdown image references are written *relative to the markdown file's
directory* so the document is portable as long as the figures travel with it.
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path


# ---- What to look for ---------------------------------------------------- #

# (filename, title) — emitted in this order if present
SF_FIGURES = [
    ("scale_factor_pt_vs_msd_2d.png",        "2D scale factor (Data eff | MC eff | SF)"),
    ("scale_factor_pt_vs_msd_2d_unc.png",    "SF + statistical uncertainty"),
    ("efficiency_vs_pt.png",                 "Efficiency vs leading-jet pT"),
    ("efficiency_vs_msd.png",                "Efficiency vs leading-jet mSD"),
    ("efficiency_vs_score.png",              "Efficiency vs leading-jet (Xbb + Xcc) score"),
    ("baseline_kinematic_distributions.png", "Baseline kinematic distributions (data vs MC)"),
]

# CSV groups: (prefix or glob, title, predicate to skip)
SF_CSV_RULES = [
    ("per_trigger_efficiency_*.csv",         "Per-trigger efficiency",      lambda n: True),
    ("per_era_efficiency_*.csv",             "Per-era efficiency",          lambda n: "_err" not in n),
    ("per_era_efficiency_*_err.csv",         "Per-era efficiency (σ table)", lambda n: True),
]

VBF_FIGURES = [
    ("vbf_trig_eff_mjj_vs_deta_data.png",    "Data efficiency vs (mjj, |Δη|)"),
    ("vbf_trig_eff_mjj_vs_deta_mc.png",      "MC efficiency vs (mjj, |Δη|)"),
    ("vbf_trig_sf_mjj_vs_deta.png",          "SF map: Data eff | MC eff | SF (mjj, |Δη|)"),
    ("vbf_trig_sf_mjj_vs_deta_unc.png",      "SF map + statistical uncertainty (mjj, |Δη|)"),
    ("vbf_trig_eff_pt1_vs_pt2_data.png",     "Data efficiency vs (pT_j1, pT_j2)"),
    ("vbf_trig_eff_pt1_vs_pt2_mc.png",       "MC efficiency vs (pT_j1, pT_j2)"),
    ("vbf_trig_sf_pt1_vs_pt2.png",           "SF map (pT_j1, pT_j2)"),
    ("vbf_trig_sf_pt1_vs_pt2_unc.png",       "SF map + uncertainty (pT_j1, pT_j2)"),
    # Legacy (mu_tag_parkingvbf method, data-only)
    ("vbf_trig_eff_mjj_vs_deta.png",         "Efficiency vs (mjj, |Δη|)"),
    ("vbf_trig_eff_pt1_vs_pt2.png",          "Efficiency vs (pT_j1, pT_j2)"),
]

VBF_CSV_RULES = [
    ("vbf_parking_trig_eff_*.csv",           "VBF trigger efficiency",   lambda n: True),
    ("vbf_parking_trig_sf_*.csv",            "VBF trigger SF",           lambda n: True),
]


# ---- Helpers ------------------------------------------------------------- #

def csv_to_markdown(path: Path) -> str:
    rows = list(csv.reader(open(path)))
    if not rows:
        return "_(empty CSV)_\n"
    header = rows[0]
    lines = ["| " + " | ".join(header) + " |",
             "|" + "|".join(["---"] * len(header)) + "|"]
    for r in rows[1:]:
        lines.append("| " + " | ".join(r) + " |")
    return "\n".join(lines)


def relpath(p: Path, anchor: Path) -> str:
    """POSIX-style relative path so MD/PDF renderers work on any OS."""
    return os.path.relpath(p, anchor).replace(os.sep, "/")


def find_csvs(directory: Path, rules):
    """Return [(path, title)] for CSVs matching any rule, in rule order."""
    found = []
    seen = set()
    if not directory.is_dir():
        return found
    for pattern, title, predicate in rules:
        for f in sorted(directory.glob(pattern)):
            if f in seen or not predicate(f.name):
                continue
            seen.add(f)
            found.append((f, title))
    return found


# ---- SF section: figures_sf_2d tree ------------------------------------- #

def collect_sf_leaves(figures_root: Path):
    """
    Yield (year, prod, region, mc_sample) for each leaf dir present.
    mc_sample is "" if the run pre-dates the side-by-side layout.
    """
    if not figures_root.is_dir():
        return
    for y in sorted(p for p in figures_root.iterdir() if p.is_dir()):
        for prod in sorted(p for p in y.iterdir() if p.is_dir()):
            for region in sorted(p for p in prod.iterdir() if p.is_dir()):
                subs = [p for p in region.iterdir() if p.is_dir()]
                if subs:
                    for mc in sorted(subs):
                        yield (y.name, prod.name, region.name, mc.name)
                else:
                    yield (y.name, prod.name, region.name, "")


def emit_sf_section(md_lines, base: Path, figures_root: Path, output_root: Path, md_dir: Path):
    leaves = list(collect_sf_leaves(figures_root))
    if not leaves:
        return False
    md_lines.append("\n## Scale-factor pipeline (`scale_factor_2d.py`)\n")

    by_year = defaultdict(list)
    for leaf in leaves:
        by_year[leaf[0]].append(leaf)

    for year in sorted(by_year):
        md_lines.append(f"\n### Year {year}\n")
        for (_, prod, region, mc) in sorted(by_year[year]):
            label = f"prod = `{prod}` · region = `{region}`"
            if mc:
                label += f" · MC = `{mc}`"
            md_lines.append(f"\n#### {label}\n")

            csv_dir = output_root / year / prod / region
            if mc:
                csv_dir = csv_dir / mc
            for csv_path, title in find_csvs(csv_dir, SF_CSV_RULES):
                md_lines.append(f"\n**{title}** — `{csv_path.name}`\n")
                md_lines.append(csv_to_markdown(csv_path))
                md_lines.append("")

            fig_dir = figures_root / year / prod / region
            if mc:
                fig_dir = fig_dir / mc
            for fname, title in SF_FIGURES:
                f = fig_dir / fname
                if f.exists():
                    md_lines.append(f"\n**{title}**\n")
                    md_lines.append(f"![{title}]({relpath(f, md_dir)})")
                    md_lines.append("")
    return True


# ---- VBF parking section ------------------------------------------------ #

def collect_vbf_leaves(figures_root: Path):
    """Yield (year, method) for each leaf under figures_parking_vbf/."""
    if not figures_root.is_dir():
        return
    for y in sorted(p for p in figures_root.iterdir() if p.is_dir()):
        subs = [p for p in y.iterdir() if p.is_dir()]
        if subs:
            for method in sorted(subs):
                yield (y.name, method.name)
        else:
            # legacy layout: figures directly under year/
            yield (y.name, "")


def emit_vbf_section(md_lines, base: Path, figures_root: Path, output_root: Path, md_dir: Path):
    leaves = list(collect_vbf_leaves(figures_root))
    if not leaves:
        return False
    md_lines.append("\n## VBF parking pipeline (`parking_vbf_eff.py`)\n")
    for year, method in leaves:
        label = f"Year {year}"
        if method:
            label += f" · method = `{method}`"
        md_lines.append(f"\n### {label}\n")

        csv_dir = output_root / year / method if method else output_root / year
        for csv_path, title in find_csvs(csv_dir, VBF_CSV_RULES):
            md_lines.append(f"\n**{title}** — `{csv_path.name}`\n")
            md_lines.append(csv_to_markdown(csv_path))
            md_lines.append("")

        fig_dir = figures_root / year / method if method else figures_root / year
        for fname, title in VBF_FIGURES:
            f = fig_dir / fname
            if f.exists():
                md_lines.append(f"\n**{title}**\n")
                md_lines.append(f"![{title}]({relpath(f, md_dir)})")
                md_lines.append("")
    return True


# ---- PDF compilation ---------------------------------------------------- #

HTML_TEMPLATE = """<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>{title}</title>
<style>
  body {{ font-family: -apple-system, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
         max-width: 1200px; margin: 2em auto; padding: 0 2em; color: #222; line-height: 1.5; }}
  h1 {{ border-bottom: 2px solid #333; padding-bottom: 6px; }}
  h2 {{ border-bottom: 1px solid #999; padding-bottom: 3px; margin-top: 2em; }}
  h3 {{ margin-top: 1.5em; }}
  h4 {{ margin-top: 1em; color: #444; }}
  table {{ border-collapse: collapse; margin: 0.5em 0 1em; font-size: 0.85em; }}
  th, td {{ border: 1px solid #bbb; padding: 4px 8px; text-align: right; }}
  th {{ background: #f0f0f0; }}
  td:first-child, th:first-child {{ text-align: left; }}
  img {{ max-width: 100%; height: auto; margin: 0.5em 0; box-shadow: 0 1px 3px rgba(0,0,0,0.2); }}
  code {{ background: #f4f4f4; padding: 1px 4px; border-radius: 3px; font-size: 0.9em; }}
  em {{ color: #666; }}
  @media print {{ body {{ max-width: none; margin: 0; padding: 1cm; }} h2 {{ page-break-before: always; }} }}
</style>
</head>
<body>
{body}
</body>
</html>
"""


def try_html(md_path: Path):
    """
    Render the markdown to a self-contained HTML file (pure Python via the
    `markdown` package). Open in any browser and 'Print → Save as PDF' for a PDF.
    """
    try:
        import markdown as _md
    except ImportError:
        sys.stderr.write(
            "  [html] python-markdown not installed; install with:\n"
            "         pip install --user markdown\n"
        )
        return None
    text = md_path.read_text()
    body = _md.markdown(text, extensions=['tables', 'fenced_code'])
    html_path = md_path.with_suffix(".html")
    html_path.write_text(HTML_TEMPLATE.format(title=md_path.stem, body=body))
    return html_path


def try_pdf(md_path: Path):
    """
    Try a few backends; return (pdf_path, cmd_used) or (None, None).
    Order is from 'no extra deps' to 'heavy LaTeX install'.
    """
    pdf = md_path.with_suffix(".pdf")
    candidates = [
        ["pandoc", str(md_path), "--pdf-engine=weasyprint", "-o", str(pdf)],
        ["pandoc", str(md_path), "--pdf-engine=wkhtmltopdf", "-o", str(pdf)],
        ["pandoc", str(md_path), "--pdf-engine=xelatex",
         "-V", "geometry:margin=1in", "-V", "graphics=true", "-o", str(pdf)],
        ["pandoc", str(md_path), "-o", str(pdf)],
    ]
    for cmd in candidates:
        if shutil.which(cmd[0]) is None:
            continue
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
            if r.returncode == 0 and pdf.exists() and pdf.stat().st_size > 0:
                return pdf, cmd
            sys.stderr.write(f"  [pdf] {' '.join(cmd[:3])} ... failed:\n{(r.stderr or '').strip()[:400]}\n")
        except subprocess.TimeoutExpired:
            sys.stderr.write(f"  [pdf] {' '.join(cmd[:3])} timed out.\n")
    return None, None


# ---- Main --------------------------------------------------------------- #

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base-outdir", default=".",
                    help="Directory containing figures_sf_2d/ and output/ (default: cwd).")
    ap.add_argument("--out", default=None,
                    help="Markdown output path (default: <base-outdir>/summary.md).")
    ap.add_argument("--pdf", action="store_true",
                    help="Also compile a PDF next to the markdown file.")
    args = ap.parse_args()

    base = Path(args.base_outdir).resolve()
    md_path = Path(args.out).resolve() if args.out else base / "summary.md"
    md_path.parent.mkdir(parents=True, exist_ok=True)
    md_dir = md_path.parent

    figs_sf  = base / "figures_sf_2d"
    figs_vbf = base / "figures_parking_vbf"
    out_dir  = base / "output"

    md = []
    md.append("# Trigger-SF Pipeline — Output Index\n")
    md.append(f"_Generated {datetime.now():%Y-%m-%d %H:%M:%S}_\n")
    md.append(f"_Source: `{base}`_\n")
    md.append(
        "\nThis document indexes every figure and CSV produced by the trigger-SF "
        "pipelines. See the project [README.md](../README.md) for the analysis-note "
        "mapping and the methodology behind each plot."
    )

    sf_present  = emit_sf_section(md, base, figs_sf, out_dir, md_dir)
    vbf_present = emit_vbf_section(md, base, figs_vbf, out_dir, md_dir)
    if not (sf_present or vbf_present):
        md.append("\n_No outputs found under this directory._")

    md_path.write_text("\n".join(md))
    print(f"Wrote {md_path}")

    if args.pdf:
        pdf, used_cmd = try_pdf(md_path)
        if pdf:
            print(f"Wrote {pdf}  (via {' '.join(used_cmd[:3])})")
        else:
            print("PDF compilation skipped: no pandoc backend found.")
            print("Falling back to HTML; open it in a browser and 'Print → Save as PDF'.")
            html = try_html(md_path)
            if html:
                print(f"Wrote {html}  (open in browser, Cmd/Ctrl-P → Save as PDF)")
            else:
                print("HTML fallback also unavailable. Install pandoc + xelatex for "
                      "native PDF, or `pip install --user markdown` for HTML.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
