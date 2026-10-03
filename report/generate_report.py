#!/usr/bin/env python3
"""Build a compact BoltzGen-style Level 2 report.

The report stage deliberately reads the canonical ``metrics/level2.csv`` and
``summary.tsv`` outputs.  It never parses model-private confidence files, so a
report can be regenerated after raw model outputs have been archived.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable


METRICS = ("plddt", "ptm", "iptm", "ranking_score", "confidence_score", "complex_plddt", "pae", "ipsae", "apae")
LABELS = {
    "plddt": "pLDDT",
    "ptm": "pTM",
    "iptm": "ipTM",
    "ranking_score": "Ranking score",
    "confidence_score": "Confidence",
    "complex_plddt": "Complex pLDDT",
    "pae": "Mean PAE",
    "ipsae": "ipSAE",
    "apae": "Interface mean PAE (aPAE)",
}


def _num(value: object) -> float | None:
    try:
        if value is None or str(value).strip() in {"", "NA", "N/A", "nan", "None"}:
            return None
        x = float(value)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def _model_name(value: str) -> str:
    return re.sub(r"_(?:apo|holo)$", "", value or "unknown")


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _read_metrics(run_dir: Path) -> list[dict[str, str]]:
    primary = run_dir / "metrics" / "level2.csv"
    paths = [primary] if primary.is_file() else sorted((run_dir / "metrics").glob("level2*.csv"))
    if not paths:
        raise FileNotFoundError(f"missing {run_dir / 'metrics' / 'level2.csv'}")
    rows: list[dict[str, str]] = []
    for path in paths:
        with path.open(encoding="utf-8", newline="") as handle:
            rows.extend(csv.DictReader(handle))
    return rows


def _score(row: dict[str, str]) -> float:
    """Return a comparable score without treating monomer ipTM=0 as failure."""
    values = []
    for key in ("plddt", "ptm", "iptm"):
        value = _num(row.get(key))
        if value is None:
            continue
        if key == "plddt":
            value /= 100.0
        # ipTM is structurally undefined for a monomer and is emitted as zero.
        if key == "iptm" and value <= 0:
            chain_text = str(row.get("chain_iptm", ""))
            chain_values = [float(x) for x in re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", chain_text)]
            if not any(abs(x) > 1e-6 for x in chain_values):
                continue
        values.append(value)
    return sum(values) / len(values) if values else float("nan")


def _prepare_rows(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    out = []
    for row in rows:
        if str(row.get("status", "")).lower() not in {"success", "ok", ""}:
            continue
        item: dict[str, object] = dict(row)
        item["design_id"] = row.get("design_id", "unknown")
        item["model"] = _model_name(row.get("model", "unknown"))
        item["state"] = "holo" if str(row.get("model", "")).endswith("_holo") else (
            "apo" if str(row.get("model", "")).endswith("_apo") else ""
        )
        item["score"] = _score(row)
        for key in METRICS:
            aliases = {
                "pae": ("pae", "complex_pde", "overall_mean_pae", "mean_pae"),
                "ipsae": ("ipsae", "overall_ipsae", "complex_ipsae"),
                "apae": ("apae", "complex_ipde", "mean_pae_at_interface", "interface_mean_pae"),
            }.get(key, (key,))
            item[key] = next((_num(row.get(alias)) for alias in aliases if _num(row.get(alias)) is not None), None)
        out.append(item)
    return out


def _write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: "" if row.get(key) is None else row.get(key, "") for key in fields})


def _summary_stats(rows: list[dict[str, object]], outdir: Path) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Write compact aggregate tables; never expand one PDF row per design."""
    import numpy as np
    summary = []
    for key in METRICS:
        vals = [float(r[key]) for r in rows if r[key] is not None]
        if not vals:
            continue
        q = np.percentile(vals, [5, 25, 50, 75, 95])
        summary.append({"metric": LABELS[key], "key": key, "count": len(vals),
                        "mean": float(np.mean(vals)), "median": float(q[2]),
                        "std": float(np.std(vals)), "min": float(np.min(vals)), "max": float(np.max(vals)),
                        "p05": float(q[0]), "p25": float(q[1]), "p75": float(q[3]), "p95": float(q[4])})
    _write_tsv(outdir / "summary_stats.tsv", summary,
               ["metric", "key", "count", "mean", "median", "std", "min", "max", "p05", "p25", "p75", "p95"])
    model_summary = []
    for model in sorted({str(r["model"]) for r in rows}):
        subset = [r for r in rows if str(r["model"]) == model]
        rec = {"model": model, "records": len(subset), "designs": len({str(r["design_id"]) for r in subset})}
        for key in METRICS:
            vals = [float(r[key]) for r in subset if r[key] is not None]
            rec[f"{key}_mean"] = float(np.mean(vals)) if vals else ""
            rec[f"{key}_median"] = float(np.median(vals)) if vals else ""
        model_summary.append(rec)
    fields = ["model", "records", "designs"] + [f"{k}_{s}" for k in METRICS for s in ("mean", "median")]
    _write_tsv(outdir / "model_stats.tsv", model_summary, fields)
    return summary, model_summary


def _plot_stats(rows: list[dict[str, object]], outdir: Path) -> list[Path]:
    # Import lazily so prepare/run do not require plotting dependencies.
    import matplotlib
    # The controller image on f101 has a broken Agg/Pillow combination while
    # its SVG backend is healthy.  Render SVG first and rasterize with the
    # system ImageMagick so reports work in both controller and GPU envs.
    matplotlib.use("svg")
    import matplotlib.pyplot as plt
    import numpy as np

    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    models = sorted({str(row["model"]) for row in rows})
    colors = {model: plt.cm.tab10(i % 10) for i, model in enumerate(models)}
    paths: list[Path] = []

    def save_figure(fig, path: Path) -> None:
        svg = path.with_suffix(".svg")
        fig.savefig(svg, bbox_inches="tight")
        try:
            subprocess.run(["convert", str(svg), str(path)], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except (OSError, subprocess.CalledProcessError):
            # Keep the vector source as a usable fallback if ImageMagick is
            # unavailable on another installation.
            path.write_bytes(svg.read_bytes())

    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5), constrained_layout=True)
    for ax, key in zip(axes.flat, METRICS):
        for model in models:
            values = [row[key] for row in rows if row["model"] == model and row[key] is not None]
            if values:
                bins = min(12, max(4, len(values)))
                ax.hist(values, bins=bins, alpha=.52, label=model, color=colors[model], edgecolor="white")
        ax.set_title(f"{LABELS[key]} distribution")
        ax.set_xlabel(LABELS[key])
        ax.set_ylabel("Count")
        ax.grid(alpha=.18)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=max(1, len(labels)), frameon=False)
    path = outdir / "metric_distributions.png"
    save_figure(fig, path)
    plt.close(fig)
    paths.append(path)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    for model in models:
        subset = [row for row in rows if row["model"] == model]
        for ax, xkey, ykey in ((axes[0], "plddt", "ptm"), (axes[1], "ptm", "iptm")):
            xy = [(row[xkey], row[ykey]) for row in subset if row[xkey] is not None and row[ykey] is not None]
            if xy:
                x, y = zip(*xy)
                axes[0 if xkey == "plddt" else 1].scatter(x, y, s=48, alpha=.8, label=model, color=colors[model])
    axes[0].set(xlabel="pLDDT", ylabel="pTM", title="Confidence relationship")
    axes[1].set(xlabel="pTM", ylabel="ipTM", title="Interface relationship")
    for ax in axes:
        ax.grid(alpha=.18)
    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=max(1, len(labels)), frameon=False)
    path = outdir / "metric_scatter.png"
    save_figure(fig, path)
    plt.close(fig)
    paths.append(path)

    # A design×model heatmap is intentionally bounded for high-throughput runs.
    # Use the best 40 designs by composite score, rather than allocating a
    # canvas proportional to thousands of input sequences.
    ranked_designs = sorted(
        ((str(d), max((float(r["score"]) for r in rows if str(r["design_id"]) == str(d) and not math.isnan(float(r["score"]))), default=-1.0))
         for d in {str(row["design_id"]) for row in rows}), key=lambda x: (-x[1], x[0])
    )
    designs = [d for d, _ in ranked_designs[:40]]
    matrix = np.full((len(designs), len(models)), np.nan)
    for i, design in enumerate(designs):
        for j, model in enumerate(models):
            values = [row["score"] for row in rows if row["design_id"] == design and row["model"] == model and not math.isnan(float(row["score"]))]
            if values:
                matrix[i, j] = float(np.mean(values))
    from matplotlib.patches import Rectangle
    fig, ax = plt.subplots(figsize=(max(5, len(models) * 1.3), max(3.5, len(designs) * .55)), constrained_layout=True)
    # Rectangle cells avoid the QuadMesh path, which is broken in the f101
    # controller's matplotlib Agg/SVG stack for object-backed arrays.
    cmap = plt.get_cmap("viridis")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix[i, j]
            face = "#E5E7EB" if np.isnan(value) else cmap(float(value))
            ax.add_patch(Rectangle((j - .5, i - .5), 1, 1, facecolor=face, edgecolor="white", linewidth=.8))
    ax.set(xticks=range(len(models)), xticklabels=models, yticks=range(len(designs)), yticklabels=designs, title="Design × model composite score")
    ax.set_xlim(-.5, len(models) - .5)
    ax.set_ylim(len(designs) - .5, -.5)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if not np.isnan(matrix[i, j]):
                ax.text(j, i, f"{matrix[i, j]:.2f}", ha="center", va="center", color="white" if matrix[i, j] < .6 else "black", fontsize=8)
    ax.text(1.0, -0.12, "cell color: composite score 0–1; gray: unavailable", transform=ax.transAxes, ha="right", va="top", fontsize=8, color="#4B5563")
    path = outdir / "model_heatmap.png"
    save_figure(fig, path)
    plt.close(fig)
    paths.append(path)

    fig, ax = plt.subplots(figsize=(11, 5.5), constrained_layout=True)
    metric_keys = ("plddt", "ptm", "iptm")
    x = np.arange(len(models))
    width = .24
    for offset, key in enumerate(metric_keys):
        means = [np.mean([row[key] for row in rows if row["model"] == model and row[key] is not None]) if any(row["model"] == model and row[key] is not None for row in rows) else np.nan for model in models]
        if key == "plddt":
            means = [v / 100 if not np.isnan(v) else v for v in means]
        ax.bar(x + (offset - 1) * width, means, width, label=LABELS[key])
    ax.set(xticks=x, xticklabels=models, ylabel="Normalized value", title="Mean model metrics")
    ax.set_ylim(0, 1.05)
    ax.grid(axis="y", alpha=.18)
    ax.legend(frameon=False, ncol=3)
    path = outdir / "model_bars.png"
    save_figure(fig, path)
    plt.close(fig)
    paths.append(path)
    return paths


def _rank_designs(rows: list[dict[str, object]], outdir: Path, top_n: int) -> list[dict[str, object]]:
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[str(row["design_id"])].append(row)
    ranked = []
    for design, values in grouped.items():
        valid = [float(row["score"]) for row in values if not math.isnan(float(row["score"]))]
        if not valid:
            continue
        best = max(values, key=lambda row: -math.inf if math.isnan(float(row["score"])) else float(row["score"]))
        ranked.append({
            "design_id": design,
            "score": sum(valid) / len(valid),
            "models_success": len(valid),
            "best_model": best["model"],
            "structure_path": best.get("structure_path", ""),
        })
    ranked.sort(key=lambda row: (-float(row["score"]), str(row["design_id"])))
    for i, row in enumerate(ranked, 1):
        row["rank"] = i
        row["selected"] = "yes" if i <= top_n else "no"
    _write_tsv(outdir / "ranking.tsv", ranked, ["rank", "design_id", "score", "models_success", "best_model", "structure_path", "selected"])
    return ranked[:top_n]


def _render_top(top: list[dict[str, object]], outdir: Path, renderer: str, resolution: str, aa_samples: int) -> list[tuple[dict[str, object], list[Path]]]:
    """Render each design and retain its ranking row beside the image."""
    rendered: list[tuple[dict[str, object], list[Path]]] = []
    for row in top:
        source = Path(str(row.get("structure_path", "")))
        if not source.is_file():
            continue
        target = outdir / "structures" / f"{int(row['rank']):03d}_{row['design_id']}"
        target.mkdir(parents=True, exist_ok=True)
        def metric_text(key: str, label: str) -> str:
            value = row.get(key)
            return f"{label}={float(value):.2f}" if value is not None else f"{label}=NA"
        subtitle = " · ".join(metric_text(k, label) for k, label in (("plddt", "pLDDT"), ("ptm", "pTM"), ("pae", "mean PAE"), ("iptm", "ipTM"), ("ipsae", "ipSAE"), ("apae", "aPAE")))
        command = [sys.executable, renderer, str(source), "--mode", "auto", "--outdir", str(target), "--views", "auto", "--panel", "--panel-cols", "2", "--image-geometry", "1000x800+14+14", "--preset", "paper", "--label", str(row["design_id"]), "--subtitle", subtitle, "--protein-opacity", "0.30", "--pocket-style", "transparent", "--ligand-rep", "stick", "--render-resolution", resolution, "--aa-samples", str(aa_samples)]
        try:
            subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        except (OSError, subprocess.CalledProcessError) as exc:
            (target / "render_error.txt").write_text(str(exc), encoding="utf-8")
            continue
        panel = sorted(target.glob("*_panel.png"))
        # Use the four source figures, not the already-compressed montage.
        # This gives the PDF one design per row: front, side, top, interface.
        views = [p for p in sorted(target.glob("*.png")) if p not in panel]
        if len(views) >= 4:
            rendered.append((row, views[:4]))
            # The renderer leaves VMD/Tachyon intermediates beside the panel;
            # keep the user-facing PNGs and discard only reproducible scratch.
            for scratch in ("*.tga", "*.dat", "*.tcl"):
                for path in target.glob(scratch):
                    path.unlink(missing_ok=True)
            shutil.rmtree(target / "_render_inputs", ignore_errors=True)
    return rendered


def _build_pdf(outdir: Path, run_dir: Path, rows: list[dict[str, object]], ranked: list[dict[str, object]], plots: list[Path], structures: list[tuple[dict[str, object], list[Path]]], pdf_name: str, summary: list[dict[str, object]], model_summary: list[dict[str, object]]) -> Path:
    from reportlab.lib import colors
    from reportlab.lib.enums import TA_CENTER
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import inch
    from reportlab.platypus import Image, PageBreak, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle

    path = outdir / pdf_name
    doc = SimpleDocTemplate(str(path), pagesize=A4, rightMargin=36, leftMargin=36, topMargin=32, bottomMargin=32)
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="ReportTitle", parent=styles["Title"], alignment=TA_CENTER, fontSize=22, leading=26, textColor=colors.HexColor("#17324D")))
    styles.add(ParagraphStyle(name="Small", parent=styles["BodyText"], fontSize=8, leading=10, textColor=colors.HexColor("#4B5563")))
    design_count = len(set(str(row["design_id"]) for row in rows))
    story = [Paragraph("ProFold Level 2 Design Report", styles["ReportTitle"]), Spacer(1, 8), Paragraph(f"Run: {run_dir}<br/>{design_count:,} designs · {len(rows):,} successful model records · aggregate report", styles["Small"]), Spacer(1, 14)]
    model_counts = defaultdict(int)
    for row in rows:
        model_counts[str(row["model"])] += 1
    table_data = [["Model", "Records", "Designs", "Mean pLDDT", "Mean pTM", "Mean ipTM"]]
    for rec in model_summary:
        def fmt(key: str) -> str:
            value = rec.get(f"{key}_mean", "")
            return f"{float(value):.3f}" if value != "" else "—"
        table_data.append([rec["model"], str(rec["records"]), str(rec["designs"]), fmt("plddt"), fmt("ptm"), fmt("iptm")])
    table = Table(table_data, repeatRows=1, colWidths=[1.5*inch, 1.25*inch, 1.1*inch, 1.1*inch, 1.1*inch])
    table.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#17324D")), ("TEXTCOLOR", (0, 0), (-1, 0), colors.white), ("GRID", (0, 0), (-1, -1), .25, colors.HexColor("#CBD5E1")), ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F1F5F9")]), ("ALIGN", (1, 1), (-1, -1), "CENTER")]))
    story += [table, Spacer(1, 14), Paragraph("Metric summary (all successful records)", styles["Heading2"])]
    stat_data = [["Metric", "N", "Mean", "Median", "SD", "P05", "P95", "Min", "Max"]]
    for rec in summary:
        stat_data.append([rec["metric"], str(rec["count"]), *(f"{float(rec[k]):.3f}" for k in ("mean", "median", "std", "p05", "p95", "min", "max"))])
    stat_table = Table(stat_data, repeatRows=1, colWidths=[1.25*inch, .45*inch, .65*inch, .65*inch, .55*inch, .55*inch, .55*inch, .6*inch, .6*inch])
    stat_table.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#17324D")), ("TEXTCOLOR", (0, 0), (-1, 0), colors.white), ("GRID", (0, 0), (-1, -1), .25, colors.HexColor("#CBD5E1")), ("FONTSIZE", (0, 0), (-1, -1), 7), ("ALIGN", (1, 1), (-1, -1), "CENTER")]))
    story += [stat_table, Spacer(1, 14), Paragraph("Top designs (bounded to requested Top-N)", styles["Heading2"])]
    rank_data = [["Rank", "Design", "Score", "Best model", "Models"]] + [[str(r["rank"]), str(r["design_id"]), f"{float(r['score']):.3f}", str(r["best_model"]), str(r["models_success"])] for r in ranked]
    rank_table = Table(rank_data, repeatRows=1, colWidths=[.55*inch, 2.0*inch, .8*inch, 1.1*inch, .75*inch])
    rank_table.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#2E6F95")), ("TEXTCOLOR", (0, 0), (-1, 0), colors.white), ("GRID", (0, 0), (-1, -1), .25, colors.HexColor("#CBD5E1")), ("ALIGN", (0, 0), (-1, -1), "CENTER")]))
    story += [rank_table, PageBreak(), Paragraph("Metric distributions", styles["Heading1"])]
    for plot in plots:
        story += [Image(str(plot), width=7.0*inch, height=4.05*inch), Spacer(1, 8)]
        if plot.name == "metric_scatter.png" or plot.name == "model_bars.png":
            story.append(PageBreak())
    if structures:
        # 4 columns × 6 rows: each cell keeps the panel's original aspect
        # ratio, so a long/short structure is never stretched or clipped.
        from PIL import Image as PILImage
        story += [PageBreak(), Paragraph("Top-N structure views (4 × 6)", styles["Heading1"])]
        cell_w, cell_h = 1.78 * inch, 1.38 * inch
        cells = []
        for row, images in structures:
            row_cells = []
            for image in images[:4]:
                with PILImage.open(image) as im:
                    iw, ih = im.size
                scale = min((cell_w - 0.08 * inch) / iw, (cell_h - 0.22 * inch) / ih)
                row_cells.append(Image(str(image), width=iw * scale, height=ih * scale))
            while len(row_cells) < 4:
                row_cells.append("")
            cells.append(row_cells)
        for start in range(0, len(cells), 6):
            chunk = cells[start:start + 6]
            while len(chunk) < 6:
                chunk.append([""] * 4)
            grid = chunk
            table = Table(grid, colWidths=[cell_w] * 4, rowHeights=[cell_h] * 6)
            table.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "MIDDLE"), ("ALIGN", (0, 0), (-1, -1), "CENTER"), ("BOX", (0, 0), (-1, -1), .25, colors.HexColor("#CBD5E1")), ("INNERGRID", (0, 0), (-1, -1), .2, colors.HexColor("#E5E7EB"))]))
            story += [table]
            if start + 24 < len(cells):
                story.append(PageBreak())
    doc.build(story)
    return path


def generate_report(run_dir: str | Path, top_n: int = 10, render_structures: bool = True, renderer: str = "/xcfhome/yhliu/14_magpcr/z_zcodex/pipeline_scripts/visualization/render_structure_panel.py", resolution: str = "1200x900", aa_samples: int = 4, pdf_name: str = "profold_level2_report.pdf") -> Path:
    run_path = Path(run_dir).resolve()
    outdir = run_path / "report"
    outdir.mkdir(parents=True, exist_ok=True)
    rows = _prepare_rows(_read_metrics(run_path))
    if not rows:
        raise ValueError(f"no successful metrics in {run_path}")
    _write_tsv(outdir / "metrics_long.tsv", rows, ["design_id", "model", "state", "status", "plddt", "ptm", "iptm", "ranking_score", "confidence_score", "complex_plddt", "score", "structure_path"])
    ranked = _rank_designs(rows, outdir, max(0, int(top_n)))
    summary, model_summary = _summary_stats(rows, outdir)
    plots = _plot_stats(rows, outdir)
    structures = _render_top(ranked, outdir, renderer, resolution, aa_samples) if render_structures and Path(renderer).is_file() else []
    report_path = _build_pdf(outdir, run_path, rows, ranked, plots, structures, pdf_name, summary, model_summary)
    (outdir / "report_manifest.json").write_text(json.dumps({"pdf": str(report_path), "top_n": top_n, "plots": [str(p) for p in plots], "structures": [str(p) for _, images in structures for p in images], "designs": len(set(str(row["design_id"]) for row in rows)), "records": len(rows), "summary_stats": str(outdir / "summary_stats.tsv"), "model_stats": str(outdir / "model_stats.tsv")}, indent=2) + "\n", encoding="utf-8")
    return report_path


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate ProFold Level 2 plots, Top-N structure panels, and PDF report")
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--top-n", type=int, default=10)
    parser.add_argument("--no-structures", action="store_true")
    parser.add_argument("--renderer", default="/xcfhome/yhliu/14_magpcr/z_zcodex/pipeline_scripts/visualization/render_structure_panel.py")
    parser.add_argument("--render-resolution", default="1200x900")
    parser.add_argument("--aa-samples", type=int, default=4)
    parser.add_argument("--pdf-name", default="profold_level2_report.pdf")
    args = parser.parse_args()
    print(generate_report(args.run_dir, args.top_n, not args.no_structures, args.renderer, args.render_resolution, args.aa_samples, args.pdf_name))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
