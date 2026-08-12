from __future__ import annotations

import math
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.font_manager import FontProperties
from matplotlib.patches import PathPatch
from matplotlib.textpath import TextPath
from matplotlib.transforms import Affine2D

COLORS = {"A": "#009E73", "C": "#0072B2", "G": "#E69F00", "T": "#D55E00"}


def logo(ax, matrix: np.ndarray, title: str) -> None:
    fp = FontProperties(family="DejaVu Sans", weight="bold")
    letters = "ACGT"
    entropy = -(matrix * np.log2(np.clip(matrix, 1e-12, 1))).sum(axis=1)
    heights = matrix * np.clip(2 - entropy, 0, 2)[:, None]
    for x, row in enumerate(heights):
        y = 0.0
        for idx in np.argsort(row):
            h = float(row[idx])
            if h <= 0:
                continue
            char = letters[idx]
            path = TextPath((0, 0), char, size=1, prop=fp)
            box = path.get_extents()
            trans = Affine2D().scale(0.82 / max(box.width, 1e-9), h / max(box.height, 1e-9)).translate(x + 0.09, y)
            ax.add_patch(PathPatch(path, transform=trans + ax.transData, color=COLORS[char], lw=0))
            y += h
    ax.set_xlim(0, matrix.shape[0])
    ax.set_ylim(0, 2.05)
    ax.set_xticks(np.arange(matrix.shape[0]) + 0.5, np.arange(1, matrix.shape[0] + 1))
    ax.set_ylabel("Information (bits)")
    ax.set_xlabel("Position")
    ax.set_title(title, loc="left", fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)


def read_streme(path: Path) -> list[tuple[str, float, np.ndarray]]:
    root = ET.parse(path).getroot()
    out = []
    for motif in root.findall(".//motif"):
        name = motif.attrib.get("id", motif.attrib.get("alt", "motif"))
        evalue = float(motif.attrib.get("test_evalue", motif.attrib.get("evalue", "nan")))
        rows = []
        for pos in motif.findall(".//pos"):
            rows.append([float(pos.attrib.get(k, 0)) for k in ("A", "C", "G", "T")])
        if rows:
            out.append((name, evalue, np.asarray(rows)))
    return out


def read_meme(path: Path, wanted: set[str]) -> dict[str, tuple[str, np.ndarray]]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    out: dict[str, tuple[str, np.ndarray]] = {}
    i = 0
    while i < len(lines):
        if lines[i].startswith("MOTIF "):
            parts = lines[i].split(maxsplit=2)
            motif_id = parts[1]
            alt = parts[2] if len(parts) > 2 else motif_id
            i += 1
            while i < len(lines) and "letter-probability matrix" not in lines[i]:
                i += 1
            i += 1
            rows = []
            while i < len(lines):
                vals = lines[i].split()
                if len(vals) != 4:
                    break
                try:
                    rows.append([float(v) for v in vals])
                except ValueError:
                    break
                i += 1
            if motif_id in wanted and rows:
                out[motif_id] = (alt, np.asarray(rows))
        i += 1
    return out


def save(fig: plt.Figure, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".png"), dpi=400, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main(run_dir: str) -> None:
    run = Path(run_dir).resolve()
    promoter = run / "08_Motif" / "Promoter"
    figures = promoter / "Figures"
    streme = read_streme(promoter / "STREME" / "streme.xml")
    for idx, (name, evalue, matrix) in enumerate(streme[:5], 1):
        fig, ax = plt.subplots(figsize=(max(6.5, matrix.shape[0] * 0.55), 3.5))
        logo(ax, matrix, f"STREME de novo candidate {idx}: {name} (E={evalue:.2g})")
        fig.text(0.01, 0.01, "Exploratory sequence motif; not evidence of direct binding.", fontsize=9, color="#555555")
        save(fig, figures / f"STREME_DeNovo_Motif_{idx}_Logo")

    if streme:
        fig, axes = plt.subplots(len(streme[:5]), 1, figsize=(10, 2.7 * len(streme[:5])), constrained_layout=True)
        axes = np.atleast_1d(axes)
        for idx, ((name, evalue, matrix), ax) in enumerate(zip(streme[:5], axes), 1):
            logo(ax, matrix, f"Candidate {idx}: {name} (E={evalue:.2g})")
        fig.suptitle("STREME de novo promoter motif candidates", fontsize=16, fontweight="bold")
        save(fig, figures / "STREME_Motif_Summary")

    ame = pd.read_csv(promoter / "AME_JASPAR" / "ame.tsv", sep="\t", comment="#")
    ame["adj_p"] = pd.to_numeric(ame["adj_p-value"])
    ame["E"] = pd.to_numeric(ame["E-value"])
    ame["label"] = ame["motif_alt_ID"].fillna(ame["motif_ID"]).astype(str)
    shown = ame.head(15).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8.5, max(4, 0.42 * len(shown) + 1.8)))
    colors = np.where(shown["E"] < 0.05, "#D55E00", "#7A7A7A")
    ax.barh(shown["label"], -np.log10(shown["adj_p"].clip(lower=1e-300)), color=colors, edgecolor="#333333", linewidth=0.4)
    for y, (_, row) in enumerate(shown.iterrows()):
        ax.text(-math.log10(max(row["adj_p"], 1e-300)) + 0.03, y, f"E={row['E']:.2g}", va="center", fontsize=8)
    ax.set_xlabel("−log10 adjusted P value")
    ax.set_ylabel("")
    ax.set_title("JASPAR known-motif enrichment candidates", loc="left", fontweight="bold", pad=34)
    ax.text(0, 1.015, "Orange: E<0.05; grey: exploratory candidate only", transform=ax.transAxes, fontsize=9, color="#555555")
    ax.spines[["top", "right"]].set_visible(False)
    save(fig, figures / "AME_Known_Motif_Enrichment")

    ame_html = """<!doctype html><html lang='en'><head><meta charset='utf-8'>
<title>AME JASPAR motif enrichment</title><style>body{font:15px/1.5 system-ui,sans-serif;max-width:1200px;margin:auto;padding:24px}table{border-collapse:collapse;width:100%;font-size:12px}th,td{border-bottom:1px solid #ddd;padding:5px;text-align:left}img{max-width:900px;width:100%}.note{border-left:4px solid #e69f00;background:#fff7e6;padding:12px}</style></head><body>
<h1>AME JASPAR known-motif enrichment</h1><p class='note'>Exploratory promoter enrichment only. This does not establish direct TF binding or causality.</p>
<p><img src='../Figures/AME_Known_Motif_Enrichment.png' alt='AME known motif enrichment'></p>
""" + ame.to_html(index=False, border=0) + "</body></html>"
    (promoter / "AME_JASPAR" / "ame.html").write_text(ame_html, encoding="utf-8")

    jaspar = Path.cwd() / "resources" / "data" / "shared" / "hg38-motif" / "JASPAR2026_CORE_vertebrates_non-redundant_pfms_meme.txt"
    matrices = read_meme(jaspar, set(ame.head(5)["motif_ID"].astype(str)))
    for _, row in ame.head(5).iterrows():
        motif_id = str(row["motif_ID"])
        if motif_id not in matrices:
            continue
        alt, matrix = matrices[motif_id]
        fig, ax = plt.subplots(figsize=(max(6.5, matrix.shape[0] * 0.55), 3.5))
        logo(ax, matrix, f"JASPAR {motif_id} {alt}: AME candidate (E={row['E']:.2g})")
        fig.text(0.01, 0.01, "Known-motif enrichment is exploratory and does not establish TF binding or causality.", fontsize=9, color="#555555")
        save(fig, figures / f"AME_{clean(alt)}_Logo")

    status = pd.DataFrame({
        "output": ["STREME logos", "STREME summary", "AME enrichment", "AME known-motif logos"],
        "status": ["generated"] * 4,
        "interpretation": ["exploratory candidate"] * 4,
    })
    status.to_csv(promoter / "Tables" / "Motif_Figure_Manifest.tsv", sep="\t", index=False)
    index_html = """<!doctype html><html lang='en'><head><meta charset='utf-8'><title>Promoter motif analysis</title><style>body{font:16px/1.5 system-ui,sans-serif;max-width:1000px;margin:auto;padding:24px}img{max-width:100%;border:1px solid #ddd}.note{border-left:4px solid #e69f00;background:#fff7e6;padding:12px}</style></head><body><h1>Promoter motif analysis</h1><p class='note'>Counts-only promoter motif results are exploratory; they are not evidence of direct binding or causality.</p><ul><li><a href='STREME/streme.html'>STREME interactive report</a></li><li><a href='AME_JASPAR/ame.html'>AME JASPAR report</a></li><li><a href='Tables/Motif_Run_Summary.tsv'>Run summary</a></li></ul><img src='Figures/STREME_Motif_Summary.png' alt='STREME motif summary'><img src='Figures/AME_Known_Motif_Enrichment.png' alt='AME enrichment'></body></html>"""
    (promoter / "index.html").write_text(index_html, encoding="utf-8")


def clean(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "_.-" else "_" for ch in value)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: render_motif_figures.py <run_dir>")
    main(sys.argv[1])
