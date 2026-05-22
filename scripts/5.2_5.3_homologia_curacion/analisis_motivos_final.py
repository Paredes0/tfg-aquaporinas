#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analisis_motivos_final.py — Huella de motivos peptídicos por sub-subfamilia (Figura 5 del TFG).

Parsea la sección MOTIF DIAGRAM de ALL_AQP.txt (MEME sobre las 121 acuaporinas curadas,
15 motivos M1-M15) y la cruza con la asignación de sub-subfamilia de tabla_aqp_ordenada.csv.
Genera el heatmap de frecuencia de motivos por sub-subfamilia e imprime la tabla de
frecuencias por subfamilia.

Este es el análisis de motivos DEFINITIVO del TFG (15 motivos sobre 121 curadas),
distinto del MEME exploratorio sobre las 258 candidatas usado durante el curado.

Inputs (data/curado/):
    ALL_AQP.txt                 MEME final (15 motivos sobre las 121 curadas)
    tabla_aqp_ordenada.csv      gene_id -> aqp_family_subfamily (sub-subfamilia)

Output (results/):
    HEATMAP_Frecuencia_Motivos_Sub-subfamilias.png

Nota tipográfica: la figura NO lleva título embebido; el título va en el pie de figura
del TFG (norma APA/UCAM).
"""
from __future__ import annotations

import csv
import os
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.common import config as _cfg

MEME_PATH = _cfg.CURADO_DIR / "ALL_AQP.txt"
TABLA_PATH = _cfg.CURADO_DIR / "tabla_aqp_ordenada.csv"
OUT_DIR = _cfg.ensure_results()

NUM_MOTIFS = 15
SUBSUB_ORDER = [
    "FaPIP1", "FaPIP2", "FaTIP1", "FaTIP2", "FaTIP3", "FaTIP4", "FaTIP5",
    "FaNIP1", "FaNIP2", "FaNIP4", "FaNIP5", "FaNIP6", "FaNIP7",
    "FaSIP1", "FaSIP2", "FaXIP1", "FaXIP2",
]
SUBFAMS = ["NIP", "TIP", "PIP", "SIP", "XIP"]
_MOTIF_RE = re.compile(r"\[(\d+)\(")


def load_gene_to_subsub(path: Path) -> dict[str, str]:
    mapping = {}
    with open(path, encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            mapping[row["gene_id"]] = row["aqp_family_subfamily"]
    return mapping


def parse_motif_diagram(path: Path) -> dict[str, str]:
    """Devuelve {seq_name: motif_diagram} de la sección MOTIF DIAGRAM de MEME."""
    content = path.read_text(encoding="utf-8")
    idx = content.find("SEQUENCE NAME            COMBINED P-VALUE  MOTIF DIAGRAM")
    diagram = content[idx:].split("\n")[2:]
    seq_motifs: dict[str, str] = {}
    current_seq, current_diagram = None, ""
    for line in diagram:
        if line.startswith("-------------"):
            break
        if not line.strip() or line.startswith("*") or line.startswith("----"):
            if current_seq:
                seq_motifs[current_seq] = current_diagram
                current_seq, current_diagram = None, ""
            continue
        if not line.startswith(" ") and not line.startswith("\t"):
            if current_seq:
                seq_motifs[current_seq] = current_diagram
            parts = line.split()
            current_seq = parts[0]
            current_diagram = parts[2] if len(parts) > 2 else ""
        else:
            current_diagram += line.strip()
        if current_diagram.endswith("\\"):
            current_diagram = current_diagram[:-1]
    if current_seq and current_seq not in seq_motifs:
        seq_motifs[current_seq] = current_diagram
    return seq_motifs


def extract_gene_id(seq_name: str) -> str:
    """De 'PIP1;1-Fxa7Dg02199' extrae 'Fxa7Dg02199'."""
    if "-" in seq_name:
        parts = seq_name.split("-", 1)
        if len(parts) == 2:
            return parts[1].split("-")[0]
    return seq_name


def main() -> None:
    gene_to_subsub = load_gene_to_subsub(TABLA_PATH)
    print(f"gene_id -> sub-subfamilia mapeados: {len(gene_to_subsub)}")

    seq_motifs = parse_motif_diagram(MEME_PATH)

    # Agrupar conjuntos de motivos por sub-subfamilia
    subsub_motifs: dict[str, list[set[int]]] = {}
    for seq, diag in seq_motifs.items():
        gid = extract_gene_id(seq)
        subsub = gene_to_subsub.get(gid)
        if subsub:
            subsub_motifs.setdefault(subsub, []).append(
                {int(m) for m in _MOTIF_RE.findall(diag)}
            )
    print(f"Secuencias matcheadas: {sum(len(v) for v in subsub_motifs.values())}")

    present = [sf for sf in SUBSUB_ORDER if sf in subsub_motifs]
    data = {}
    for sf in present:
        seqs = subsub_motifs[sf]
        n = len(seqs)
        data[sf] = {f"M{m}": (sum(1 for s in seqs if m in s) / n * 100 if n else 0)
                    for m in range(1, NUM_MOTIFS + 1)}
    df = pd.DataFrame(data).T
    df.index.name = "Sub-subfamilia"
    df.index = [sf.replace("Fa", "") for sf in df.index]

    # Heatmap (sin título embebido: el título va en el pie de figura del TFG)
    plt.figure(figsize=(14, 8))
    sns.heatmap(df, cmap="YlGnBu", annot=True, fmt=".0f", linewidths=0.5,
                cbar_kws={"label": "% de secuencias con el motivo"}, vmin=0, vmax=100)
    plt.xlabel("Motivo MEME")
    plt.ylabel("Sub-subfamilia")
    plt.tight_layout()
    out_path = OUT_DIR / "HEATMAP_Frecuencia_Motivos_Sub-subfamilias.png"
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Heatmap guardado: {out_path}")

    # Tabla de frecuencias por subfamilia (a stdout)
    print("\nFrecuencia de cada motivo por subfamilia (%):")
    sub_seqs = {s: [] for s in SUBFAMS}
    for seq, diag in seq_motifs.items():
        for sf in SUBFAMS:
            if seq.startswith(sf):
                sub_seqs[sf].append({int(m) for m in _MOTIF_RE.findall(diag)})
                break
    header = "Motif   " + "".join(f"{sf:>8}" for sf in SUBFAMS) + f"{'GLOBAL':>9}"
    print(header)
    for m in range(1, NUM_MOTIFS + 1):
        line = f"M{m:<7}"
        tot_with = tot_seqs = 0
        for sf in SUBFAMS:
            has = sum(1 for s in sub_seqs[sf] if m in s)
            n = len(sub_seqs[sf])
            line += f"{(has / n * 100 if n else 0):>7.1f}%"
            tot_with += has
            tot_seqs += n
        line += f"{(tot_with / tot_seqs * 100 if tot_seqs else 0):>8.1f}%"
        print(line)


if __name__ == "__main__":
    main()
