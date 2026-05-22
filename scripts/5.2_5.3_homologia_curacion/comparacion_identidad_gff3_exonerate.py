#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
comparacion_identidad_gff3_exonerate.py — Identidad peptidica GFF3 oficial vs Exonerate.

Reproduce la figura de validacion del curado (§6.1 del TFG): para cada gen de
acuaporina con anotacion en ambas fuentes, alinea globalmente la proteina del
GFF3 'Benihoppe' v1 frente a la de Exonerate y calcula el % de identidad.

Panel A: identidad por gen, ordenada (discrepantes vs identicos).
Panel B: distribucion de la identidad (histograma) con la media y los identicos.

El script anterior que generaba esta figura era temporal y no se conservo; este
lo sustituye de forma reproducible (lee de data/curado via config.py).

Nota tipografica: la figura NO lleva titulo embebido; el titulo va en el pie de
figura del TFG (norma APA/UCAM). Nombres de especie con x latina.

Inputs (data/curado/):
    aquaporin_peptides.fasta            proteinas del GFF3 oficial (clave = mRNA_gff_ID)
    exonerate_genes_aqp.fasta           proteinas de Exonerate    (clave = mRNA_exonerate_id)
    tabla_aquaporinas_traduccion.tabular  gene_id <-> mRNA_gff_ID <-> mRNA_exonerate_id

Output (results/):
    identidad_gff3_vs_exonerate.{png,pdf}
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.common import config

FASTA_GFF3 = config.CURADO_DIR / "aquaporin_peptides.fasta"
FASTA_EXO = config.CURADO_DIR / "exonerate_genes_aqp.fasta"
TABLA = config.CURADO_DIR / "tabla_aquaporinas_traduccion.tabular"

C_IDENT = "#3498DB"   # identicos (azul, = TIP del TFG)
C_DISCR = "#E74C3C"   # discrepantes (rojo, = PIP del TFG)

_aligner = PairwiseAligner()
_aligner.mode = "global"
_aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
_aligner.open_gap_score = -10.0
_aligner.extend_gap_score = -0.5


def load_fasta(path: Path) -> dict[str, str]:
    return {r.id: str(r.seq).rstrip("*") for r in SeqIO.parse(str(path), "fasta")}


def pct_identity(s1: str, s2: str) -> float:
    """% de identidad sobre la longitud del alineamiento global (columnas)."""
    if s1 == s2 and s1:
        return 100.0
    aln = _aligner.align(s1, s2)[0]
    matches = 0
    for (a0, a1), (b0, b1) in zip(*aln.aligned):
        matches += sum(x == y for x, y in zip(s1[a0:a1], s2[b0:b1]))
    return matches / aln.length * 100.0


def main() -> None:
    seqs_gff3 = load_fasta(FASTA_GFF3)
    seqs_exo = load_fasta(FASTA_EXO)
    df = pd.read_csv(TABLA, sep="\t")
    pairs = df.dropna(subset=["mRNA_gff_ID", "mRNA_exonerate_id"])

    rows = []
    for _, r in pairs.iterrows():
        s1 = seqs_gff3.get(str(r["mRNA_gff_ID"]), "")
        s2 = seqs_exo.get(str(r["mRNA_exonerate_id"]), "")
        if not s1 or not s2:
            continue
        rows.append({"gene_id": r["gene_id"], "pid": pct_identity(s1, s2)})

    res = pd.DataFrame(rows).sort_values("pid").reset_index(drop=True)
    ident = res[res["pid"] >= 99.999]
    discr = res[res["pid"] < 99.999]
    n_ident, n_discr, n_tot = len(ident), len(discr), len(res)
    mean_discr = discr["pid"].mean()
    print(f"Genes comparados: {n_tot}  |  identicos: {n_ident}  |  discrepantes: {n_discr}")
    print(f"Identidad media de los discrepantes: {mean_discr:.1f}%  "
          f"(rango {discr['pid'].min():.1f}-{discr['pid'].max():.1f}%, nunca 100%); "
          f"los {n_ident} identicos son 100% por igualdad exacta de secuencia.")

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(18, 7.2), dpi=150)

    # ── Panel A: identidad por gen, ordenada ─────────────────────────────────
    colors = [C_IDENT if p >= 99.999 else C_DISCR for p in res["pid"]]
    axA.bar(range(n_tot), res["pid"], color=colors, width=1.0, edgecolor="none")
    axA.axhline(100, color="#2C3E50", linestyle="--", linewidth=1.0, alpha=0.6)
    axA.axvline(n_discr - 0.5, color="#7f8c8d", linestyle=":", linewidth=1.4, alpha=0.8)
    axA.set_xlim(-0.5, n_tot - 0.5)
    axA.set_ylim(0, 105)
    axA.set_xlabel("Genes (ordenados por % identidad)", fontsize=12)
    axA.set_ylabel("% Identidad secuencia (alineamiento global)", fontsize=12)
    axA.text(0.0, 1.02, "(A)", transform=axA.transAxes, fontsize=15,
             fontweight="bold", va="bottom", ha="left", color="#1f3a5f")
    axA.legend(handles=[
        Patch(facecolor=C_DISCR, label=f"Discrepantes (n={n_discr})"),
        Patch(facecolor=C_IDENT, label=f"Idénticos (n={n_ident}, 100% identidad)"),
    ], loc="lower right", fontsize=10, framealpha=0.95)

    # ── Panel B: distribución de la identidad ────────────────────────────────
    bins = np.arange(0, 105, 5)
    axB.hist(discr["pid"], bins=bins, color=C_DISCR, edgecolor="white", linewidth=0.6)
    axB.hist(ident["pid"], bins=bins, color=C_IDENT, edgecolor="white", linewidth=0.6,
             alpha=0.85)
    axB.axvline(mean_discr, color="#C0392B", linestyle="--", linewidth=2.0,
                label="Media discrepantes: " + f"{mean_discr:.1f}%".replace(".", ","))
    axB.plot([], [], color=C_IDENT, linewidth=8, alpha=0.85,
             label=f"100% idénticos (n={n_ident})")
    axB.set_xlim(0, 102)
    axB.set_xlabel("% Identidad secuencia", fontsize=12)
    axB.set_ylabel("Nº de genes", fontsize=12)
    axB.text(0.0, 1.02, "(B)", transform=axB.transAxes, fontsize=15,
             fontweight="bold", va="bottom", ha="left", color="#1f3a5f")
    axB.legend(loc="upper left", fontsize=10, framealpha=0.95)

    fig.tight_layout()
    out_dir = config.ensure_results()
    for ext in ("png", "pdf"):
        out = out_dir / f"identidad_gff3_vs_exonerate.{ext}"
        fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
        print(f"Guardado: {out}")
    plt.close()


if __name__ == "__main__":
    main()
