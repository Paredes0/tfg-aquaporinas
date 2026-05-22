#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
regenerar_figuras.py — Regenera TODAS las figuras y tablas del TFG desde data/.

Esto NO es el pipeline pesado de RNA-seq (descarga SRA, fastp, HISAT2,
featureCounts, DESeq2): ese se ejecuta en Galaxy / servidor y produce las matrices
que ya están guardadas en `data/rna_seq/`. Este script regenera, a partir de los
datos derivados incluidos en `data/` (vía scripts/common/config.py), todas las
figuras del cuerpo (4-12), el visor eFP y las figuras/tablas de los anexos
reproducibles.

`config.py` por sí solo NO ejecuta nada: solo define las rutas. Para regenerar,
ejecuta este script:

    python scripts/regenerar_figuras.py

Salidas en `results/` (y, para los anexos, en su propia carpeta).
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PY = sys.executable
# Forzar utf-8 en los hijos para que sus prints (✓, →, etc.) no fallen al ir a un pipe.
ENV = {**os.environ, "PYTHONIOENCODING": "utf-8", "PYTHONUTF8": "1"}

# (etiqueta, ruta relativa al repo). El orden importa: profiling (Fig 6) genera
# PCA_Coordenadas_Finales.csv que necesita el Anexo C.
SCRIPTS = [
    ("Figura 4  — longitudes Rosaceae", "scripts/5.2_5.3_homologia_curacion/generar_visualizaciones_tfg.py"),
    ("Figura 5  — identidad GFF3/Exonerate", "scripts/5.2_5.3_homologia_curacion/comparacion_identidad_gff3_exonerate.py"),
    ("Figura 6  — PCA fisicoquímico (+HTML)", "scripts/5.2_5.3_homologia_curacion/profiling_final_integrated.py"),
    ("Figura 8  — MEME motivos", "scripts/5.2_5.3_homologia_curacion/analisis_motivos_final.py"),
    ("Figura 9  — TPM basal por subfamilia", "scripts/5.5_rna_seq/5.5.2_de_abundancia/figures_composition/compose_fig6_basal_subfamilia.py"),
    ("Figura 10 — expresión diferencial", "scripts/5.5_rna_seq/5.5.2_de_abundancia/figures_composition/compose_fig_de_subfamilia.py"),
    ("Figura 11 — grupos homeólogos", "scripts/5.5_rna_seq/5.5.2_de_abundancia/figures_composition/compose_fig_homeologos_basal.py"),
    ("Figura 12 — tándems NIP1", "scripts/5.5_rna_seq/5.5.2_de_abundancia/figures_composition/compose_fig_tandems_schema.py"),
    ("Figura 13 — visor eFP homeólogos (HTML)", "scripts/5.5_rna_seq/5.5.3_homeologos/15_homeolog_efp_viewer.py"),
    ("Anexo C — robustez del PCA", "annexes/C_pca_robustez/Anexo_C_script_reproducible.py"),
    ("Anexo E — soportes filogenéticos", "annexes/E_soportes_filo/Anexo_E_script_reproducible.py"),
    ("Anexo J — parámetros fisicoquímicos", "annexes/J_parametros_fisicoquimicos/Anexo_J_script_reproducible.py"),
]


def main() -> None:
    results = []
    for label, rel in SCRIPTS:
        print(f"\n{'='*72}\n>>> {label}\n    {rel}\n{'='*72}")
        proc = subprocess.run([PY, str(REPO / rel)], cwd=str(REPO),
                              capture_output=True, text=True,
                              encoding="utf-8", errors="replace", env=ENV)
        ok = proc.returncode == 0
        stream = (proc.stdout if ok else proc.stderr) or ""
        tail = (stream.strip().splitlines() or ["(sin salida)"])[-1]
        print(("    OK  " if ok else "    FALLO  ") + tail)
        results.append((label, ok))

    # ── Reunir las figuras finales en UNA carpeta, numeradas como en el TFG ──
    # (cada script las deja en su subcarpeta; aquí se copian a un índice limpio)
    dest = REPO / "results" / "figuras_TFG"
    if dest.exists():
        shutil.rmtree(dest)
    dest.mkdir(parents=True, exist_ok=True)
    coleccion = [
        ("Figura_04_longitudes_Rosaceae", ["results/visualizaciones_tfg/histograma_longitudes_aqp.png",
                                           "results/visualizaciones_tfg/histograma_longitudes_aqp.pdf"]),
        ("Figura_05_identidad_GFF3_Exonerate", ["results/identidad_gff3_vs_exonerate.png",
                                                "results/identidad_gff3_vs_exonerate.pdf"]),
        ("Figura_06_PCA_fisicoquimico", ["results/profiling_aqp_motifs_final/PCA_FINAL_INTEGRADO.png"]),
        ("Figura_08_MEME_motivos", ["results/HEATMAP_Frecuencia_Motivos_Sub-subfamilias.png"]),
        ("Figura_09_TPM_basal", ["results/figuras_rnaseq/figura6_perfiles_subfamilia.png",
                                 "results/figuras_rnaseq/figura6_perfiles_subfamilia.pdf"]),
        ("Figura_10_expresion_diferencial", ["results/figuras_rnaseq/figura_de_stripplot.png",
                                             "results/figuras_rnaseq/figura_de_stripplot.pdf"]),
        ("Figura_11_homeologos", ["results/figuras_rnaseq/figura_homeologos_basal.png",
                                  "results/figuras_rnaseq/figura_homeologos_basal.pdf"]),
        ("Figura_12_tandems_NIP1", ["results/figuras_rnaseq/figura_tandems_schema.png",
                                    "results/figuras_rnaseq/figura_tandems_schema.pdf"]),
        ("Figura_13_eFP_homeologos", ["results/efp_viewer_homeologs.html"]),
        ("PCA_interactivo", ["results/profiling_aqp_motifs_final/PCA_INTERACTIVO_FINAL.html"]),
    ]
    for nombre, fuentes in coleccion:
        for src in fuentes:
            p = REPO / src
            if p.exists():
                shutil.copy2(p, dest / f"{nombre}{p.suffix}")
    print(f"\n>>> Figuras del TFG reunidas y numeradas en: {dest}")
    print("    (la Figura 7 — árbol — se hace en iTOL; las 1-3 son de otros papers)")

    print(f"\n{'='*72}\nRESUMEN\n{'='*72}")
    for label, ok in results:
        print(f"  [{'OK   ' if ok else 'FALLO'}] {label}")
    n_ok = sum(ok for _, ok in results)
    print(f"\n{n_ok}/{len(results)} scripts ejecutados correctamente.")
    sys.exit(0 if n_ok == len(results) else 1)


if __name__ == "__main__":
    main()
