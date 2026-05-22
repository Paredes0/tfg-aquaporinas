#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Anexo_I_script_reproducible.py — Distribución de soportes filogenéticos.

Recorre el árbol filogenético definitivo (Q.PLANT+R6, 282 hojas, log L = -45.149,26)
y extrae los tres estadísticos de soporte que IQ-TREE anota en cada nodo interno:
SH-aLRT (1000 iteraciones, %), aBayes y ultrafast bootstrap UFBoot (1000 iter, %).

Genera tres salidas en la misma carpeta:
  - Anexo_I_tabla_277_nodos_soportes.csv      (un nodo por fila + indicadores de umbral)
  - Anexo_I_figura_histograma_soportes.png    (3 paneles con líneas de umbral)
  - Anexo_I_tabla_nodos_subfamilias.csv       (soporte del nodo que define cada subfamilia)

El árbol se localiza vía la variable de entorno TFG_RNA_SEQ_ROOT o, en su defecto,
en las rutas locales conocidas. Formato de soporte en el .treefile:
    )SH-aLRT/aBayes/UFboot:branch_length     (p. ej.  )98.9/1/100:0.23 )

Umbrales convencionales de confianza alta:
    SH-aLRT >= 80 %    |    aBayes >= 0,95    |    UFBoot >= 95 %
"""
from __future__ import annotations

import csv
import os
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Localización del árbol final ─────────────────────────────────────────────
# El árbol de referencia del TFG (final_without_partials: 282 hojas, 430 sitios,
# Q.PLANT+R6, log L = -45.149,26 → 277 nodos internos) viene incluido en el repo,
# en data/filogenia/; las cifras de soporte de §6.2.2 proceden de ahí.
# Override opcional con $TFG_TREEFILE para apuntar a otra copia.
_REPO_TREE = str(Path(__file__).resolve().parents[2] / "data" / "filogenia" / "arbol_acuaporinas.treefile")
_CANDIDATE_TREES = [
    os.environ.get("TFG_TREEFILE", ""),
    _REPO_TREE,  # copia incluida en el repo (data/filogenia/) → autorreproducible
]

# Mapeo gene_id -> subfamilia para las hojas de Fragaria (homeolog_groups.tsv del repo).
_REPO_HOMEOLOGS = str(Path(__file__).resolve().parents[2] / "data" / "rna_seq" / "homeologos" / "homeolog_groups.tsv")
_CANDIDATE_HOMEOLOGS = [
    os.environ.get("TFG_HOMEOLOGS", ""),
    _REPO_HOMEOLOGS,  # incluido en el repo → autorreproducible
]

OUT_DIR = Path(__file__).resolve().parent

# Umbrales convencionales
THR_SHALRT = 80.0   # %
THR_ABAYES = 0.95
THR_UFBOOT = 95.0   # %

SUBFAMILIES = ["PIP", "TIP", "NIP", "SIP", "XIP"]
# Paleta oficial del TFG (para la barra de la subtabla, coherencia visual)
SF_COLORS = {"PIP": "#E74C3C", "TIP": "#3498DB", "NIP": "#2ECC71",
             "SIP": "#F39C12", "XIP": "#9B59B6"}


def _first_existing(paths: list[str]) -> str | None:
    for p in paths:
        if os.path.exists(p):
            return p
    return None


def parse_node_supports(newick: str) -> list[tuple[float, float, float]]:
    """Extrae (SH-aLRT, aBayes, UFBoot) de cada nodo interno anotado."""
    pattern = re.compile(r"\)([0-9.]+)/([0-9.]+)/([0-9.]+):")
    out = []
    for sh, ab, uf in pattern.findall(newick):
        out.append((float(sh), float(ab), float(uf)))
    return out


def write_node_table(supports: list[tuple[float, float, float]], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow([
            "node_id", "SH_aLRT", "aBayes", "UFBoot",
            "pasa_SH_aLRT_80", "pasa_aBayes_0.95", "pasa_UFBoot_95", "triple_soporte_alto",
        ])
        for i, (sh, ab, uf) in enumerate(supports, start=1):
            p_sh = sh >= THR_SHALRT
            p_ab = ab >= THR_ABAYES
            p_uf = uf >= THR_UFBOOT
            w.writerow([
                f"N{i:03d}", f"{sh:g}", f"{ab:g}", f"{uf:g}",
                int(p_sh), int(p_ab), int(p_uf), int(p_sh and p_ab and p_uf),
            ])


def plot_histograms(supports: list[tuple[float, float, float]], path: Path) -> None:
    sh = [s[0] for s in supports]
    ab = [s[1] for s in supports]
    uf = [s[2] for s in supports]

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))

    panels = [
        (axes[0], sh, "SH-aLRT (%)", THR_SHALRT, "#34495E", (0, 100)),
        (axes[1], ab, "aBayes", THR_ABAYES, "#16A085", (0, 1)),
        (axes[2], uf, "UFBoot (%)", THR_UFBOOT, "#C0392B", (0, 100)),
    ]
    for ax, data, label, thr, color, xlim in panels:
        ax.hist(data, bins=20, range=xlim, color=color, alpha=0.85, edgecolor="white")
        ax.axvline(thr, color="black", linestyle="--", linewidth=1.3)
        n_pass = sum(1 for v in data if v >= thr)
        ax.set_title(f"{label}\n{n_pass}/{len(data)} nodos por encima del umbral", fontsize=10)
        ax.set_xlabel(label)
        ax.set_ylabel("Nº de nodos internos")
        ax.set_xlim(xlim)

    n_triple = sum(1 for s, a, u in supports if s >= THR_SHALRT and a >= THR_ABAYES and u >= THR_UFBOOT)
    fig.suptitle(
        f"Distribución de soportes sobre {len(supports)} nodos internos "
        f"(triple soporte alto: {n_triple} nodos, {100*n_triple/len(supports):.1f} %)",
        fontsize=11, y=1.02,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def load_leaf_subfamilies(homeolog_path: str | None) -> dict[str, str]:
    """Mapea cada hoja del árbol a su subfamilia.

    - Referencias (At/Os/Md/Hb/...): subfamilia por regex en el nombre.
    - Fragaria (Fxa...): subfamilia desde homeolog_groups.tsv.
    """
    mapping: dict[str, str] = {}
    if homeolog_path and os.path.exists(homeolog_path):
        with open(homeolog_path, encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                gid = row.get("gene_id", "").strip()
                # La columna 'subfamily' contiene la sub-subfamilia (p. ej. PIP2,
                # TIP1, NIP1); se extrae la subfamilia canónica de los 3 primeros
                # caracteres.
                raw_sf = row.get("subfamily", "").strip().upper()
                m = re.match(r"(PIP|TIP|NIP|SIP|XIP)", raw_sf)
                if gid and m:
                    mapping[gid] = m.group(1)
    return mapping


def leaf_subfamily(leaf: str, fxa_map: dict[str, str]) -> str | None:
    """Devuelve la subfamilia de una hoja, o None si no se puede asignar."""
    # Hoja de Fragaria: puede venir como 'mRNA_xxxx-FxaYYYY' o 'FxaYYYY'
    m = re.search(r"(Fxa[0-9][A-Za-z0-9]+)", leaf)
    if m:
        return fxa_map.get(m.group(1))
    # Hoja de referencia: subfamilia en el propio nombre
    m = re.search(r"(PIP|TIP|NIP|SIP|XIP)", leaf)
    return m.group(1) if m else None


def write_subfamily_table(tree_path: str, fxa_map: dict[str, str], path: Path) -> None:
    """Soporte del nodo que define (MRCA) cada subfamilia canónica."""
    from Bio import Phylo  # import local: solo se necesita aquí

    tree = Phylo.read(tree_path, "newick")
    leaves = tree.get_terminals()

    # Agrupar hojas por subfamilia
    by_sf: dict[str, list] = {sf: [] for sf in SUBFAMILIES}
    for lf in leaves:
        sf = leaf_subfamily(lf.name, fxa_map)
        if sf in by_sf:
            by_sf[sf].append(lf)

    rows = []
    for sf in SUBFAMILIES:
        members = by_sf[sf]
        if len(members) < 2:
            rows.append([sf, len(members), "n/d", "n/d", "n/d", "subfamilia con <2 hojas mapeadas"])
            continue
        mrca = tree.common_ancestor(members)
        # En Bio.Phylo el label compuesto SH/aB/UF queda en confidence o name
        raw = mrca.confidence if mrca.confidence is not None else mrca.name
        raw = str(raw) if raw is not None else ""
        m = re.match(r"([0-9.]+)/([0-9.]+)/([0-9.]+)", raw)
        if m:
            sh, ab, uf = m.group(1), m.group(2), m.group(3)
            triple = (float(sh) >= THR_SHALRT and float(ab) >= THR_ABAYES and float(uf) >= THR_UFBOOT)
            nota = "MRCA monofilético" if mrca.count_terminals() == len(members) else \
                   f"clado MRCA contiene {mrca.count_terminals()} hojas ({len(members)} de la subfamilia)"
            rows.append([sf, len(members), sh, ab, uf, ("triple alto" if triple else nota)])
        elif mrca is tree.root:
            rows.append([sf, len(members), "n/d", "n/d", "n/d",
                         "clado basal: su MRCA coincide con la raíz del árbol no "
                         "enraizado, punto que no porta soporte"])
        else:
            rows.append([sf, len(members), "n/d", "n/d", "n/d", "nodo sin soporte anotado"])

    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["subfamilia", "n_hojas_mapeadas", "SH_aLRT", "aBayes", "UFBoot", "observacion"])
        w.writerows(rows)


def main() -> None:
    tree_path = _first_existing(_CANDIDATE_TREES)
    if tree_path is None:
        raise SystemExit(
            "No se encontró el .treefile final. Define TFG_RNA_SEQ_ROOT o ajusta las rutas."
        )
    print(f"[Anexo E] Árbol: {tree_path}")
    newick = Path(tree_path).read_text(encoding="utf-8")

    supports = parse_node_supports(newick)
    print(f"[Anexo E] Nodos internos con triple soporte: {len(supports)}")

    write_node_table(supports, OUT_DIR / "Anexo_I_tabla_277_nodos_soportes.csv")
    plot_histograms(supports, OUT_DIR / "Anexo_I_figura_histograma_soportes.png")

    fxa_map = load_leaf_subfamilies(_first_existing(_CANDIDATE_HOMEOLOGS))
    print(f"[Anexo E] Hojas Fragaria mapeadas a subfamilia: {len(fxa_map)}")
    write_subfamily_table(tree_path, fxa_map, OUT_DIR / "Anexo_I_tabla_nodos_subfamilias.csv")

    print("[Anexo E] Salidas generadas en:", OUT_DIR)


if __name__ == "__main__":
    main()
