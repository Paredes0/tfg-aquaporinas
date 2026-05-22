#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Anexo_J_script_reproducible.py — Tabla de parámetros fisicoquímicos por acuaporina.

Genera la tabla por gen de las 121 acuaporinas funcionales con los parámetros
citados en el TFG (pI, peso molecular, nº de hélices transmembrana) y las
variantes de los dos motivos NPA (bucles B y E), que sustentan, entre otras, la
discusión de las SIP (NPL/NPT/NPS en el bucle B).

Fuente: data/curado/tabla_aqp_ordenada.csv (salida consolidada del curado).

Output (en la carpeta del anexo):
    Anexo_J_parametros_fisicoquimicos.csv
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.common import config

SRC = config.CURADO_DIR / "tabla_aqp_ordenada.csv"
OUT = Path(__file__).resolve().parent / "Anexo_J_parametros_fisicoquimicos.csv"

SUBFAM_ORDER = {"PIP": 0, "TIP": 1, "NIP": 2, "SIP": 3, "XIP": 4}

COLUMNS = [
    ("gene_id", "gene_id"),
    ("subfamilia", "subfamilia_phylo"),
    ("sub_subfamilia", "aqp_family_subfamily"),
    ("longitud_aa", "longitud_aa"),
    ("pI", "pI"),
    ("Mw_kDa", "Mw_kDa"),
    ("TMHs", "TMHs"),
    ("NPA_bucle_B", "motivo_B"),
    ("NPA_bucle_E", "motivo_E"),
    ("localizacion", "localizacion"),
]


def main() -> None:
    rows = list(csv.DictReader(open(SRC, encoding="utf-8")))
    rows.sort(key=lambda r: (SUBFAM_ORDER.get(r["subfamilia_phylo"], 9),
                             r["aqp_family_subfamily"], r["gene_id"]))

    with open(OUT, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([dst for dst, _ in COLUMNS])
        for r in rows:
            w.writerow([r.get(src, "") for _, src in COLUMNS])

    # Resumen a stdout
    from collections import Counter
    sub = Counter(r["subfamilia_phylo"] for r in rows)
    npb = Counter(r["motivo_B"] for r in rows if r["subfamilia_phylo"] == "SIP")
    print(f"Acuaporinas en la tabla: {len(rows)}  ({dict(sub)})")
    print(f"Motivo NPA del bucle B en SIP: {dict(npb)}")
    print(f"Tabla guardada: {OUT}")


if __name__ == "__main__":
    main()
