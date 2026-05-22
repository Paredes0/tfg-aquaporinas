#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Anexo_H_script_reproducible.py — Tabla de parámetros fisicoquímicos por acuaporina.

Genera la tabla por gen de las 121 acuaporinas funcionales con los parámetros
citados en el TFG (pI, peso molecular, nº de hélices transmembrana) y las
variantes de los dos motivos NPA (bucles B y E), que sustentan, entre otras, la
discusión de las SIP (NPL/NPT/NPS en el bucle B).

Fuente: data/curado/tabla_aqp_ordenada.csv (salida consolidada del curado).

Output (en la carpeta del anexo):
    Anexo_H_parametros_fisicoquimicos.csv
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.common import config

SRC = config.CURADO_DIR / "tabla_aqp_ordenada.csv"
HOM = config.RNASEQ_HOM_DIR / "homeolog_groups.tsv"   # isoforma (sub_subfamily con ;N)
OUT = Path(__file__).resolve().parent / "Anexo_H_parametros_fisicoquimicos.csv"

SUBFAM_ORDER = {"PIP": 0, "TIP": 1, "NIP": 2, "SIP": 3, "XIP": 4}
MODELO = {"GFF3": "GFF3 oficial", "EXONERATE": "Exonerate (sustituido)"}

# (columna_salida, columna_origen en tabla_aqp_ordenada)
COLUMNS = [
    ("gene_id", "gene_id"),
    ("subfamilia", "subfamilia_phylo"),
    ("isoforma", None),          # se rellena desde homeolog_groups (sub_subfamily)
    ("modelo", None),            # se rellena desde fuente_seq
    ("mRNA_modelo", None),       # ID del mRNA realmente usado (GFF3 o Exonerate)
    ("longitud_aa", "longitud_aa"),
    ("pI", "pI"),
    ("Mw_kDa", "Mw_kDa"),
    ("TMHs", "TMHs"),
    ("NPA_bucle_B", "motivo_B"),
    ("NPA_bucle_E", "motivo_E"),
    ("localizacion", "localizacion"),
]


def load_isoformas(path: Path) -> dict[str, str]:
    """gene_id -> isoforma (sub_subfamily, p. ej. SIP1;3) desde homeolog_groups.tsv."""
    iso = {}
    with open(path, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            iso[r["gene_id"]] = r.get("sub_subfamily", "")
    return iso


def main() -> None:
    rows = list(csv.DictReader(open(SRC, encoding="utf-8")))
    iso_map = load_isoformas(HOM)

    rows.sort(key=lambda r: (SUBFAM_ORDER.get(r["subfamilia_phylo"], 9),
                             r["aqp_family_subfamily"], r["gene_id"]))

    def value(r, dst, src):
        if dst == "isoforma":
            # isoforma del grupo homeólogo; si falta, sub-subfamilia de la tabla
            return iso_map.get(r["gene_id"]) or r["aqp_family_subfamily"]
        if dst == "modelo":
            return MODELO.get(r["fuente_seq"], r["fuente_seq"])
        if dst == "mRNA_modelo":
            return r["mRNA_exonerate_id"] if r["fuente_seq"] == "EXONERATE" else r["mRNA_gff_ID"]
        return r.get(src, "")

    with open(OUT, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([dst for dst, _ in COLUMNS])
        for r in rows:
            w.writerow([value(r, dst, src) for dst, src in COLUMNS])

    # Resumen / verificación a stdout
    from collections import Counter
    sub = Counter(r["subfamilia_phylo"] for r in rows)
    fuente = Counter(r["fuente_seq"] for r in rows)
    npb_sip = Counter(r["motivo_B"] for r in rows if r["subfamilia_phylo"] == "SIP")
    npe_sip = Counter(r["motivo_E"] for r in rows if r["subfamilia_phylo"] == "SIP")
    sin_iso = [r["gene_id"] for r in rows if not iso_map.get(r["gene_id"])]
    print(f"Acuaporinas en la tabla: {len(rows)}  ({dict(sub)})")
    print(f"Modelo de secuencia: {dict(fuente)}  (EXONERATE = sustituido)")
    print(f"NPA bucle B en SIP: {dict(npb_sip)}  |  NPA bucle E en SIP: {dict(npe_sip)}")
    print(f"Genes sin isoforma en homeolog_groups (usan sub-subfamilia): {len(sin_iso)} {sin_iso}")
    print(f"Tabla guardada: {OUT}")


if __name__ == "__main__":
    main()
