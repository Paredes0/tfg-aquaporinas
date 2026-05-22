#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Genera un archivo FASTA de consenso (consenso_aqp.fasta) seleccionando 
la mejor secuencia para cada gen basándose en la columna 'fuente_seq' 
de la tabla de traducción.

Reglas de selección y renombrado:
  - GFF3          -> Usa seq de aquaporin_peptides.fasta. Header: >gene_id
  - EXONERATE     -> Usa seq de exonerate_genes_aqp.fasta. Header: >mRNA_exonerate_id-gene_id
  - GFF3_FALLBACK -> Usa seq de aquaporin_peptides.fasta. Header: >gene_id-partial
"""

import pandas as pd
from Bio import SeqIO
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.common import config

# Entradas en el repo (data/curado/); salida en results/. Override con $TFG_DATA_ROOT.
TABLA_FILE = config.CURADO_DIR / 'tabla_aquaporinas_traduccion.tabular'
FASTA_GFF3 = config.CURADO_DIR / 'aquaporin_peptides.fasta'
FASTA_EXO  = config.CURADO_DIR / 'exonerate_genes_aqp.fasta'
FASTA_MAKER = config.CURADO_DIR / 'missing_aquaporins.fasta'

# Archivo de salida
OUTPUT_FILE = config.ensure_results() / 'consenso_aqp.fasta'

def main():
    print(f"Leyendo tabla: {TABLA_FILE}...")
    df = pd.read_csv(TABLA_FILE, sep='\t')
    
    # Verificar valores únicos en fuente_seq
    print("Valores encontrados en 'fuente_seq':", df['fuente_seq'].unique())

    # Cargar secuencias GFF3
    print(f"Cargando secuencias GFF3 de {FASTA_GFF3}...")
    seqs_gff3 = {}
    for record in SeqIO.parse(FASTA_GFF3, "fasta"):
        # El ID en el fasta suele ser 'Fxa1Ag01329-mRNA-1'
        # Lo usaremos tal cual para buscar, asumiendo que coincide con mRNA_gff_ID
        seqs_gff3[record.id] = str(record.seq)

    # Cargar secuencias EXONERATE
    print(f"Cargando secuencias Exonerate de {FASTA_EXO}...")
    seqs_exo = {}
    for record in SeqIO.parse(FASTA_EXO, "fasta"):
        # El ID suele ser 'mRNA_375'
        seqs_exo[record.id] = str(record.seq)

    # Cargar secuencias MAKER (GDR)
    print(f"Cargando secuencias MAKER de {FASTA_MAKER}...")
    seqs_maker = {}
    if os.path.exists(FASTA_MAKER):
        for record in SeqIO.parse(FASTA_MAKER, "fasta"):
            gid = record.id.replace('-mRNA-1_Benihoppe_v1', '')
            seqs_maker[gid] = str(record.seq)

    print(f"Generando {OUTPUT_FILE}...")
    
    count_gff3 = 0
    count_exo = 0
    count_fallback = 0
    count_maker = 0
    count_missing = 0
    
    with open(OUTPUT_FILE, 'w') as f_out:
        for index, row in df.iterrows():
            gene_id = str(row['gene_id'])
            fuente = str(row['fuente_seq']).strip()
            
            # IDs de referencia
            mrna_gff_id = str(row['mRNA_gff_ID']) if pd.notna(row['mRNA_gff_ID']) else None
            mrna_exo_id = str(row['mRNA_exonerate_id']) if pd.notna(row['mRNA_exonerate_id']) else None
            
            sequence = None
            header = None
            
            if fuente == 'GFF3':
                if mrna_gff_id and mrna_gff_id in seqs_gff3:
                    sequence = seqs_gff3[mrna_gff_id]
                    header = gene_id
                    count_gff3 += 1
                else:
                    print(f"[WARNING] Secuencia GFF3 no encontrada para {gene_id} (ID buscado: {mrna_gff_id})")

            elif fuente == 'EXONERATE':
                if mrna_exo_id and mrna_exo_id in seqs_exo:
                    sequence = seqs_exo[mrna_exo_id]
                    header = f"{mrna_exo_id}-{gene_id}"
                    count_exo += 1
                else:
                    print(f"[WARNING] Secuencia Exonerate no encontrada para {gene_id} (ID buscado: {mrna_exo_id})")

            elif fuente == 'GFF3_FALLBACK':
                if mrna_gff_id and mrna_gff_id in seqs_gff3:
                    sequence = seqs_gff3[mrna_gff_id]
                    header = f"{gene_id}-partial"
                    count_fallback += 1
                else:
                    print(f"[WARNING] Secuencia GFF3_FALLBACK no encontrada para {gene_id} (ID buscado: {mrna_gff_id})")
            
            elif fuente == 'MAKER_GFF3':
                if gene_id in seqs_maker:
                    sequence = seqs_maker[gene_id]
                    tmhs = row.get('TMHs', 0)
                    if pd.notna(tmhs) and int(tmhs) >= 6:
                        header = f"{gene_id}-maker"
                    else:
                        header = f"{gene_id}-partial-maker"
                    count_maker += 1
                else:
                    print(f"[WARNING] Secuencia MAKER no encontrada para {gene_id}")

            else:
                # Si hay algún otro estado o nan, informar
                pass

            if sequence and header:
                f_out.write(f">{header}\n{sequence}\n")
            else:
                count_missing += 1

    print("-" * 30)
    print("Resumen de secuencias escritas:")
    print(f"  GFF3:          {count_gff3}")
    print(f"  EXONERATE:     {count_exo}")
    print(f"  GFF3_FALLBACK: {count_fallback}")
    print(f"  MAKER_GFF3:    {count_maker}")
    print(f"  TOTAL:         {count_gff3 + count_exo + count_fallback + count_maker}")
    print("-" * 30)
    if count_missing > 0:
        print(f"[ATENCIÓN] Hubo {count_missing} filas sin asignar (puede ser por fuente desconocida o secuencia no encontrada).")

if __name__ == "__main__":
    main()
