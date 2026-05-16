#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Auditoria de diferencias entre GFF3 Original vs Exonerate.
Compara:
1. Identidad exacta de secuencia proteica (FASTA).
2. Identidad de coordenadas genómicas (GFF3).
"""

import pandas as pd
from Bio import SeqIO
import os

# Ajusta las rutas si es necesario (asumo ejecución desde carpeta 'analisis proteinas aquaporina')
FASTA_GFF3 = 'aquaporin_peptides.fasta'
FASTA_EXO  = 'exonerate_genes_aqp.fasta'
TABLA_FILE = 'tabla_Aquaporinas_traduccion.tabular'

# Los GFF3 están en la carpeta padre
PARENT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GFF3_MAKER = os.path.join(PARENT_DIR, 'gff3_Aqp.gff3')
GFF3_EXO   = os.path.join(PARENT_DIR, 'exonerate_Aqp.gff3')

def load_fasta_to_dict(filepath):
    s = {}
    for r in SeqIO.parse(filepath, "fasta"):
        s[r.id] = str(r.seq)
    return s

def extract_cds_coords(gff_path, feature_type_filter=None):
    """
    Parsea GFF3 y extrae coordenadas de CDS agrupadas por ID padre (mRNA).
    Retorna: dict { mrna_id: [(start, end), (start, end)...] } ordenados
    """
    coords = {}
    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9: 
                continue
            
            ftype = parts[2]
            # En exonerate a veces es 'cds', en maker 'CDS'
            if ftype.lower() != 'cds':
                continue
            
            start = int(parts[3])
            end = int(parts[4])
            
            # Extraer Parent del atributo
            attr = parts[8]
            parent = None
            if 'Parent=' in attr:
                parent = attr.split('Parent=')[1].split(';')[0].split(',')[0]
            
            if parent:
                if parent not in coords:
                    coords[parent] = []
                coords[parent].append( (start, end) )
    
    # Ordenar coordenadas para comparación consistente
    for k in coords:
        coords[k].sort()
    return coords

def main():
    print("=== AUDITORIA GFF3 vs EXONERATE ===")
    
    # 1. Cargar Tabla de Pares
    print(f"Leyendo tabla: {TABLA_FILE}")
    df = pd.read_csv(TABLA_FILE, sep='\t')
    
    # Filtrar solo las que tienen ambos IDs (que son las comparables)
    pairs = df.dropna(subset=['mRNA_gff_ID', 'mRNA_exonerate_id'])
    print(f"Comparando {len(pairs)} pares secuencia/estructura.")
    
    # 2. Cargar FASTAS
    print("Cargando secuencias...")
    seqs_maker = load_fasta_to_dict(FASTA_GFF3)
    seqs_exo   = load_fasta_to_dict(FASTA_EXO)
    
    # 3. Cargar Coordenadas GFF3
    print("Extrayendo coordenadas genómicas...")
    coords_maker = extract_cds_coords(GFF3_MAKER)
    coords_exo   = extract_cds_coords(GFF3_EXO)
    
    # Contadores
    seq_match = 0
    struct_match = 0
    total = 0
    
    print("\nDetalle de diferencias (Primeros 10 casos discordantes):")
    print(f"{'Gen ID':<15} | {'Secuencia AA':<12} | {'Estructura GFF':<15} | {'Detalle'}")
    print("-" * 80)
    
    shown = 0
    
    for _, row in pairs.iterrows():
        gff_id = row['mRNA_gff_ID']
        exo_id = str(row['mRNA_exonerate_id']) # Ensure string
        
        # --- Comparar Secuencia ---
        s1 = seqs_maker.get(gff_id, "")
        s2 = seqs_exo.get(exo_id, "")
        
        # Limpiar asteriscos finales si existen
        s1 = s1.rstrip('*')
        s2 = s2.rstrip('*')
        
        is_seq_ident = (s1 == s2) and (len(s1) > 0)
        
        # --- Comparar Estructura (CDS Coords) ---
        c1 = coords_maker.get(gff_id, [])
        c2 = coords_exo.get(exo_id, [])
        
        is_struct_ident = (c1 == c2) and (len(c1) > 0)
        
        if is_seq_ident: seq_match += 1
        if is_struct_ident: struct_match += 1
        total += 1
        
        status_seq = "IGUAL" if is_seq_ident else "DIF"
        status_stru = "IGUAL" if is_struct_ident else "DIF"
        
        # Mostrar casos interesantes (Secuencia igual pero estructura diferente)
        if is_seq_ident and not is_struct_ident:
            if shown < 10:
                print(f"{row['gene_id']:<15} | {status_seq:<12} | {status_stru:<15} | Exones Maker: {c1}")
                print(f"{'':<15} | {'':<12} | {'':<15} | Exones Exo:   {c2}")
                print("-" * 80)
                shown += 1
        elif not is_seq_ident:
             # Si son diferentes en secuencia, seguro son diferentes en estructura
             pass

    
    print("\n=== RESUMEN FINAL ===")
    print(f"Total pares analizados: {total}")
    print(f"Identidad de Secuencia Proteica (FASTA): {seq_match} ({seq_match/total*100:.1f}%)")
    print(f"Identidad de Estructura Genómica  (GFF3): {struct_match} ({struct_match/total*100:.1f}%)")
    
    if seq_match > struct_match:
        print("\nCONCLUSIÓN:")
        print("Hay muchas proteínas idénticas cuyas coordenadas genómicas GFF3 NO son idénticas.")
        print("Esto explica por qué gffcompare reporta 'o' (overlap) y no '='.")
        print("Posibles causas:")
        print("1. Diferencias en UTRs (Maker suele incluirlos, Exonerate protein2genome a veces no).")
        print("2. Pequeños desplazamientos en fronteras intrón-exón que no cambian el marco de lectura (fase).")
        print("3. Exonerate fusionando exones pequeños cercanos.")

if __name__ == "__main__":
    main()
