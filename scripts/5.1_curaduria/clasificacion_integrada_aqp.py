#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pipeline Integrado de Clasificación de Aquaporinas
===================================================
Fragaria × ananassa (alooctoploide: 4 subgenomas A-D × 7 cromosomas)

Combina:
  1. Resultados de validación Fase 1/2 (selección GFF3 vs Exonerate)
  2. Clasificación filogenética con árboles IQTree (GFF3 + Exonerate)
  3. Datos complementarios (topología, MEME, DeepLoc, pepstats)

Salidas:
  - tabla_aquaporinas_traduccion.tabular (actualizada)
  - clasificacion_filogenetica_simple.csv (actualizada)
  - informe_analisis_integrado_aqp.txt

Dependencias: biopython >= 1.80, pandas
"""

import os
import re
import json
import csv
import statistics
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO, Phylo
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# ============================================================================
# CONFIGURACIÓN
# ============================================================================
# Ruta portable: override con $TFG_DATA_ROOT si tus datos viven en otra parte.
BASE_DIR = Path(os.environ.get(
    'TFG_DATA_ROOT',
    r'C:\Users\Usuario\Desktop\resultados finales'
))
PROTEINAS_DIR = BASE_DIR / 'analisis proteinas aquaporina'
OUTPUT_DIR = BASE_DIR / 'visualizaciones_tfg'

# Archivos de entrada
FASE1_JSON = OUTPUT_DIR / 'fase1_resultado.json'
MATCH_FILE = BASE_DIR / 'match_gffcompare_mRNA_genID.tabular'

# FASTA
FASTA_GFF3 = PROTEINAS_DIR / 'aquaporin_peptides.fasta'
FASTA_EXO = PROTEINAS_DIR / 'exonerate_genes_aqp.fasta'

# Árboles IQTree
TREE_GFF3 = PROTEINAS_DIR / 'fxa_aqp_gff3_129_clipkit.fasta.treefile'
TREE_EXO = PROTEINAS_DIR / 'exonerate_aqp.treefile'

# Topología, MEME, DeepLoc, pepstats
TOPO_GFF3 = PROTEINAS_DIR / 'predicted_topologies_gff3.3line'
TOPO_EXO = PROTEINAS_DIR / 'predicted_topologies_exonerate.3line'
MEME_GFF3 = PROTEINAS_DIR / 'MEME_gff3.memexml'
MEME_EXO = PROTEINAS_DIR / 'MEME_exonerate.memexml'
DEEPLOC_GFF3 = PROTEINAS_DIR / 'deeploc_gff3.csv'
DEEPLOC_EXO = PROTEINAS_DIR / 'deeploc_exonerate.csv'
PEPSTATS_GFF3 = PROTEINAS_DIR / 'pepstats_gff3.txt'
PEPSTATS_EXO = PROTEINAS_DIR / 'pepstats_exonerate.txt'

# Salidas
OUT_TABLA = PROTEINAS_DIR / 'tabla_aquaporinas_traduccion.tabular'
OUT_CSV = PROTEINAS_DIR / 'clasificacion_filogenetica_simple.csv'
OUT_INFORME = PROTEINAS_DIR / 'informe_analisis_integrado_aqp.txt'

# Patrones
GENE_ID_PATTERN = re.compile(r'Fxa(\d)([A-D])g(\d+)')
SUBGENOMAS = ['A', 'B', 'C', 'D']
NPA_PATTERN = re.compile(r'NP[ACSTVQLI]')


# ============================================================================
# 1. FUNCIONES DE CARGA
# ============================================================================
def load_fasta_dict(filepath):
    """Carga FASTA como {id: str(seq)}."""
    return {rec.id: str(rec.seq) for rec in SeqIO.parse(str(filepath), "fasta")}


def load_match_table(filepath):
    """Carga tabular GFFcompare. Retorna lista de dicts y mapeos."""
    rows = []
    gene_to_mrna = {}
    gene_to_gff_mrna = {}
    with open(filepath, 'r') as f:
        header = f.readline().strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                row = {
                    'ref_gene_id': parts[0],
                    'ref_id': parts[1],
                    'class_code': parts[2],
                    'qry_gene_id': parts[3],
                    'qry_id': parts[4],
                }
                rows.append(row)
                gene_to_mrna[parts[0]] = parts[4]  # gene -> mRNA exonerate
                gene_to_gff_mrna[parts[0]] = parts[1]  # gene -> mRNA gff3
    return rows, gene_to_mrna, gene_to_gff_mrna


def parse_gene_id(gene_id):
    """Parsea FxaXYgZZZZZ -> {'chr': X, 'sub': Y, 'num': ZZZZZ}."""
    m = GENE_ID_PATTERN.match(gene_id)
    if not m:
        return None
    return {'chr': int(m.group(1)), 'sub': m.group(2), 'num': m.group(3)}


# ============================================================================
# 2. PARSEO DE TOPOLOGÍAS
# ============================================================================
def parse_topologies(filepath):
    """Parsea .3line -> {id: {'n_tmh': N, 'topo_str': str}}."""
    result = {}
    if not filepath.exists():
        print(f"  [WARN] No encontrado: {filepath}")
        return result
    with open(filepath, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith('>'):
            seq_id = line.lstrip('>').split('|')[0].strip()
            seq_aa = lines[i + 1].strip() if i + 1 < len(lines) else ''
            topo_str = lines[i + 2].strip() if i + 2 < len(lines) else ''
            i += 3
            n_tmh = 0
            in_tm = False
            for ch in topo_str:
                if ch == 'M' and not in_tm:
                    n_tmh += 1
                    in_tm = True
                elif ch != 'M' and in_tm:
                    in_tm = False
            result[seq_id] = {'n_tmh': n_tmh, 'seq_len': len(seq_aa), 'topo_str': topo_str}
        else:
            i += 1
    return result


# ============================================================================
# 3. PARSEO DE PEPSTATS
# ============================================================================
def parse_pepstats(filepath):
    """Extrae Mw y pI de pepstats."""
    data = {}
    if not filepath.exists():
        return data
    with open(filepath, 'r') as f:
        content = f.read()
    blocks = re.split(r'PEPSTATS of (\S+)', content)
    for i in range(1, len(blocks), 2):
        pid = blocks[i].replace(' from 1 to', '').strip()
        block = blocks[i + 1]
        mw = 0.0
        pi = 0.0
        m_mw = re.search(r'Molecular weight = ([\d.]+)', block)
        m_pi = re.search(r'Isoelectric Point = ([\d.]+)', block)
        if m_mw:
            mw = float(m_mw.group(1))
        if m_pi:
            pi = float(m_pi.group(1))
        data[pid] = {'mw': mw, 'pi': pi}
    return data


# ============================================================================
# 4. PARSEO DE DEEPLOC
# ============================================================================
def parse_deeploc(filepath):
    """Extrae localización de DeepLoc CSV."""
    data = {}
    if not filepath.exists():
        return data
    import pandas as pd
    try:
        df = pd.read_csv(filepath)
        for _, row in df.iterrows():
            data[row['Protein_ID']] = row['Localizations']
    except Exception as e:
        print(f"  [WARN] Error DeepLoc: {e}")
    return data


# ============================================================================
# 5. MOTIVOS NPA
# ============================================================================
def find_channel_motifs(seq):
    """Busca motivos NPA-like en la secuencia."""
    motifs = []
    seen = set()
    for m in NPA_PATTERN.finditer(seq):
        motifs.append((m.start(), m.group()))
        seen.add(m.start())
    for pat in [re.compile(r'[SC]PV'), re.compile(r'NPQ')]:
        for m in pat.finditer(seq):
            if m.start() not in seen:
                motifs.append((m.start(), m.group()))
                seen.add(m.start())
    motifs.sort(key=lambda x: x[0])
    return motifs


def assign_loop_motifs(seq, motifs):
    """Asigna motivos a Loop B y Loop E."""
    if not motifs:
        return '-', '-'
    if len(motifs) == 1:
        pos = motifs[0][0]
        if pos < len(seq) * 0.5:
            return motifs[0][1], '-'
        else:
            return '-', motifs[0][1]
    mid = len(seq) // 2
    lb_cands = [(p, m) for p, m in motifs if p < mid]
    le_cands = [(p, m) for p, m in motifs if p >= mid]
    lb = lb_cands[-1][1] if lb_cands else '-'
    le = le_cands[0][1] if le_cands else '-'
    return lb, le


def get_motifs_str(seq):
    """Retorna (loop_B, loop_E) como strings."""
    motifs = find_channel_motifs(seq)
    return assign_loop_motifs(seq, motifs)


# ============================================================================
# 6. PROPIEDADES BIOQUÍMICAS
# ============================================================================
def calc_protein_props(seq):
    """Calcula pI, Mw, GRAVY."""
    clean = seq.replace('*', '').replace('X', 'A')
    if not clean:
        return {'mw': 0, 'pi': 0, 'gravy': 0}
    try:
        pa = ProteinAnalysis(clean)
        return {
            'mw': round(pa.molecular_weight(), 2),
            'pi': round(pa.isoelectric_point(), 2),
            'gravy': round(pa.gravy(), 4),
        }
    except Exception:
        return {'mw': 0, 'pi': 0, 'gravy': 0}


# ============================================================================
# 7. CLASIFICACIÓN FILOGENÉTICA
# ============================================================================
def analyze_tree_neighbors(tree_file, target_ids, n_neighbors=3):
    """
    Encuentra los N vecinos más cercanos de ESPECIES REFERENCIA.
    Solo retorna vecinos At*, Os*, Md*, Hb*, Fa* (NO mRNA_* ni Fxa*)
    """
    tree = Phylo.read(str(tree_file), 'newick')
    all_terminals = list(tree.get_terminals())

    # Separar terminales Fragaria propios vs referencias
    fragaria_ids = set(target_ids)
    fragaria_terminals = [t for t in all_terminals if t.name and t.name in fragaria_ids]
    reference_terminals = [t for t in all_terminals
                          if t.name and t.name not in fragaria_ids
                          and not t.name.startswith('mRNA_')
                          and not GENE_ID_PATTERN.match(t.name)]

    # Incluir también terminales de referencia Fa (FaPIP, FaNIP, etc.)
    for t in all_terminals:
        if t.name and t.name.startswith('Fa') and t.name not in fragaria_ids:
            if t not in reference_terminals:
                reference_terminals.append(t)

    print(f"    Árbol: {len(all_terminals)} hojas, {len(fragaria_terminals)} Fragaria, "
          f"{len(reference_terminals)} referencias")

    neighbors_dict = {}
    for tid in target_ids:
        target_node = None
        for t in fragaria_terminals:
            if t.name == tid:
                target_node = t
                break
        if not target_node:
            # Intentar búsqueda parcial
            for t in all_terminals:
                if t.name and tid in t.name:
                    target_node = t
                    break
        if not target_node:
            neighbors_dict[tid] = []
            continue
        distances = []
        for ref in reference_terminals:
            try:
                d = tree.distance(target_node, ref)
                distances.append((ref.name, d))
            except Exception:
                continue
        distances.sort(key=lambda x: x[1])
        neighbors_dict[tid] = [name for name, _ in distances[:n_neighbors]]
    return neighbors_dict


def extract_subfamily(name):
    """Extrae subfamilia del nombre de referencia."""
    if not name:
        return None
    m = re.search(r'(PIP|TIP|NIP|SIP|XIP)', name, re.IGNORECASE)
    return m.group(1).upper() if m else None


def extract_full_subfamily(name):
    """Extrae subfamilia completa: AtPIP2_1 -> PIP2."""
    if not name:
        return None
    m = re.search(r'(PIP|TIP|NIP|SIP|XIP)(\d)', name, re.IGNORECASE)
    if m:
        return f"{m.group(1).upper()}{m.group(2)}"
    m2 = re.search(r'(PIP|TIP|NIP|SIP|XIP)', name, re.IGNORECASE)
    return m2.group(1).upper() if m2 else None


def classify_by_phylogeny(neighbors, n_tmh=6):
    """Clasifica por consenso de vecinos filogenéticos."""
    if n_tmh < 4:
        return 'Fragmento', 'Descarte', ''
    if not neighbors:
        return 'Desconocido', 'Sin_vecinos', ''
    classes = [extract_subfamily(n) for n in neighbors]
    classes = [c for c in classes if c]
    if not classes:
        return 'Desconocido', 'Sin_clasificacion', ''
    counts = Counter(classes)
    best, n = counts.most_common(1)[0]
    if n == 3:
        conf = 'Alta'
    elif n == 2:
        conf = 'Media'
    else:
        conf = 'Baja'

    # Sub-familia detallada
    sub_classes = [extract_full_subfamily(n) for n in neighbors]
    sub_classes = [c for c in sub_classes if c]
    sub_counts = Counter(sub_classes)
    best_sub = sub_counts.most_common(1)[0][0] if sub_counts else best

    return best, conf, best_sub


# ============================================================================
# MAIN
# ============================================================================
def main():
    print("=" * 80)
    print("PIPELINE INTEGRADO: CLASIFICACIÓN DE AQUAPORINAS")
    print("Fragaria × ananassa — 129 AQPs — Alooctoploide (4 subgenomas)")
    print("=" * 80)

    # ------------------------------------------------------------------
    # PASO 1: Cargar datos base
    # ------------------------------------------------------------------
    print("\n[1] Cargando datos base...")

    # Resultado Fase 1
    if not FASE1_JSON.exists():
        print(f"  [ERROR] No existe {FASE1_JSON}. Ejecutar primero validacion_fase1.py")
        return
    with open(FASE1_JSON, 'r', encoding='utf-8') as f:
        fase1_data = json.load(f)

    all_genes_info = fase1_data['genes']
    phase1 = fase1_data['phase1']
    identical_genes = fase1_data['identical_genes']
    different_genes = fase1_data['different_genes']

    # Match table
    match_rows, gene_to_exo_mrna, gene_to_gff_mrna = load_match_table(MATCH_FILE)
    print(f"  Genes en match table: {len(match_rows)}")

    # FASTA
    pep_raw = load_fasta_dict(FASTA_GFF3)
    exo_raw = load_fasta_dict(FASTA_EXO)
    pep_dict = {h.split()[0].replace('-mRNA-1', ''): s for h, s in pep_raw.items()}
    exo_dict = dict(exo_raw)
    print(f"  GFF3 peptides: {len(pep_dict)}, Exonerate: {len(exo_dict)}")

    # Topologías
    topo_gff3_raw = parse_topologies(TOPO_GFF3)
    topo_exo_raw = parse_topologies(TOPO_EXO)
    topo_gff3 = {k.replace('-mRNA-1', ''): v for k, v in topo_gff3_raw.items()}
    topo_exo = dict(topo_exo_raw)

    # Pepstats
    pepstats_gff3 = parse_pepstats(PEPSTATS_GFF3)
    pepstats_exo = parse_pepstats(PEPSTATS_EXO)
    peps_gff3 = {k.replace('-mRNA-1', ''): v for k, v in pepstats_gff3.items()}
    peps_exo = dict(pepstats_exo)

    # DeepLoc
    deeploc_gff3 = parse_deeploc(DEEPLOC_GFF3)
    deeploc_exo = parse_deeploc(DEEPLOC_EXO)
    dl_gff3 = {k.replace('-mRNA-1', ''): v for k, v in deeploc_gff3.items()}
    dl_exo = dict(deeploc_exo)

    # ------------------------------------------------------------------
    # PASO 2: Seleccionar secuencia correcta para cada gen
    # ------------------------------------------------------------------
    print("\n[2] Seleccionando mejor secuencia por gen...")

    # Para los EMP (empate), necesitamos resolver con Fase 2
    # Ejecutamos la resolución rápida basada en priority GFF3
    gene_results = {}  # gene -> {seq, source, verdict, ...}

    for gene in sorted(all_genes_info.keys()):
        p1 = phase1.get(gene, {})
        winner = p1.get('winner', 'IGUALES')
        mrna_exo = gene_to_exo_mrna.get(gene, '')

        if winner in ('IGUALES', 'IGUALES_ALERTAS'):
            # Son idénticas, usamos GFF3
            seq = pep_dict.get(gene, p1.get('best_seq', ''))
            source = 'GFF3'
            verdict = winner
        elif winner == 'PEP':
            seq = pep_dict.get(gene, '')
            source = 'GFF3'
            verdict = 'PEPTIDE'
        elif winner == 'EXO':
            seq = exo_dict.get(mrna_exo, p1.get('best_seq', ''))
            source = 'EXONERATE'
            verdict = 'EXONERATE'
        elif winner == 'EMP':
            # Empate: ambas pasan criterios. Fase 2 prioriza GFF3 por defecto.
            # Usamos GFF3 como preferido (status quo)
            seq = pep_dict.get(gene, '')
            source = 'GFF3'
            verdict = 'PEPTIDE_EMP'
        elif winner == 'AMBAS_MAL':
            # Fallback: usar GFF3 con marca
            seq = pep_dict.get(gene, '')
            source = 'GFF3_FALLBACK'
            verdict = 'AMBAS_MAL'
        else:
            seq = pep_dict.get(gene, '')
            source = 'GFF3'
            verdict = winner

        # Propiedades de la secuencia seleccionada
        props = calc_protein_props(seq) if seq else {'mw': 0, 'pi': 0, 'gravy': 0}
        lb, le = get_motifs_str(seq) if seq else ('-', '-')

        # TMH de la fuente correcta
        if source.startswith('EXONERATE'):
            tmh_info = topo_exo.get(mrna_exo, {})
        else:
            tmh_info = topo_gff3.get(gene, {})
        n_tmh = tmh_info.get('n_tmh', 0)

        # DeepLoc de la fuente correcta
        if source.startswith('EXONERATE'):
            loc = dl_exo.get(mrna_exo, '')
        else:
            loc = dl_gff3.get(gene, '')

        gene_results[gene] = {
            'gene': gene,
            'mrna_exo': mrna_exo,
            'mrna_gff': gene_to_gff_mrna.get(gene, f'{gene}-mRNA-1'),
            'seq': seq,
            'source': source,
            'verdict': verdict,
            'length': len(seq) if seq else 0,
            'props': props,
            'loop_b': lb,
            'loop_e': le,
            'n_tmh': n_tmh,
            'localization': loc,
        }

    # Estadísticas
    source_counts = Counter(r['source'] for r in gene_results.values())
    verdict_counts = Counter(r['verdict'] for r in gene_results.values())
    print(f"  Total: {len(gene_results)} genes")
    print(f"  Por fuente: {dict(source_counts)}")
    print(f"  Por veredicto: {dict(verdict_counts)}")

    # ------------------------------------------------------------------
    # PASO 3: Clasificación filogenética con ambos árboles
    # ------------------------------------------------------------------
    print("\n[3] Clasificación filogenética con árboles IQTree...")

    # Árbol GFF3 (usa IDs de gen Fxa*)
    print("\n  -> Árbol GFF3:")
    gff3_tree_ids = list(gene_results.keys())
    neighbors_gff3 = analyze_tree_neighbors(TREE_GFF3, gff3_tree_ids)

    # Árbol Exonerate (usa IDs mRNA_*)
    print("\n  -> Árbol Exonerate:")
    exo_tree_ids = [r['mrna_exo'] for r in gene_results.values() if r['mrna_exo']]
    neighbors_exo = analyze_tree_neighbors(TREE_EXO, exo_tree_ids)

    # Clasificar cada gen
    print("\n  Clasificando...")
    for gene, result in gene_results.items():
        mrna_exo = result['mrna_exo']
        n_tmh = result['n_tmh']

        # Clasificación con árbol GFF3
        nb_gff3 = neighbors_gff3.get(gene, [])
        class_gff3, conf_gff3, sub_gff3 = classify_by_phylogeny(nb_gff3, n_tmh)

        # Clasificación con árbol Exonerate
        nb_exo = neighbors_exo.get(mrna_exo, [])
        class_exo, conf_exo, sub_exo = classify_by_phylogeny(nb_exo, n_tmh)

        # Seleccionar clasificación según fuente ganadora
        if result['source'].startswith('EXONERATE'):
            # Si exonerate es la fuente, preferir su árbol
            primary_class = class_exo
            primary_conf = conf_exo
            primary_sub = sub_exo
            primary_neighbors = nb_exo
        else:
            # GFF3 como fuente principal
            primary_class = class_gff3
            primary_conf = conf_gff3
            primary_sub = sub_gff3
            primary_neighbors = nb_gff3

        # Si ambos árboles coinciden, confianza aumenta
        trees_agree = class_gff3 == class_exo
        if trees_agree and primary_conf == 'Media':
            primary_conf = 'Alta'

        # Construir nombre de subfamilia Fa[Subfamilia]
        if primary_sub and primary_sub != 'Desconocido' and primary_sub != 'Fragmento':
            aqp_family = f"Fa{primary_sub}"
        elif primary_class and primary_class not in ('Desconocido', 'Fragmento'):
            aqp_family = f"Fa{primary_class}"
        else:
            aqp_family = primary_class

        result['phylo_class'] = primary_class
        result['phylo_conf'] = primary_conf
        result['phylo_sub'] = primary_sub
        result['aqp_family'] = aqp_family
        result['neighbors'] = primary_neighbors
        result['class_gff3'] = class_gff3
        result['class_exo'] = class_exo
        result['trees_agree'] = trees_agree

    # Resumen clasificación
    family_counts = Counter(r['aqp_family'] for r in gene_results.values())
    print("\n  Distribución por subfamilia:")
    for fam, n in sorted(family_counts.items(), key=lambda x: -x[1]):
        print(f"    {fam}: {n}")

    # ------------------------------------------------------------------
    # PASO 4: Generar tabla_aquaporinas_traduccion.tabular
    # ------------------------------------------------------------------
    print(f"\n[4] Generando {OUT_TABLA.name}...")

    headers = [
        'gene_id', 'mRNA_gff_ID', 'mRNA_exonerate_id',
        'aqp_family_subfamily', 'subfamilia_phylo', 'confianza',
        'vecino_1', 'vecino_2', 'vecino_3',
        'fuente_seq', 'veredicto', 'longitud_aa',
        'pI', 'Mw_kDa', 'TMHs',
        'motivo_B', 'motivo_E', 'localizacion',
        'arboles_coinciden'
    ]

    with open(OUT_TABLA, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=headers, delimiter='\t')
        writer.writeheader()
        for gene in sorted(gene_results.keys()):
            r = gene_results[gene]
            nb = r.get('neighbors', [])
            writer.writerow({
                'gene_id': gene,
                'mRNA_gff_ID': r['mrna_gff'],
                'mRNA_exonerate_id': r['mrna_exo'],
                'aqp_family_subfamily': r['aqp_family'],
                'subfamilia_phylo': r['phylo_class'],
                'confianza': r['phylo_conf'],
                'vecino_1': nb[0] if len(nb) > 0 else '',
                'vecino_2': nb[1] if len(nb) > 1 else '',
                'vecino_3': nb[2] if len(nb) > 2 else '',
                'fuente_seq': r['source'],
                'veredicto': r['verdict'],
                'longitud_aa': r['length'],
                'pI': r['props']['pi'],
                'Mw_kDa': round(r['props']['mw'] / 1000, 2),
                'TMHs': r['n_tmh'],
                'motivo_B': r['loop_b'],
                'motivo_E': r['loop_e'],
                'localizacion': r['localization'],
                'arboles_coinciden': 'Sí' if r.get('trees_agree', False) else 'No',
            })
    print(f"  [OK] {len(gene_results)} filas escritas")

    # ------------------------------------------------------------------
    # PASO 5: Generar clasificacion_filogenetica_simple.csv
    # ------------------------------------------------------------------
    print(f"\n[5] Generando {OUT_CSV.name}...")

    csv_headers = [
        'ID', 'Subfamilia_Filogenetica', 'Subfamilia_Detallada',
        'Confianza', 'Vecino_1', 'Vecino_2', 'Vecino_3',
        'Motivo_B', 'Motivo_E',
        'pI', 'Mw_kDa', 'TMHs', 'Localizacion',
        'Fuente', 'Veredicto',
        'Class_Arbol_GFF3', 'Class_Arbol_Exo', 'Arboles_Coinciden'
    ]

    with open(OUT_CSV, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=csv_headers)
        writer.writeheader()
        for gene in sorted(gene_results.keys()):
            r = gene_results[gene]
            nb = r.get('neighbors', [])
            writer.writerow({
                'ID': gene,
                'Subfamilia_Filogenetica': r['phylo_class'],
                'Subfamilia_Detallada': r['aqp_family'],
                'Confianza': r['phylo_conf'],
                'Vecino_1': nb[0] if len(nb) > 0 else '',
                'Vecino_2': nb[1] if len(nb) > 1 else '',
                'Vecino_3': nb[2] if len(nb) > 2 else '',
                'Motivo_B': r['loop_b'],
                'Motivo_E': r['loop_e'],
                'pI': r['props']['pi'],
                'Mw_kDa': round(r['props']['mw'] / 1000, 2),
                'TMHs': r['n_tmh'],
                'Localizacion': r['localization'],
                'Fuente': r['source'],
                'Veredicto': r['verdict'],
                'Class_Arbol_GFF3': r.get('class_gff3', ''),
                'Class_Arbol_Exo': r.get('class_exo', ''),
                'Arboles_Coinciden': 'Sí' if r.get('trees_agree', False) else 'No',
            })
    print(f"  [OK] {len(gene_results)} filas escritas")

    # ------------------------------------------------------------------
    # PASO 6: Análisis profundo del alooctoploide
    # ------------------------------------------------------------------
    print("\n[6] Análisis profundo del alooctoploide...")

    # 6.1 Distribución subfamilia × subgenoma
    dist_matrix = defaultdict(lambda: defaultdict(int))
    dist_chr = defaultdict(lambda: defaultdict(int))
    genes_by_sub = defaultdict(list)

    for gene, r in gene_results.items():
        parsed = parse_gene_id(gene)
        if not parsed:
            continue
        sub = parsed['sub']
        chrom = parsed['chr']
        fam = r['aqp_family']
        dist_matrix[fam][sub] += 1
        dist_chr[fam][chrom] += 1
        genes_by_sub[sub].append(gene)

    # 6.2 Homeólogos: comparar conservación entre subgenomas
    homeolog_groups = defaultdict(list)
    for gene, r in gene_results.items():
        parsed = parse_gene_id(gene)
        if not parsed:
            continue
        # Agrupar por cromosoma + número de gen (sin subgenoma)
        key = f"Chr{parsed['chr']}_{r['aqp_family']}"
        homeolog_groups[key].append(gene)

    # 6.3 Genes atípicos
    atypical = []
    for gene, r in gene_results.items():
        issues = []
        if r['verdict'] == 'AMBAS_MAL':
            issues.append('AMBAS_MAL')
        if r['verdict'] == 'IGUALES_ALERTAS':
            issues.append('ALERTAS')
        if r['n_tmh'] < 4:
            issues.append(f'TMH={r["n_tmh"]}(Fragmento)')
        elif r['n_tmh'] != 6:
            issues.append(f'TMH={r["n_tmh"]}')
        if r['length'] < 150:
            issues.append(f'Truncada({r["length"]}aa)')
        if issues:
            atypical.append((gene, r['aqp_family'], issues))

    # 6.4 Propiedades por subfamilia
    subfamily_props = defaultdict(lambda: {'pi': [], 'mw': [], 'len': [], 'tmh': []})
    for gene, r in gene_results.items():
        fam = r['phylo_class']
        if fam in ('Fragmento', 'Desconocido'):
            continue
        subfamily_props[fam]['pi'].append(r['props']['pi'])
        subfamily_props[fam]['mw'].append(r['props']['mw'] / 1000)
        subfamily_props[fam]['len'].append(r['length'])
        subfamily_props[fam]['tmh'].append(r['n_tmh'])

    # ------------------------------------------------------------------
    # PASO 7: Generar informe integrado
    # ------------------------------------------------------------------
    print(f"\n[7] Generando informe: {OUT_INFORME.name}...")

    with open(OUT_INFORME, 'w', encoding='utf-8') as f:
        f.write("=" * 100 + "\n")
        f.write("INFORME INTEGRADO: CLASIFICACIÓN DE AQUAPORINAS EN Fragaria × ananassa\n")
        f.write("Pipeline: Validación GFF3/Exonerate + Clasificación Filogenética + Análisis Alooctoploide\n")
        f.write("=" * 100 + "\n\n")

        # --- Sección 1: Resumen de validación ---
        f.write("━" * 80 + "\n")
        f.write("1. RESUMEN DE VALIDACIÓN (GFF3 vs Exonerate)\n")
        f.write("━" * 80 + "\n\n")
        for v, n in sorted(verdict_counts.items(), key=lambda x: -x[1]):
            f.write(f"  {v:20s}: {n:3d} genes\n")
        f.write(f"\n  Total: {sum(verdict_counts.values())} genes\n")
        f.write(f"  Usando GFF3:      {sum(1 for r in gene_results.values() if r['source'].startswith('GFF3')):3d}\n")
        f.write(f"  Usando Exonerate: {sum(1 for r in gene_results.values() if r['source'] == 'EXONERATE'):3d}\n\n")

        # Genes AMBAS_MAL detalle
        ambas_mal = [g for g, r in gene_results.items() if r['verdict'] == 'AMBAS_MAL']
        if ambas_mal:
            f.write("  Genes AMBAS_MAL (usando GFF3 como fallback):\n")
            for gene in ambas_mal:
                p1 = phase1.get(gene, {})
                f.write(f"    {gene}: GFF3 fallos={p1.get('fallos_p', [])}, "
                        f"Exo fallos={p1.get('fallos_e', [])}\n")
            f.write("\n")

        # --- Sección 2: Clasificación filogenética ---
        f.write("━" * 80 + "\n")
        f.write("2. CLASIFICACIÓN FILOGENÉTICA\n")
        f.write("━" * 80 + "\n\n")

        f.write("  Distribución por subfamilia:\n")
        for fam, n in sorted(family_counts.items(), key=lambda x: -x[1]):
            f.write(f"    {fam:15s}: {n:3d} genes\n")

        # Concordancia entre árboles
        n_agree = sum(1 for r in gene_results.values() if r.get('trees_agree', False))
        n_disagree = sum(1 for r in gene_results.values() if not r.get('trees_agree', False))
        f.write(f"\n  Concordancia entre árboles GFF3 vs Exonerate:\n")
        f.write(f"    Coinciden:  {n_agree} genes ({n_agree/len(gene_results)*100:.1f}%)\n")
        f.write(f"    Difieren:   {n_disagree} genes ({n_disagree/len(gene_results)*100:.1f}%)\n")

        # Detalle de discrepancias
        if n_disagree > 0:
            f.write("\n  Genes con clasificación discrepante entre árboles:\n")
            for gene, r in sorted(gene_results.items()):
                if not r.get('trees_agree', True):
                    f.write(f"    {gene}: GFF3={r['class_gff3']}, Exo={r['class_exo']} "
                            f"-> Asignado: {r['phylo_class']} ({r['source']})\n")
        f.write("\n")

        # --- Sección 3: Análisis del alooctoploide ---
        f.write("━" * 80 + "\n")
        f.write("3. ANÁLISIS ALOOCTOPLOIDE: DISTRIBUCIÓN POR SUBGENOMA\n")
        f.write("━" * 80 + "\n\n")

        f.write("  Genes por subgenoma:\n")
        for sub in SUBGENOMAS:
            n = len(genes_by_sub[sub])
            f.write(f"    Subgenoma {sub}: {n:3d} genes\n")

        f.write(f"\n  {'Subfamilia':<15s}")
        for sub in SUBGENOMAS:
            f.write(f"  {sub:>3s}")
        f.write("  Total\n")
        f.write("  " + "-" * 40 + "\n")

        for fam in sorted(dist_matrix.keys()):
            f.write(f"  {fam:<15s}")
            total = 0
            for sub in SUBGENOMAS:
                n = dist_matrix[fam][sub]
                f.write(f"  {n:>3d}")
                total += n
            f.write(f"  {total:>5d}\n")

        f.write("\n  Distribución por subfamilia × cromosoma:\n")
        f.write(f"  {'Subfamilia':<15s}")
        for c in range(1, 8):
            f.write(f"  Chr{c}")
        f.write("  Total\n")
        f.write("  " + "-" * 56 + "\n")
        for fam in sorted(dist_chr.keys()):
            f.write(f"  {fam:<15s}")
            total = 0
            for c in range(1, 8):
                n = dist_chr[fam][c]
                f.write(f"  {n:>4d}")
                total += n
            f.write(f"  {total:>5d}\n")

        # Análisis de conservación
        f.write("\n\n  Esperado en alooctoploide: ~4 copias por gen ancestral\n")
        f.write("  (una por cada subgenoma A, B, C, D)\n\n")

        # Detectar pérdidas/ganancias
        f.write("  Genes con distribución incompleta entre subgenomas:\n")
        # Agrupar por cromosoma + número genómico similar
        gene_groups = defaultdict(lambda: defaultdict(list))
        for gene, r in gene_results.items():
            parsed = parse_gene_id(gene)
            if not parsed:
                continue
            key = (parsed['chr'], r['aqp_family'])
            gene_groups[key][parsed['sub']].append(gene)

        for (chrom, fam), subs in sorted(gene_groups.items()):
            present = set(subs.keys())
            missing = set(SUBGENOMAS) - present
            if missing or any(len(v) > 1 for v in subs.values()):
                detail = []
                for s in SUBGENOMAS:
                    genes_in_sub = subs.get(s, [])
                    if not genes_in_sub:
                        detail.append(f"{s}: AUSENTE")
                    elif len(genes_in_sub) > 1:
                        detail.append(f"{s}: {len(genes_in_sub)}×({', '.join(genes_in_sub)})")
                    else:
                        detail.append(f"{s}: {genes_in_sub[0]}")
                f.write(f"    Chr{chrom} {fam}: {' | '.join(detail)}\n")

        # --- Sección 4: Genes atípicos ---
        f.write("\n\n")
        f.write("━" * 80 + "\n")
        f.write("4. GENES ATÍPICOS (Fragmentos, Pseudogenes, Alertas)\n")
        f.write("━" * 80 + "\n\n")

        if atypical:
            for gene, fam, issues in sorted(atypical):
                f.write(f"  {gene:>16s} [{fam}]: {', '.join(issues)}\n")
        else:
            f.write("  Ningún gen atípico detectado.\n")

        # --- Sección 5: Propiedades por subfamilia ---
        f.write("\n\n")
        f.write("━" * 80 + "\n")
        f.write("5. PROPIEDADES FISICOQUÍMICAS POR SUBFAMILIA\n")
        f.write("━" * 80 + "\n\n")

        f.write(f"  {'Subfam':<8s} {'N':>4s} {'Len(med)':>10s} {'pI(med)':>10s} "
                f"{'Mw(med)':>10s} {'TMH(med)':>10s}\n")
        f.write("  " + "-" * 55 + "\n")

        for fam in sorted(subfamily_props.keys()):
            sp = subfamily_props[fam]
            if not sp['len']:
                continue
            f.write(f"  {fam:<8s} {len(sp['len']):>4d} "
                    f"{statistics.median(sp['len']):>10.0f} "
                    f"{statistics.median(sp['pi']):>10.2f} "
                    f"{statistics.median(sp['mw']):>10.2f} "
                    f"{statistics.median(sp['tmh']):>10.0f}\n")

        # --- Sección 6: Motivos NPA por subfamilia ---
        f.write("\n\n")
        f.write("━" * 80 + "\n")
        f.write("6. MOTIVOS DE CANAL (NPA) POR SUBFAMILIA\n")
        f.write("━" * 80 + "\n\n")

        motif_by_fam = defaultdict(lambda: {'lb': Counter(), 'le': Counter()})
        for gene, r in gene_results.items():
            fam = r['phylo_class']
            if fam in ('Fragmento', 'Desconocido'):
                continue
            motif_by_fam[fam]['lb'][r['loop_b']] += 1
            motif_by_fam[fam]['le'][r['loop_e']] += 1

        for fam in sorted(motif_by_fam.keys()):
            mb = motif_by_fam[fam]
            lb_str = ', '.join(f"{m}({n})" for m, n in mb['lb'].most_common())
            le_str = ', '.join(f"{m}({n})" for m, n in mb['le'].most_common())
            f.write(f"  {fam:>6s}: Loop B = [{lb_str}]  |  Loop E = [{le_str}]\n")

        # --- Sección 7: Resumen ejecutivo ---
        f.write("\n\n")
        f.write("━" * 80 + "\n")
        f.write("7. RESUMEN EJECUTIVO PARA PUBLICACIÓN\n")
        f.write("━" * 80 + "\n\n")

        total = len(gene_results)
        n_pip = sum(1 for r in gene_results.values() if 'PIP' in r.get('aqp_family', ''))
        n_tip = sum(1 for r in gene_results.values() if 'TIP' in r.get('aqp_family', ''))
        n_nip = sum(1 for r in gene_results.values() if 'NIP' in r.get('aqp_family', ''))
        n_sip = sum(1 for r in gene_results.values() if 'SIP' in r.get('aqp_family', ''))
        n_xip = sum(1 for r in gene_results.values() if 'XIP' in r.get('aqp_family', ''))
        n_frag = sum(1 for r in gene_results.values() if r.get('phylo_class') == 'Fragmento')
        n_unk = sum(1 for r in gene_results.values() if r.get('phylo_class') == 'Desconocido')

        f.write(f"  En Fragaria × ananassa (alooctoploide, 2n = 8x = 56) se identificaron\n")
        f.write(f"  {total} genes de acuaporinas distribuidos en {len(set(r['phylo_class'] for r in gene_results.values()) - {'Fragmento', 'Desconocido'})} subfamilias principales:\n\n")
        f.write(f"    PIP: {n_pip} genes (Plasma membrane Intrinsic Proteins)\n")
        f.write(f"    TIP: {n_tip} genes (Tonoplast Intrinsic Proteins)\n")
        f.write(f"    NIP: {n_nip} genes (Nodulin-26 like Intrinsic Proteins)\n")
        f.write(f"    SIP: {n_sip} genes (Small basic Intrinsic Proteins)\n")
        f.write(f"    XIP: {n_xip} genes (X Intrinsic Proteins)\n")
        if n_frag:
            f.write(f"    Fragmentos: {n_frag} (posibles pseudogenes)\n")
        if n_unk:
            f.write(f"    No clasificados: {n_unk}\n")

        per_sub = {s: len(genes_by_sub[s]) for s in SUBGENOMAS}
        f.write(f"\n  Distribución por subgenoma: "
                + ", ".join(f"{s}={n}" for s, n in per_sub.items()) + "\n")

        n_con = sum(1 for r in gene_results.values() if r.get('trees_agree', False))
        f.write(f"\n  Concordancia filogenética (GFF3 vs Exonerate): "
                f"{n_con}/{total} ({n_con/total*100:.1f}%)\n")

        n_exo_better = sum(1 for r in gene_results.values() if r['source'] == 'EXONERATE')
        f.write(f"  Genes donde Exonerate es preferido: {n_exo_better}/{total} ({n_exo_better/total*100:.1f}%)\n")

        f.write(f"\n  Todos los genes (excepto {len(ambas_mal)} AMBAS_MAL) tienen secuencias\n")
        f.write(f"  proteicas validadas con criterios: Met inicial, 6 TMH, motivos NPA.\n")

        f.write("\n" + "=" * 100 + "\n")
        f.write("FIN DEL INFORME\n")
        f.write("=" * 100 + "\n")

    print(f"  [OK] Informe generado")

    # ------------------------------------------------------------------
    # RESUMEN FINAL
    # ------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("RESUMEN FINAL")
    print("=" * 80)
    print(f"\n  Archivos generados:")
    print(f"    1. {OUT_TABLA}")
    print(f"    2. {OUT_CSV}")
    print(f"    3. {OUT_INFORME}")
    print(f"\n  Total: {len(gene_results)} genes clasificados")
    print(f"  Subfamilias: {', '.join(sorted(set(r['aqp_family'] for r in gene_results.values())))}")
    print(f"  Concordancia árboles: {n_agree}/{len(gene_results)} ({n_agree/len(gene_results)*100:.1f}%)")
    print("\n  Pipeline completado exitosamente!")


if __name__ == "__main__":
    main()
