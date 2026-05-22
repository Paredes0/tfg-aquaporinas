#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
    PROFILING AQUAPORINAS: INTEGRATED GOLD STANDARD (v4.0)
    - Física Profunda (TMHMM/Pepstats)
    - Huella Estadística MEME (M1-M12) desde FASTA combinado unificado
    - Elipses de Confianza Geométricamente Exactas
===============================================================================
"""

import warnings
warnings.filterwarnings('ignore')

import sys, io, os, re, json
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Ellipse
from scipy.stats import zscore

from sklearn.preprocessing import StandardScaler, MultiLabelBinarizer
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import plotly.express as px

# ─── CONFIGURACIÓN ───────────────────────────────────────────────────────────
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))
from scripts.common import config as _cfg

# Inputs: data/curado/ del repo (los datos GDR viven en la misma carpeta).
BASE_DIR = str(_cfg.CURADO_DIR)
GDR_DIR = str(_cfg.GDR_DIR)
# Outputs: results/ del repo (no versionado).
OUT = str(_cfg.RESULTS_DIR / 'profiling_aqp_motifs_final')
os.makedirs(OUT, exist_ok=True)
COLOR_MAP = {
    'PIP': '#E74C3C', 'TIP': '#3498DB', 'NIP': '#2ECC71',
    'SIP': '#F39C12', 'XIP': '#9B59B6', 'Fragmento': '#95A5A6',
    'PARTIAL': '#D3D3D3', 'CURATED': '#808B96'
}

def banner(title):
    print("\n" + "="*80 + "\n" + f" {title} ".center(80, " ") + "\n" + "="*80)

# ─── PASO 1: MAPEO Y HUELLA DE MOTIVOS ───────────────────────────────────────
banner("PASO 1: GENERACIÓN DE HUELLA ESTRUCTURAL (M1-M12)")

df_tab = pd.read_csv(os.path.join(BASE_DIR, 'tabla_aquaporinas_traduccion.tabular'), sep='\t')
exo_to_gene = {str(r['mRNA_exonerate_id']): r['gene_id'] for _, r in df_tab.iterrows() if pd.notna(r.get('mRNA_exonerate_id'))}

# ── Construir el mapeo: para cada gene_id, cuál es la secuencia seleccionada ──
# En el FASTA combinado (exonerate_gff3_aqp.fasta) hay 258 secuencias:
#   - 129 con IDs FxaXXX (del GFF3 oficial)
#   - 129 con IDs mRNA_XXX (de Exonerate)
# El archivo MEME analiza todas las 258 juntas, asegurando motivos comparables.
# Para cada gen, solo tomamos los motivos de la secuencia que fuente_seq indica.
selected_fasta_id = {}  # gene_id -> ID tal cual aparece en el FASTA/MEME
for _, r in df_tab.iterrows():
    gid = r['gene_id']
    src = r['fuente_seq']
    eid = str(r['mRNA_exonerate_id']) if pd.notna(r.get('mRNA_exonerate_id')) else None
    if src == 'EXONERATE' and eid:
        selected_fasta_id[gid] = eid       # ej: mRNA_54367
    else:
        selected_fasta_id[gid] = gid        # ej: Fxa2Ag01701

# Conjunto inverso: FASTA_id -> gene_id (para las seleccionadas)
fasta_to_gene = {}
for gid, fid in selected_fasta_id.items():
    fasta_to_gene[fid] = gid

def parse_meme_combined_block(filepath):
    """
    Parsea la tabla 'Combined block diagrams' del archivo MEME.
    Devuelve un dict: {fasta_id: set(motifs)} para TODAS las 258 secuencias.
    No traduce IDs aquí; eso se hace después con el mapeo de fuente_seq.
    """
    if not os.path.exists(filepath):
        print(f"  [WARN] No se encontró: {filepath}")
        return {}
    res = {}
    in_block = False
    current_sid = None
    accumulated_diagram = ""
    
    with open(filepath, 'r') as f:
        for line in f:
            if 'Combined block diagrams: non-overlapping sites with p-value < 0.0001' in line:
                in_block = True
                continue
            if in_block:
                if line.startswith('---'): continue
                if line.startswith('SEQUENCE NAME'): continue
                if line.startswith('***') or 'CPU:' in line: break
                
                line = line.strip()
                if not line: continue
                
                parts = line.split()
                # Nueva línea de secuencia (nombre + p-value + diagrama)
                if len(parts) >= 3 and not line.startswith('['):
                    sid = parts[0].replace('-mRNA-1', '').replace('-partial', '')
                    current_sid = sid
                    diagram_part = parts[2]
                else:
                    # Línea de continuación
                    diagram_part = parts[0]
                
                if diagram_part.endswith('\\'):
                    accumulated_diagram += diagram_part[:-1]
                else:
                    accumulated_diagram += diagram_part
                    if current_sid:
                        motifs = set(f"M{m}" for m in re.findall(r'\[(\d+)\(', accumulated_diagram))
                        if current_sid not in res:
                            res[current_sid] = set()
                        res[current_sid].update(motifs)
                    current_sid = None
                    accumulated_diagram = ""
    return res

# ── Parsear el archivo MEME del FASTA combinado ──
# IMPORTANTE: Este archivo proviene de correr MEME sobre exonerate_gff3_aqp.fasta
# que contiene las 258 secuencias (129 GFF3 + 129 Exonerate).
def parse_fimo_missing(filepath):
    if not os.path.exists(filepath): return {}
    res = {}
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('motif_id') or line.startswith('#') or not line.strip(): continue
            parts = line.split('\t')
            if len(parts) >= 3:
                alt_id = parts[1]
                seq_name = parts[2].replace('-mRNA-1_Benihoppe_v1', '')
                m_match = re.search(r'MEME-(\d+)', alt_id)
                if m_match:
                    motif = f"M{m_match.group(1)}"
                    if seq_name not in res: res[seq_name] = set()
                    res[seq_name].add(motif)
    return res

MEME_COMBINED = os.path.join(BASE_DIR, 'MEME_exonerate_gff3_aqp.txt')
all_motifs = parse_meme_combined_block(MEME_COMBINED)
FIMO_MISSING = os.path.join(GDR_DIR, 'motifs_aqp_missing.tsv')
all_motifs.update(parse_fimo_missing(FIMO_MISSING))
print(f"  Secuencias parseadas del MEME (+FIMO): {len(all_motifs)}")

# ── Seleccionar solo los motivos de la secuencia elegida por fuente_seq ──
presence_map = {}
mapped_count = 0
missing_ids = []
for gid, fasta_id in selected_fasta_id.items():
    if fasta_id in all_motifs:
        presence_map[gid] = all_motifs[fasta_id]
        mapped_count += 1
    else:
        missing_ids.append((gid, fasta_id))
        presence_map[gid] = set()  # Sin motivos asignados

print(f"  Genes mapeados correctamente: {mapped_count}/{len(selected_fasta_id)}")
if missing_ids:
    print(f"  [WARN] Genes sin match en MEME ({len(missing_ids)}):")
    for gid, fid in missing_ids[:10]:
        print(f"    {gid} -> esperaba '{fid}' en MEME")

motif_cols = [f"M{i}" for i in range(1, 13)]
motif_rows = []
for gid, found in presence_map.items():
    row = {'ID': gid}
    for m in motif_cols: row[m] = 1 if m in found else 0
    row['Total_Motivos'] = sum(row[m] for m in motif_cols)
    # Score Diagnóstico: Suma de motivos clave por clado
    row['Score_Diagnostico'] = sum(row[m] for m in ['M5', 'M7', 'M9', 'M10'])
    motif_rows.append(row)
df_motifs = pd.DataFrame(motif_rows)

# ── Informe rápido de distribución de motivos ──
print(f"\n  Distribución de motivos en las {len(df_motifs)} secuencias seleccionadas:")
for m in motif_cols:
    count = df_motifs[m].sum()
    print(f"    {m}: {int(count)} secuencias ({count/len(df_motifs)*100:.0f}%)")

# ─── PASO 2: FÍSICA Y TOPOLOGÍA (TMHMM + Pepstats + DeepLoc) ─────────────────
banner("PASO 2: INTEGRACIÓN DE FÍSICA, TOPOLOGÍA Y LOCALIZACIÓN PROFUNDA")

# ── TMHMM: 5 features de topología ──
def parse_tmhmm(path):
    if not os.path.exists(path): return {}
    res = {}
    with open(path, 'r') as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        if lines[i].startswith('>'):
            gid = lines[i].split()[0].replace('>', '').replace('-mRNA-1', '').replace('_Benihoppe_v1', '')
            gid = exo_to_gene.get(gid, gid)
            if i + 2 < len(lines):
                topo = lines[i+2].strip()
                segments = re.findall(r'M+', topo)
                res[gid] = {
                    'TMHMM_Segments': len(segments),
                    'TMHMM_Avg_TM_Len': round(np.mean([len(s) for s in segments]), 1) if segments else 0,
                    'TMHMM_Frac_Inside': round(topo.count('I') / len(topo), 4) if topo else 0,
                    'TMHMM_Frac_Outside': round(topo.count('O') / len(topo), 4) if topo else 0,
                    'TMHMM_Frac_TM': round(topo.count('M') / len(topo), 4) if topo else 0
                }
            i += 3
        else:
            i += 1
    return res

# ── Pepstats EMBOSS: 13 features de composición ──
def parse_pepstats(path):
    if not os.path.exists(path): return {}
    res = {}
    with open(path, 'r') as f:
        content = f.read()
    blocks = re.split(r'PEPSTATS of (\S+)', content)
    for i in range(1, len(blocks)-1, 2):
        sid = blocks[i].replace('-mRNA-1', '').replace('_Benihoppe_v1', '')
        gid = exo_to_gene.get(sid, sid)
        block = blocks[i+1]
        data = {}
        m = re.search(r'Charge\s*=\s*([-\d.]+)', block)
        data['PS_Charge'] = float(m.group(1)) if m else 0
        m = re.search(r'A280 Extinction Coefficient 1mg/ml\s*=\s*([\d.]+)', block)
        data['PS_A280_1mg'] = float(m.group(1)) if m else 0
        m = re.search(r'Improbability of expression in inclusion bodies\s*=\s*([\d.]+)', block)
        data['PS_Inclusion_Body'] = float(m.group(1)) if m else 0
        m = re.search(r'Average Residue Weight\s*=\s*([\d.]+)', block)
        data['PS_Avg_Residue_Wt'] = float(m.group(1)) if m else 0
        for prop in ['Tiny', 'Small', 'Aliphatic', 'Aromatic', 'Non-polar',
                     'Polar', 'Charged', 'Basic', 'Acidic']:
            m = re.search(rf'{prop}\s+\([^)]+\)\s+\d+\s+([\d.]+)', block)
            data[f'PS_{prop}_Pct'] = float(m.group(1)) if m else 0
        res[gid] = data
    return res

tm_data_gff3 = {}
for f in ['predicted_topologies_gff3.3line']:
    fp = os.path.join(BASE_DIR, f)
    if os.path.exists(fp):
        tm_data_gff3.update(parse_tmhmm(fp))

tm_data_exo = {}
for f in ['predicted_topologies_exonerate.3line']:
    fp = os.path.join(BASE_DIR, f)
    if os.path.exists(fp):
        tm_data_exo.update(parse_tmhmm(fp))

# Seleccionar TMHMM según fuente_seq del gen
tm_data_gdr = {}
fp_gdr = os.path.join(GDR_DIR, 'missing_aqp_gff3.3line')
if os.path.exists(fp_gdr): tm_data_gdr.update(parse_tmhmm(fp_gdr))

tm_data = {}
for _, r in df_tab.iterrows():
    gid = r['gene_id']
    src = r['fuente_seq']
    if src == 'MAKER_GFF3' and gid in tm_data_gdr:
        tm_data[gid] = tm_data_gdr[gid]
    elif src == 'EXONERATE' and gid in tm_data_exo:
        tm_data[gid] = tm_data_exo[gid]
    elif src == 'GFF3_FALLBACK' and gid in tm_data_gff3:
        tm_data[gid] = tm_data_gff3[gid]
    elif src == 'GFF3' and gid in tm_data_gff3:
        tm_data[gid] = tm_data_gff3[gid]
    elif gid in tm_data_gff3:
        tm_data[gid] = tm_data_gff3[gid]
    elif gid in tm_data_exo:
        tm_data[gid] = tm_data_exo[gid]
print(f"  TMHMM: {len(tm_data)} genes con topologia (5 features)")

ps_data_gff3 = {}
for f in ['pepstats_gff3.txt']:
    fp = os.path.join(BASE_DIR, f)
    if os.path.exists(fp):
        ps_data_gff3.update(parse_pepstats(fp))

ps_data_exo = {}
for f in ['pepstats_exonerate.txt']:
    fp = os.path.join(BASE_DIR, f)
    if os.path.exists(fp):
        ps_data_exo.update(parse_pepstats(fp))

# Seleccionar Pepstats según fuente_seq del gen
ps_data_gdr = {}
fp_gdr2 = os.path.join(GDR_DIR, 'missing_aqp_pepstats.pepstats')
if os.path.exists(fp_gdr2): ps_data_gdr.update(parse_pepstats(fp_gdr2))

ps_data = {}
for _, r in df_tab.iterrows():
    gid = r['gene_id']
    src = r['fuente_seq']
    if src == 'MAKER_GFF3' and gid in ps_data_gdr:
        ps_data[gid] = ps_data_gdr[gid]
    elif src == 'EXONERATE' and gid in ps_data_exo:
        ps_data[gid] = ps_data_exo[gid]
    elif src == 'GFF3_FALLBACK' and gid in ps_data_gff3:
        ps_data[gid] = ps_data_gff3[gid]
    elif src == 'GFF3' and gid in ps_data_gff3:
        ps_data[gid] = ps_data_gff3[gid]
    elif gid in ps_data_gff3:
        ps_data[gid] = ps_data_gff3[gid]
    elif gid in ps_data_exo:
        ps_data[gid] = ps_data_exo[gid]
print(f"  Pepstats: {len(ps_data)} genes con composicion (13 features)")

# ── DeepLoc: 10 probabilidades por compartimento ──
DEEPLOC_COMPS = ['Cytoplasm', 'Nucleus', 'Extracellular', 'Cell membrane',
                 'Mitochondrion', 'Plastid', 'Endoplasmic reticulum',
                 'Lysosome/Vacuole', 'Golgi apparatus', 'Peroxisome']

def load_deeploc(path):
    if not os.path.exists(path): return {}
    dl = pd.read_csv(path)
    result = {}
    for _, row in dl.iterrows():
        gid = str(row['Protein_ID']).replace('-mRNA-1', '').replace('_Benihoppe_v1', '')
        gid = exo_to_gene.get(gid, gid)
        result[gid] = {c: row.get(c, 0) for c in DEEPLOC_COMPS}
    return result

dl_gff3 = load_deeploc(os.path.join(BASE_DIR, 'deeploc_gff3.csv'))
dl_exo  = load_deeploc(os.path.join(BASE_DIR, 'deeploc_exonerate.csv'))
dl_gdr  = load_deeploc(os.path.join(GDR_DIR, 'missing_aqp_deeploc.csv'))
print(f"  DeepLoc: GFF3={len(dl_gff3)}, Exonerate={len(dl_exo)}, GDR={len(dl_gdr)} genes (10 compartimentos)")

# ─── PASO 3: UNIFICACIÓN Y BIOPYTHON ─────────────────────────────────────────
banner("PASO 3: UNIFICACIÓN FINAL Y FEATURE ENGINEERING")

df = df_tab.copy()
df = df.rename(columns={
    'gene_id': 'ID',
    'subfamilia_phylo': 'Subfamilia_Filogenetica',
    'confianza': 'Confianza',
    'localizacion': 'Localizacion'
})

# Carga de secuencias
seq_dict_gff3 = {}
path_gff3 = os.path.join(BASE_DIR, 'aquaporin_peptides.fasta')
if os.path.exists(path_gff3):
    for rec in SeqIO.parse(path_gff3, 'fasta'):
        gid = rec.id.replace('-mRNA-1', '')
        seq_dict_gff3[gid] = str(rec.seq)

seq_dict_exo = {}
path_exo = os.path.join(BASE_DIR, 'exonerate_genes_aqp.fasta')
if os.path.exists(path_exo):
    for rec in SeqIO.parse(path_exo, 'fasta'):
        seq_dict_exo[rec.id] = str(rec.seq)

seq_dict_gdr = {}
path_gdr = os.path.join(GDR_DIR, 'missing_aquaporins.fasta')
if os.path.exists(path_gdr):
    for rec in SeqIO.parse(path_gdr, 'fasta'):
        seq_dict_gdr[rec.id.replace('-mRNA-1_Benihoppe_v1', '')] = str(rec.seq)

def get_final_seq(row):
    gid = row['ID']
    src = row.get('fuente_seq', 'GFF3')
    eid = row.get('mRNA_exonerate_id')
    if src == 'MAKER_GFF3':
        return seq_dict_gdr.get(gid, "")
    elif src == 'EXONERATE' and pd.notna(eid) and str(eid) in seq_dict_exo:
        return seq_dict_exo[str(eid)]
    elif src == 'GFF3_FALLBACK':
        return seq_dict_gff3.get(gid, "")
    return seq_dict_gff3.get(gid, seq_dict_exo.get(str(eid), ""))

df['Sequence'] = df.apply(get_final_seq, axis=1)

# Recalcular pI y Mw_kDa desde la secuencia real seleccionada (no del CSV estático)
def recalc_pI_Mw(seq):
    if not seq or len(seq) < 10:
        return pd.Series({'pI': 0, 'Mw_kDa': 0})
    try:
        clean = ''.join(c for c in seq.upper() if c in 'ACDEFGHIKLMNPQRSTVWY')
        pa = ProteinAnalysis(clean)
        return pd.Series({'pI': round(pa.isoelectric_point(), 2),
                          'Mw_kDa': round(pa.molecular_weight() / 1000, 2)})
    except:
        return pd.Series({'pI': 0, 'Mw_kDa': 0})

recalc = df['Sequence'].apply(recalc_pI_Mw)
df['pI'] = recalc['pI']
df['Mw_kDa'] = recalc['Mw_kDa']
print(f"  pI y Mw_kDa: recalculados desde secuencia real (fuente_seq)")

# Biopython: 9 features fisicoquímicas (expandido)
def compute_bio(seq):
    default = {'GRAVY': 0, 'Instability': 0, 'Aromaticity': 0,
               'Helix_Frac': 0, 'Sheet_Frac': 0, 'Coil_Frac': 0,
               'Seq_Length': 0, 'Charged_Frac': 0, 'Polar_Frac': 0}
    if not seq or len(seq) < 10: return default
    try:
        clean = ''.join(c for c in seq.upper() if c in 'ACDEFGHIKLMNPQRSTVWY')
        if len(clean) < 10: return default
        pa = ProteinAnalysis(clean)
        ss = pa.secondary_structure_fraction()
        aa_pct = pa.get_amino_acids_percent()
        charged = sum(aa_pct.get(a, 0) for a in 'RKDE')
        polar = sum(aa_pct.get(a, 0) for a in 'STNQCRKDEH')
        return {
            'GRAVY': round(pa.gravy(), 4),
            'Instability': round(pa.instability_index(), 2),
            'Aromaticity': round(pa.aromaticity(), 4),
            'Helix_Frac': round(ss[0], 4),
            'Sheet_Frac': round(ss[1], 4),
            'Coil_Frac': round(ss[2], 4),
            'Seq_Length': len(clean),
            'Charged_Frac': round(charged, 4),
            'Polar_Frac': round(polar, 4)
        }
    except: return default

bio_cols_expand = ['GRAVY', 'Instability', 'Aromaticity', 'Helix_Frac', 'Sheet_Frac',
                   'Coil_Frac', 'Seq_Length', 'Charged_Frac', 'Polar_Frac']
bio_df = pd.DataFrame(df['Sequence'].apply(compute_bio).tolist())
for c in bio_cols_expand:
    if c in df.columns: df = df.drop(columns=[c])
df = pd.concat([df, bio_df], axis=1)
print(f"  Biopython: 9 features calculadas")

# Merge TMHMM (5 cols)
if tm_data:
    tm_df = pd.DataFrame.from_dict(tm_data, orient='index').reset_index().rename(columns={'index': 'ID'})
    df = df.merge(tm_df, on='ID', how='left')
    # Sobreescribir TMHs del CSV con el valor real de TMHMM (fuente_seq-aware)
    df['TMHs'] = df['TMHMM_Segments'].fillna(df['TMHs'])

# Merge Pepstats (13 cols)
if ps_data:
    ps_df = pd.DataFrame.from_dict(ps_data, orient='index').reset_index().rename(columns={'index': 'ID'})
    df = df.merge(ps_df, on='ID', how='left')

# DeepLoc (segun fuente_seq)
def get_deeploc_probs(row):
    gid = row['ID']
    src = row.get('fuente_seq', 'GFF3')
    if src == 'MAKER_GFF3':
        dl = dl_gdr.get(gid, {})
    elif src == 'EXONERATE':
        dl = dl_exo.get(gid, {})
    elif src == 'GFF3_FALLBACK':
        dl = dl_gff3.get(gid, {})
    else:
        dl = dl_gff3.get(gid, dl_exo.get(gid, {}))
    return {f'DL_{c.replace("/", "_").replace(" ", "_")}': dl.get(c, 0) for c in DEEPLOC_COMPS}

dl_feats = df.apply(get_deeploc_probs, axis=1)
dl_df_out = pd.DataFrame(dl_feats.tolist())
df = pd.concat([df.reset_index(drop=True), dl_df_out], axis=1)
print(f"  DeepLoc: {len(dl_df_out.columns)} probabilidades por compartimento")

# Merge motivos (metadatos, no entran en PCA)
df = df.merge(df_motifs, on='ID', how='left').fillna(0)

mlb = MultiLabelBinarizer()
loc_lists = df['Localizacion'].fillna('Unknown').apply(lambda x: [i.strip() for i in x.split('|')]).tolist()
loc_encoded = mlb.fit_transform(loc_lists)
loc_df = pd.DataFrame(loc_encoded, columns=[f'Loc_{l}' for l in mlb.classes_])
df = pd.concat([df, loc_df], axis=1)

# ─── PASO 4: FILTRADO Y CATEGORIZACIÓN (PARTIAL / CLEAN) ───────────────────
banner("PASO 4: FILTRADO Y CATEGORIZACIÓN ESTRUCTURAL")

# Definir Criterios de Limpieza para el Modelado
numeric_crit = ['pI', 'Mw_kDa', 'GRAVY', 'Seq_Length']
# Primero calculamos outliers sobre TODAS las candidatas razonables
# Primero calculamos outliers sobre TODAS las candidatas razonables (Train_Cat level: >= 4 TMHs)
candidatos_base = df[(df['fuente_seq'] != 'GFF3_FALLBACK') & (df['TMHs'] >= 4) & (df['Confianza'] != 'Descarte')].copy()
z_scores = np.abs(zscore(candidatos_base[numeric_crit]))
outlier_ids = candidatos_base[ (z_scores > 3).any(axis=1) ]['ID'].tolist()

# --- DOBLE CATEGORIZACIÓN ---
# Train_Cat: define qué secuencias ENTRENAN el PCA y el Random Forest (exactamente 6 TMHs)
#            Las aquaporinas completas siempre tienen 6 hélices transmembrana.
# Plot_Cat:  define cómo se COLOREAN los puntos en la gráfica (= Train_Cat)
#            Las de != 6 TMHs se muestran como PARTIAL (semitransparentes).

def get_train_cat(r):
    if (r['fuente_seq'] == 'GFF3_FALLBACK' or 
        r['fuente_seq'] == 'MAKER_GFF3' or 
        r['veredicto'] == 'AMBAS_MAL' or 
        r['Subfamilia_Filogenetica'] == 'Fragmento' or 
        r['Confianza'] == 'Descarte' or 
        r['TMHs'] != 6 or
        r['ID'] in outlier_ids):
        return 'PARTIAL'
    return r['Subfamilia_Filogenetica']

def get_plot_cat(r):
    # Unificamos todo lo "no canónico" como PARTIAL para el plot
    if (r['fuente_seq'] == 'GFF3_FALLBACK' or 
        r['fuente_seq'] == 'MAKER_GFF3' or 
        r['veredicto'] == 'AMBAS_MAL' or 
        r['Subfamilia_Filogenetica'] == 'Fragmento' or 
        r['Confianza'] == 'Descarte' or 
        r['TMHs'] != 6 or
        r['ID'] in outlier_ids):
        return 'PARTIAL'
    return r['Subfamilia_Filogenetica']

df['Train_Cat'] = df.apply(get_train_cat, axis=1)
df['Plot_Cat'] = df.apply(get_plot_cat, axis=1)

# Dataset para Visualización (genes restantes > 4 TMHMM)
df_plot = df.copy()

# Crear columna de perfil de motivos como string resumen para hover
def make_perfil_motivos(row):
    present = [m for m in motif_cols if row.get(m, 0) == 1]
    return ', '.join(present) if present else 'Ninguno'

df_plot['Perfil_Motivos'] = df_plot.apply(make_perfil_motivos, axis=1)

# Definir Orden de la Leyenda (PARTIAL al final)
counts = df_plot['Plot_Cat'].value_counts()
legend_order = [c for c in counts.index if c != 'PARTIAL'] + (['PARTIAL'] if 'PARTIAL' in counts.index else [])

# Dataset para Modelado (Solo los "Clean" segun TRAMO)
df_final = df_plot[ (df_plot['Train_Cat'] != 'PARTIAL') & (df_plot['Train_Cat'] != 'Fragmento') ].copy()
print(f"Dataset de modelado (Clean orig, para PCA): {len(df_final)} genes")
print(f"Genes marcados como PARTIAL/Excluidos: {len(df) - len(df_final)}")

# ── Calcular huella de motivos por subfamilia para la leyenda ──
SUBFAMILIAS_PLOT = ['PIP', 'TIP', 'NIP', 'SIP', 'XIP']
freq_sf = {}
for sf in SUBFAMILIAS_PLOT:
    sf_ids = df_final[df_final['Subfamilia_Filogenetica'] == sf]['ID'].tolist()
    n = len(sf_ids)
    if n == 0: continue
    freqs = {}
    for m in motif_cols:
        count = sum(1 for gid in sf_ids if df_final.loc[df_final['ID']==gid, m].values[0] == 1)
        freqs[m] = round(count / n * 100, 0)
    freq_sf[sf] = freqs

# Generar texto de huella para la leyenda
motif_legend_lines = []
for sf in SUBFAMILIAS_PLOT:
    if sf not in freq_sf: continue
    exclusive = []
    enriched = []
    absent = []
    low = []  # Motivos presentes pero raros
    for m in motif_cols:
        f_this = freq_sf[sf][m]
        f_others = [freq_sf[osf][m] for osf in SUBFAMILIAS_PLOT if osf != sf and osf in freq_sf]
        max_others = max(f_others) if f_others else 0
        avg_others = np.mean(f_others) if f_others else 0
        if f_this >= 50 and max_others < 15:
            exclusive.append(m)
        elif f_this >= 60 and avg_others < 30:
            enriched.append(m)
        # Ausente: <5% aquí y al menos otra subfamilia lo tiene >40%
        if f_this < 5 and max_others >= 40:
            absent.append(m)
        # Bajo: 5-25% aquí y al menos otra subfamilia lo tiene >50%
        elif 5 <= f_this <= 25 and max_others >= 50:
            low.append(m)
    parts = []
    if exclusive:
        parts.append(u'\u2605 ' + ','.join(exclusive))
    if enriched:
        parts.append(u'\u2191 ' + ','.join(enriched))
    if absent:
        parts.append(u'\u2717 ' + ','.join(absent))
    if low:
        parts.append(u'\u2193 ' + ','.join(low))
    motif_legend_lines.append(f"{sf}: {' | '.join(parts)}" if parts else f"{sf}: (conservados)")

# Construir cabecera solo con los símbolos que realmente aparecen
all_parts_text = '\n'.join(motif_legend_lines)
header_symbols = []
if u'\u2605' in all_parts_text: header_symbols.append(u'\u2605' + '=exclusivo')
if u'\u2191' in all_parts_text: header_symbols.append(u'\u2191' + '=enriquecido')
if u'\u2717' in all_parts_text: header_symbols.append(u'\u2717' + '=ausente')
if u'\u2193' in all_parts_text: header_symbols.append(u'\u2193' + '=bajo')
motif_box_text = 'HUELLA DE MOTIVOS\n' + '  '.join(header_symbols) + '\n' + u'\u2500' * 38 + '\n' + '\n'.join(motif_legend_lines)
print("\n" + motif_box_text)

# ─── PASO 5: MODELADO PCA Y RESULTADOS ───────────────────────────────────────
banner("PASO 5: PCA INTEGRADO Y RANKING DE IMPORTANCIA")

# ── Construcción de la matriz de features ──
physico_cols = ['pI', 'Mw_kDa', 'TMHs']
biopython_cols = ['GRAVY', 'Instability', 'Aromaticity', 'Helix_Frac',
                  'Sheet_Frac', 'Coil_Frac', 'Seq_Length', 'Charged_Frac', 'Polar_Frac']
tmhmm_cols = [c for c in df.columns if c.startswith('TMHMM_')]
pepstats_cols = [c for c in df.columns if c.startswith('PS_')]
deeploc_cols = [c for c in df.columns if c.startswith('DL_')]
loc_cols = [c for c in df.columns if c.startswith('Loc_')]

# feature_cols: TODA la física + topología + composición + localización (SIN motivos binarios)
# Los motivos M1-M12 se conservan como metadatos enriquecidos para el hover.
feature_cols = physico_cols + biopython_cols + tmhmm_cols + pepstats_cols + loc_cols
# NOTA: deeploc_cols excluidos temporalmente para probar separación sin ellos

print(f"  Features ANTES del filtro de correlación: {len(feature_cols)}")
print(f"    Fisicoquímicas básicas:  {len(physico_cols)}")
print(f"    Biopython:               {len(biopython_cols)}")
print(f"    TMHMM topología:         {len(tmhmm_cols)}")
print(f"    Pepstats composición:    {len(pepstats_cols)}")
print(f"    DeepLoc compartimentos:  {len(deeploc_cols)} (excluidos del PCA)")
print(f"    Localización binarizada: {len(loc_cols)}")

# ── FILTRO DE MULTICOLINEALIDAD (|r| > 0.85) ──────────────────────────────────
# Calcula la matriz de correlación absoluta de Pearson sobre el dataset de
# entrenamiento (df_final) y elimina una variable de cada par altamente
# correlacionado. Se conserva la variable con menor correlación media respecto
# al resto de features para maximizar la información retenida.
CORR_THRESHOLD = 0.85

corr_matrix = df_final[feature_cols].fillna(0).corr(method='pearson').abs()

# Máscara triangular superior (sin diagonal) para no duplicar pares
upper_tri = corr_matrix.where(
    np.triu(np.ones(corr_matrix.shape, dtype=bool), k=1)
)

# Identificar pares con |r| > umbral
high_corr_pairs = []
for col in upper_tri.columns:
    correlated = upper_tri.index[upper_tri[col] > CORR_THRESHOLD].tolist()
    for row_name in correlated:
        high_corr_pairs.append((row_name, col, corr_matrix.loc[row_name, col]))

# Para cada par, eliminar la variable con mayor correlación media al resto
to_drop = set()
for v1, v2, r_val in high_corr_pairs:
    if v1 in to_drop or v2 in to_drop:
        continue  # Ya se marcó una del par
    mean_corr_v1 = corr_matrix[v1].drop(v1).mean()
    mean_corr_v2 = corr_matrix[v2].drop(v2).mean()
    drop = v1 if mean_corr_v1 >= mean_corr_v2 else v2
    keep = v2 if drop == v1 else v1
    to_drop.add(drop)
    print(f"    [CORR] {v1} ↔ {v2}  (r={r_val:.3f})  → Descartada: {drop}  (conservada: {keep})")

# Actualizar feature_cols eliminando las variables redundantes
feature_cols_original = list(feature_cols)
feature_cols = [c for c in feature_cols if c not in to_drop]

print(f"\n  ── Resumen del filtro de correlación (umbral |r| > {CORR_THRESHOLD}) ──")
print(f"    Variables descartadas: {len(to_drop)} de {len(feature_cols_original)}")
if to_drop:
    print(f"    Lista descartada: {sorted(to_drop)}")
print(f"    Features finales para PCA/RF: {len(feature_cols)}")
# ──────────────────────────────────────────────────────────────────────────────

X_fit = df_final[feature_cols].fillna(0).values
scaler = StandardScaler()
X_fit_scaled = scaler.fit_transform(X_fit)
X_all_scaled = scaler.transform(df_plot[feature_cols].fillna(0).values)

pca = PCA(n_components=2)
X_pca_fit = pca.fit_transform(X_fit_scaled)
X_pca_all = pca.transform(X_all_scaled)

# Exportar Coordenadas Finales y Metadatos para análisis externo
df_pca = pd.DataFrame(X_pca_all, columns=['PC1', 'PC2'])
df_pca['ID'] = df_plot['ID'].values
df_pca['Subfamilia'] = df_plot['Subfamilia_Filogenetica'].values
df_pca['Total_Motivos'] = df_plot['Total_Motivos'].values
df_pca['Score_Diagnostico'] = df_plot['Score_Diagnostico'].values
df_pca.to_csv(os.path.join(OUT, 'PCA_Coordenadas_Finales.csv'), index=False)

rf = RandomForestClassifier(n_estimators=500, random_state=42)
rf.fit(X_fit_scaled, df_final['Subfamilia_Filogenetica'])
feat_imp = pd.DataFrame({'Feature': feature_cols, 'Importance': rf.feature_importances_}).sort_values('Importance', ascending=False)
feat_imp.to_csv(os.path.join(OUT, 'RANKING_FINAL_INTEGRADO.csv'), index=False)

# ─── PREDICCIÓN DE MAKER_GFF3 CON EL MODELO INTEGRADO ───
mask_maker = df_plot['fuente_seq'] == 'MAKER_GFF3'
if mask_maker.any():
    print(f"\n  [PREDICCIÓN] Reclasificando {mask_maker.sum()} secuencias MAKER_GFF3 con el RF limpio...")
    X_maker = df_plot.loc[mask_maker, feature_cols].fillna(0).values
    X_maker_scaled = scaler.transform(X_maker)
    preds = rf.predict(X_maker_scaled)
    
    # 1. Actualizar memoria para las gráficas que se generarán después
    df_plot.loc[mask_maker, 'Subfamilia_Filogenetica'] = preds
    
    # 2. Actualizar copia en el sistema de archivos
    contador = {}
    for i, idx in enumerate(df_plot[mask_maker].index):
        gid = df_plot.loc[idx, 'ID']
        pred_sf = preds[i]
        contador[pred_sf] = contador.get(pred_sf, 0) + 1
        
        # En tabular (ya cargado como df_tab al principio del script)
        tab_mask = df_tab['gene_id'] == gid
        if tab_mask.any():
            df_tab.loc[tab_mask, 'aqp_family_subfamily'] = pred_sf
            df_tab.loc[tab_mask, 'subfamilia_phylo'] = pred_sf

    print(f"    Distribución de subfamilias predichas: {contador}")
    df_tab.to_csv(os.path.join(BASE_DIR, 'tabla_aquaporinas_traduccion.tabular'), sep='\t', index=False)
    print("    → Archivo TABULAR actualizado y sincronizado en disco.")
# ────────────────────────────────────────────────────────

# PLOTS
plt.figure(figsize=(10, 8))
sns.barplot(data=feat_imp.head(25), x='Importance', y='Feature', palette='viridis')
# Título en el pie de figura (norma APA/UCAM), no embebido en la imagen
plt.tight_layout()
plt.savefig(os.path.join(OUT, 'FINAL_IMPORTANCE.png'))

varianza_total = (pca.explained_variance_ratio_[0] + pca.explained_variance_ratio_[1]) * 100

def plot_static_pca(include_partial=True, filename='PCA_FINAL_INTEGRADO.png'):
    fig_st, ax_st = plt.subplots(figsize=(12, 10))
    import matplotlib.colors as mcolors

    # Extraer número de sub-subfamilia (ej: FaPIP2 → '2', FaTIP1 → '1')
    def get_subsf_num(aqp_detail):
        if pd.isna(aqp_detail) or not aqp_detail:
            return ''
        m = re.search(r'(\d+)$', str(aqp_detail))
        return m.group(1) if m else ''

    df_plot['SubSF_Num'] = df_plot['aqp_family_subfamily'].apply(get_subsf_num)

    # Primero plotear las subfamilias "clean" (no PARTIAL)
    for cat in legend_order:
        if cat == 'PARTIAL':
            continue
        mask = df_plot['Plot_Cat'] == cat
        ax_st.scatter(X_pca_all[mask, 0], X_pca_all[mask, 1], c=COLOR_MAP.get(cat, 'gray'), 
                   label=f"{cat} (n={sum(mask)})", s=120, alpha=0.75, edgecolors='black', linewidth=0.5, zorder=2)
        # Números de sub-subfamilia dentro de cada punto
        for idx in df_plot.index[mask]:
            pos = df_plot.index.get_loc(idx)
            num = df_plot.loc[idx, 'SubSF_Num']
            if num:
                ax_st.text(X_pca_all[pos, 0], X_pca_all[pos, 1], num,
                          ha='center', va='center', fontsize=6, fontweight='bold',
                          color='white', zorder=3)

    # PARTIAL: color subfamilia transparente + borde gris fino
    if include_partial and 'PARTIAL' in legend_order:
        partial_mask = df_plot['Plot_Cat'] == 'PARTIAL'
        partial_subfams = df_plot.loc[partial_mask, 'Subfamilia_Filogenetica'].values
        partial_colors = [COLOR_MAP.get(s, '#808080') for s in partial_subfams]
        partial_rgba = [mcolors.to_rgba(c, alpha=0.25) for c in partial_colors]
        ax_st.scatter(X_pca_all[partial_mask, 0], X_pca_all[partial_mask, 1], c=partial_rgba,
                   s=120, edgecolors='#666666', linewidth=1, zorder=2, label=f"PARTIAL (n={sum(partial_mask)})")
        # Números también para PARTIAL
        for idx in df_plot.index[partial_mask]:
            pos = df_plot.index.get_loc(idx)
            num = df_plot.loc[idx, 'SubSF_Num']
            if num:
                ax_st.text(X_pca_all[pos, 0], X_pca_all[pos, 1], num,
                          ha='center', va='center', fontsize=6, fontweight='bold',
                          color='#999999', zorder=3)

    def draw_correct_ellipse(ax, x, y, color):
        if len(x) < 5: return
        cov = np.cov(x, y)
        vals, vecs = np.linalg.eigh(cov)
        order = vals.argsort()[::-1]
        vals, vecs = vals[order], vecs[:, order]
        angle = np.degrees(np.arctan2(vecs[1,0], vecs[0,0]))
        width, height = 2 * 2.447 * np.sqrt(np.maximum(vals, 0))
        ellipse = Ellipse(xy=(np.mean(x), np.mean(y)), width=width, height=height, 
                         angle=angle, facecolor=color, alpha=0.1, edgecolor=color, 
                         linestyle='--', linewidth=2, zorder=1)
        ax.add_patch(ellipse)

    for subfam in df_final['Subfamilia_Filogenetica'].unique():
        mask = df_final['Subfamilia_Filogenetica'] == subfam
        draw_correct_ellipse(ax_st, X_pca_fit[mask, 0], X_pca_fit[mask, 1], COLOR_MAP[subfam])

    # AÑADIR VECTORES (LOADINGS) DE LAS VARIABLES MÁS IMPORTANTES
    # Biplot clásico: las flechas parten desde el origen (0, 0) del espacio PCA.
    top_vars = feat_imp.head(6)['Feature'].tolist()
    scale_arrows = (np.max(np.abs(X_pca_all[:, :2])) * 0.6) / np.max(np.abs(pca.components_[:2, :]))

    for i, var in enumerate(feature_cols):
        if var in top_vars:
            arrow_x = pca.components_[0, i] * scale_arrows
            arrow_y = pca.components_[1, i] * scale_arrows
            # Z-order 1 las pone por debajo de los puntos (que tienen zorder 2)
            ax_st.arrow(0, 0, arrow_x, arrow_y, color='#B0B0B0', alpha=0.9,
                        head_width=0.15, head_length=0.2, linewidth=1.5, zorder=1)
            clean_name = var.replace('Loc_', '').replace('PS_', '')
            ax_st.text(arrow_x * 1.15, arrow_y * 1.15, clean_name, color='#707070',
                       ha='center', va='center', fontsize=10, fontweight='bold',
                       bbox=dict(facecolor='white', alpha=0.4, edgecolor='none', pad=0.5), zorder=1)

    ax_st.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    ax_st.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    varianza_total = (pca.explained_variance_ratio_[0] + pca.explained_variance_ratio_[1]) * 100
    # Título en el pie de figura (norma APA/UCAM); la varianza explicada se indica en el pie
    ax_st.legend(fontsize=9, framealpha=0.9)
    ax_st.grid(alpha=0.2)
    ax_st.margins(0.1)

    # Recuadro de huella de motivos retirado: era metadato MEME (M1-M12 del análisis
    # exploratorio) que no aporta a la figura del PCA y resultaba inconsistente con la
    # numeración M1-M15 del análisis de motivos definitivo (Figura 5).

    plt.tight_layout()
    fig_st.savefig(os.path.join(OUT, filename), bbox_inches='tight')
    plt.close(fig_st)

print("  Generando gráficos PCA estáticos...")
plot_static_pca(include_partial=True, filename='PCA_FINAL_INTEGRADO.png')
plot_static_pca(include_partial=False, filename='PCA_FINAL_INTEGRADO_CLEAN.png')

import plotly.express as px
import plotly.graph_objects as go


# 4. Plotly Interactivo (HTML)
print("[PLOT 3/3] Generando explorador interactivo (HTML)...")
df_pca_html = pd.DataFrame(X_pca_all, columns=['PC1', 'PC2'])

# Identificadores estéticos (igual que en multidimensional)
def make_display_id(row):
    gid = row['ID']
    eid = row.get('mRNA_exonerate_id')
    src = row.get('fuente_seq', 'GFF3')
    ver = row.get('veredicto', '')
    # Solo marcar partial_manual si es GFF3_FALLBACK
    if src == 'GFF3_FALLBACK' and ver == 'MANUAL_CURATED': return f"{gid}-partial_manual"
    if src == 'EXONERATE' and pd.notna(eid): return f"{eid}-{gid}"
    if src == 'GFF3_FALLBACK': return f"{gid}-partial"
    return gid

df_pca_html['Display_ID'] = df_plot.apply(make_display_id, axis=1).values
df_pca_html['Cat'] = df_plot['Plot_Cat'].values
df_pca_html['Subfamilia'] = df_plot['Subfamilia_Filogenetica'].values
df_pca_html['AQP_Detail'] = df_plot['aqp_family_subfamily'].values
df_pca_html['Localizacion'] = df_plot['Localizacion'].values
df_pca_html['pI'] = df_plot['pI'].round(2).values
df_pca_html['Mw'] = df_plot['Mw_kDa'].round(2).values
df_pca_html['Total_Motivos'] = df_plot['Total_Motivos'].values
df_pca_html['Score_Diag'] = df_plot['Score_Diagnostico'].values
df_pca_html['Perfil_Motivos'] = df_plot['Perfil_Motivos'].values

# Top features para hover (Añadirlas al DF de Plotly)
top_features = feat_imp.head(5)['Feature'].tolist()
for f in top_features:
    df_pca_html[f] = df_plot[f].values

hover_cols = ['AQP_Detail', 'Localizacion', 'pI', 'Mw', 'Total_Motivos', 'Score_Diag', 'Perfil_Motivos'] + top_features

# Para Plotly: crear una columna de color real (subfamilia) para los PARTIAL
df_pca_html['Real_Subfam'] = df_plot['Subfamilia_Filogenetica'].values

# Separar PARTIAL y no-PARTIAL para plotearlos con estilos distintos
df_pca_html['SubSF_Num'] = df_plot['SubSF_Num'].values
df_clean = df_pca_html[df_pca_html['Cat'] != 'PARTIAL'].copy()
df_partial = df_pca_html[df_pca_html['Cat'] == 'PARTIAL'].copy()

fig_int = go.Figure()

# Plotear primero las subfamilias clean
for cat in legend_order:
    if cat == 'PARTIAL': continue
    sub = df_clean[df_clean['Cat'] == cat]
    if len(sub) == 0: continue
    fig_int.add_trace(go.Scatter(
        x=sub['PC1'], y=sub['PC2'], mode='markers+text',
        name=f"{cat} (n={len(sub)})",
        text=sub['SubSF_Num'],
        textfont=dict(size=8, color='white', family='Arial Black'),
        textposition='middle center',
        customdata=sub[hover_cols].values,
        hovertext=sub['Display_ID'],
        hovertemplate='<b>%{hovertext}</b><br>' + '<br>'.join(f'{c}: %{{customdata[{i}]}}' for i, c in enumerate(hover_cols)) + '<extra></extra>',
        marker=dict(size=14, color=COLOR_MAP.get(cat, 'gray'), opacity=0.8,
                    line=dict(width=1, color='black'))
    ))

# Plotear PARTIAL: color de su subfamilia real pero muy transparente + borde negro grueso
if len(df_partial) > 0:
    # Un unico trace para la leyenda, pero con colores individuales por subfamilia
    partial_colors = [COLOR_MAP.get(s, '#808080') for s in df_partial['Real_Subfam']]
    fig_int.add_trace(go.Scatter(
        x=df_partial['PC1'], y=df_partial['PC2'], mode='markers+text',
        name=f"PARTIAL (n={len(df_partial)})",
        text=df_partial['SubSF_Num'],
        textfont=dict(size=8, color='#999999', family='Arial Black'),
        textposition='middle center',
        customdata=df_partial[hover_cols].values,
        hovertext=df_partial['Display_ID'],
        hovertemplate='<b>%{hovertext}</b><br>' + '<br>'.join(f'{c}: %{{customdata[{i}]}}' for i, c in enumerate(hover_cols)) + '<extra></extra>',
        marker=dict(size=14, color=partial_colors, opacity=0.25,
                    line=dict(width=1.5, color='#666666'))
    ))

fig_int.update_layout(
    template='plotly_white',
    title=f'<b>Explorador Multidimensional: {varianza_total:.1f}% Varianza Explicada</b>',
)

# AÑADIR ELIPSES INTERACTIVAS
def get_plotly_ellipse(x, y, color, name):
    if len(x) < 5: return None
    cov = np.cov(x, y)
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    
    # Parametric ellipse
    theta = np.linspace(0, 2*np.pi, 100)
    # 2.447 for 95% confidence
    a, b = 2.447 * np.sqrt(np.maximum(vals[0],0)), 2.447 * np.sqrt(np.maximum(vals[1],0))
    
    # Rotation matrix
    rotation_matrix = vecs
    ellipse_base = np.array([a * np.cos(theta), b * np.sin(theta)])
    ellipse_final = np.dot(rotation_matrix, ellipse_base)
    
    return go.Scatter(
        x=ellipse_final[0, :] + np.mean(x),
        y=ellipse_final[1, :] + np.mean(y),
        mode='lines',
        line=dict(color=color, width=2, dash='dash'),
        fill='toself',
        fillcolor=color,
        opacity=0.1,
        name=f"95% CI {name}",
        showlegend=False,
        hoverinfo='skip'
    )

for subfam in df_final['Subfamilia_Filogenetica'].unique():
    mask = df_final['Subfamilia_Filogenetica'] == subfam
    e_trace = get_plotly_ellipse(X_pca_fit[mask, 0], X_pca_fit[mask, 1], COLOR_MAP.get(subfam, 'gray'), subfam)
    if e_trace: fig_int.add_trace(e_trace)

# AÑADIR VECTORES (LOADINGS) AL HTML (biplot clásico: desde el origen)
scale_arrows_plotly = (np.max(np.abs(X_pca_all[:, :2])) * 0.6) / np.max(np.abs(pca.components_[:2, :]))

for i, var in enumerate(feature_cols):
    if var in top_features[:6]:
        arrow_x = pca.components_[0, i] * scale_arrows_plotly
        arrow_y = pca.components_[1, i] * scale_arrows_plotly
        clean_name = var.replace('Loc_', '').replace('PS_', '')

        # Text label
        fig_int.add_trace(go.Scatter(
            x=[arrow_x * 1.1], y=[arrow_y * 1.1],
            mode='text', text=[f"<b>{clean_name}</b>"],
            textposition='middle center',
            textfont=dict(color='#888888', size=12),
            showlegend=False, hoverinfo='skip'
        ))

        # Arrow annotation (desde el origen)
        fig_int.add_annotation(
            x=arrow_x, y=arrow_y,
            ax=0, ay=0,
            xref="x", yref="y", axref="x", ayref="y",
            showarrow=True, arrowhead=2, arrowsize=1, arrowwidth=2, arrowcolor="#B0B0B0", opacity=0.8
        )

# Recuadro de huella de motivos retirado tambien del HTML interactivo:
# los motivos MEME no son features del PCA (son metadato), igual que se
# retiro de la figura estatica. El explorador se centra en las coordenadas.

fig_int.update_layout(
    font_family="DejaVu Sans", title_font_size=20,
    legend_title_text='Categor\u00eda',
    hoverlabel=dict(bgcolor="white", font_size=13)
)

fig_int.write_html(os.path.join(OUT, 'PCA_INTERACTIVO_FINAL.html'))
print(f"    → PCA_INTERACTIVO_FINAL.html")

banner("PIPELINE FINALIZADO CON ÉXITO")
print(f"Resultados en: {OUT}")
