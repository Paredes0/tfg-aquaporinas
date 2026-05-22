#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
    ANÁLISIS DE MOTIVOS MEME POR SUBFAMILIA (v2.0 - MEME Unificado)
    - Usa MEME_exonerate_gff3_aqp.txt (FASTA combinado, 258 secuencias)
    - Selecciona la secuencia correcta según fuente_seq del tabular
    - Frecuencia, exclusividad, ausencias y visualizaciones
===============================================================================
"""

import warnings
warnings.filterwarnings('ignore')

import sys, io, os, re
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

# ─── Configuración ───────────────────────────────────────────────────────────
plt.rcParams['figure.dpi'] = 200
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'DejaVu Sans'

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))
from scripts.common import config as _cfg

BASE = str(_cfg.CURADO_DIR)
OUT  = str(_cfg.RESULTS_DIR / 'analisis_motivos_unificado')
os.makedirs(OUT, exist_ok=True)

COLOR_MAP = {
    'PIP': '#E74C3C', 'TIP': '#3498DB', 'NIP': '#2ECC71',
    'SIP': '#F39C12', 'XIP': '#9B59B6'
}
SUBFAMILIAS = ['PIP', 'TIP', 'NIP', 'SIP', 'XIP']

def banner(title):
    print("\n" + "=" * 80 + "\n" + f" {title} ".center(80, " ") + "\n" + "=" * 80)

# =============================================================================
# PASO 1: CARGA DE DATOS
# =============================================================================
banner("PASO 1: CARGA DE DATOS Y MAPEO")

df_tab = pd.read_csv(os.path.join(BASE, 'tabla_aquaporinas_traduccion.tabular'), sep='\t')
df_class = pd.read_csv(os.path.join(BASE, 'clasificacion_filogenetica_simple.csv'))

# Merge para tener subfamilia + fuente_seq
df = df_class.merge(
    df_tab[['gene_id', 'fuente_seq', 'veredicto', 'mRNA_exonerate_id', 'subfamilia_phylo', 'aqp_family_subfamily']],
    left_on='ID', right_on='gene_id', how='left'
)
df['Subfamilia_Filogenetica'] = df['subfamilia_phylo'].fillna(df['Subfamilia_Filogenetica'])

# Filtrar solo genes "clean" (no fragmentos, no descartes, TMHs >= 4)
df_clean = df[
    (df['fuente_seq'] != 'GFF3_FALLBACK') &
    (df['veredicto'] != 'MANUAL_CURATED') &
    (df['veredicto'] != 'AMBAS_MAL') &
    (df['Subfamilia_Filogenetica'] != 'Fragmento') &
    (df['Confianza'] != 'Descarte') &
    (df['TMHs'] >= 4)
].copy()

print(f"[INFO] Dataset limpio: {len(df_clean)} aquaporinas válidas")
for sf in SUBFAMILIAS:
    n = len(df_clean[df_clean['Subfamilia_Filogenetica'] == sf])
    print(f"  {sf}: {n}")

# Mapeo: para cada gene_id, su ID en el FASTA/MEME
selected_fasta_id = {}
for _, r in df_tab.iterrows():
    gid = r['gene_id']
    src = r['fuente_seq']
    eid = str(r['mRNA_exonerate_id']) if pd.notna(r.get('mRNA_exonerate_id')) else None
    if src == 'EXONERATE' and eid:
        selected_fasta_id[gid] = eid
    else:
        selected_fasta_id[gid] = gid

# =============================================================================
# PASO 2: PARSEO DEL MEME UNIFICADO
# =============================================================================
banner("PASO 2: PARSEO DE MEME UNIFICADO")

MEME_FILE = os.path.join(BASE, 'MEME_exonerate_gff3_aqp.txt')

# Extraer nombres de motivos del MEME
motif_names = {}
with open(MEME_FILE, 'r') as f:
    for line in f:
        m = re.match(r'MOTIF (\S+) MEME-(\d+)', line)
        if m:
            motif_names[f"M{m.group(2)}"] = m.group(1)

print("[INFO] Motivos encontrados en MEME:")
for mid, mname in sorted(motif_names.items(), key=lambda x: int(x[0][1:])):
    print(f"  {mid}: {mname}")

# Parsear Combined block diagrams
def parse_meme_combined_block(filepath):
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
                if len(parts) >= 3 and not line.startswith('['):
                    sid = parts[0].replace('-mRNA-1', '').replace('-partial', '')
                    current_sid = sid
                    diagram_part = parts[2]
                else:
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

all_motifs_data = parse_meme_combined_block(MEME_FILE)
print(f"\n[INFO] Total secuencias parseadas: {len(all_motifs_data)}")

# Asignar motivos a cada gen según fuente_seq
motif_cols = [f"M{i}" for i in range(1, 13)]
gene_motifs = {}  # gene_id -> set of motifs

for _, row in df_clean.iterrows():
    gid = row['ID']
    fasta_id = selected_fasta_id.get(gid, gid)
    if fasta_id in all_motifs_data:
        gene_motifs[gid] = all_motifs_data[fasta_id]
    else:
        gene_motifs[gid] = set()

mapped = sum(1 for v in gene_motifs.values() if len(v) > 0)
print(f"[INFO] Genes con motivos asignados: {mapped}/{len(gene_motifs)}")

# =============================================================================
# PASO 3: MATRIZ DE FRECUENCIA POR SUBFAMILIA
# =============================================================================
banner("PASO 3: MATRIZ DE FRECUENCIA POR SUBFAMILIA")

# Calcular frecuencia de cada motivo en cada subfamilia
freq_matrix = {}
for sf in SUBFAMILIAS:
    sf_genes = df_clean[df_clean['Subfamilia_Filogenetica'] == sf]['ID'].tolist()
    n_total = len(sf_genes)
    if n_total == 0:
        continue
    freqs = {}
    for m in motif_cols:
        count = sum(1 for gid in sf_genes if m in gene_motifs.get(gid, set()))
        freqs[m] = round(count / n_total * 100, 1)
    freq_matrix[sf] = freqs

# Imprimir tabla
header = "Subfamilia   " + "  ".join(f"{m:>4}" for m in motif_cols)
print(header)
print("-" * len(header))
for sf in SUBFAMILIAS:
    if sf in freq_matrix:
        vals = "  ".join(f"{freq_matrix[sf][m]:4.0f}%" for m in motif_cols)
        n = len(df_clean[df_clean['Subfamilia_Filogenetica'] == sf])
        print(f"{sf:<13}{vals}  (n={n})")

# Guardar como CSV
freq_df = pd.DataFrame(freq_matrix).T
freq_df.index.name = 'Subfamilia'
freq_df.to_csv(os.path.join(OUT, 'FRECUENCIA_Motivos_por_Subfamilia.csv'))
print(f"\n[SAVED] FRECUENCIA_Motivos_por_Subfamilia.csv")

# =============================================================================
# PASO 4: ANÁLISIS DE EXCLUSIVIDAD Y AUSENCIAS
# =============================================================================
banner("PASO 4: ANÁLISIS DE EXCLUSIVIDAD Y AUSENCIAS")

exclusivity_report = []
absence_report = []

for sf in SUBFAMILIAS:
    if sf not in freq_matrix:
        continue
    
    n_prots = len(df_clean[df_clean['Subfamilia_Filogenetica'] == sf])
    print(f"\n[{sf}] (n={n_prots} proteínas)")
    
    universales = []
    exclusivos = []
    ausentes = []
    enriquecidos = []
    
    for m in motif_cols:
        freq_this = freq_matrix[sf][m]
        freq_others = [freq_matrix[osf][m] for osf in SUBFAMILIAS if osf != sf and osf in freq_matrix]
        avg_others = np.mean(freq_others) if freq_others else 0
        max_others = max(freq_others) if freq_others else 0
        
        # Universales: > 90% en esta subfamilia
        if freq_this >= 90:
            universales.append((m, freq_this))
        
        # Exclusivos: > 50% aquí y < 15% en TODAS las demás
        if freq_this >= 50 and max_others < 15:
            exclusivos.append((m, freq_this, avg_others))
            exclusivity_report.append({
                'Subfamilia': sf, 'Motivo': m,
                'Nombre_Motivo': motif_names.get(m, '?'),
                'Freq_%': freq_this, 'Avg_Otras_%': round(avg_others, 1),
                'Tipo': 'EXCLUSIVO'
            })
        
        # Enriquecidos: > 60% aquí y < 30% promedio en otras
        elif freq_this >= 60 and avg_others < 30:
            enriquecidos.append((m, freq_this, avg_others))
            exclusivity_report.append({
                'Subfamilia': sf, 'Motivo': m,
                'Nombre_Motivo': motif_names.get(m, '?'),
                'Freq_%': freq_this, 'Avg_Otras_%': round(avg_others, 1),
                'Tipo': 'ENRIQUECIDO'
            })
        
        # Ausentes: < 10% aquí pero > 50% promedio en otras
        if freq_this < 10 and avg_others >= 50:
            ausentes.append((m, freq_this, avg_others))
            absence_report.append({
                'Subfamilia_sin_motivo': sf, 'Motivo': m,
                'Nombre_Motivo': motif_names.get(m, '?'),
                'Freq_en_Subfamilia_%': freq_this,
                'Avg_en_Otras_%': round(avg_others, 1)
            })
    
    # Imprimir resultados
    if universales:
        print(f"  > Motivos Universales (>90%): {', '.join(m for m, _ in universales)}")
    else:
        print(f"  > Motivos Universales: Ninguno")
    
    if exclusivos:
        print(f"  > Motivos EXCLUSIVOS:")
        for m, freq, avg_o in exclusivos:
            print(f"    ★ {m} ({motif_names.get(m, '?')[:40]}): {freq:.0f}% aquí vs {avg_o:.0f}% otras")
    
    if enriquecidos:
        print(f"  > Motivos Enriquecidos:")
        for m, freq, avg_o in enriquecidos:
            print(f"    ↑ {m} ({motif_names.get(m, '?')[:40]}): {freq:.0f}% aquí vs {avg_o:.0f}% otras")
    
    if ausentes:
        print(f"  > Motivos AUSENTES (<10% aquí, >50% en otras):")
        for m, freq, avg_o in ausentes:
            print(f"    ✗ {m} ({motif_names.get(m, '?')[:40]}): {freq:.0f}% aquí vs {avg_o:.0f}% otras")

# Guardar reportes
if exclusivity_report:
    pd.DataFrame(exclusivity_report).to_csv(os.path.join(OUT, 'MOTIVOS_Exclusivos_Enriquecidos.csv'), index=False)
    print(f"\n[SAVED] MOTIVOS_Exclusivos_Enriquecidos.csv")

if absence_report:
    pd.DataFrame(absence_report).to_csv(os.path.join(OUT, 'MOTIVOS_Ausentes_por_Subfamilia.csv'), index=False)
    print(f"[SAVED] MOTIVOS_Ausentes_por_Subfamilia.csv")

# =============================================================================
# PASO 5: INFORME TEXTUAL COMPLETO
# =============================================================================
banner("PASO 5: GENERANDO INFORME TEXTUAL")

informe_lines = []
informe_lines.append("=" * 80)
informe_lines.append("INFORME ESTADÍSTICO: HUELLAS DACTILARES DE MOTIVOS MEME (v2.0)")
informe_lines.append("Fuente: MEME_exonerate_gff3_aqp.txt (FASTA combinado unificado)")
informe_lines.append("=" * 80)
informe_lines.append("")
informe_lines.append("Este informe analiza la distribución de los 12 motivos MEME identificados en las")
informe_lines.append("acuaporinas de Fragaria x ananassa, usando el análisis MEME unificado.")
informe_lines.append("")
informe_lines.append("-" * 60)
informe_lines.append("1. CORRESPONDENCIA DE MOTIVOS")
informe_lines.append("-" * 60)
for mid in motif_cols:
    informe_lines.append(f"  {mid}: {motif_names.get(mid, '?')}")

informe_lines.append("")
informe_lines.append("-" * 60)
informe_lines.append("2. MATRIZ DE FRECUENCIA POR SUBFAMILIA (%)")
informe_lines.append("-" * 60)
header_inf = f"{'Subfamilia':<12}" + " ".join(f"{m:>5}" for m in motif_cols)
informe_lines.append(header_inf)
informe_lines.append("-" * len(header_inf))
for sf in SUBFAMILIAS:
    if sf in freq_matrix:
        vals = " ".join(f"{freq_matrix[sf][m]:4.0f}%" for m in motif_cols)
        informe_lines.append(f"{sf:<12}{vals}")

informe_lines.append("")
informe_lines.append("-" * 60)
informe_lines.append("3. HALLAZGOS POR SUBFAMILIA (DIAGNÓSTICOS)")
informe_lines.append("-" * 60)

for sf in SUBFAMILIAS:
    if sf not in freq_matrix:
        continue
    n_prots = len(df_clean[df_clean['Subfamilia_Filogenetica'] == sf])
    informe_lines.append(f"\n[{sf}] (n={n_prots})")
    
    uni = [m for m in motif_cols if freq_matrix[sf][m] >= 90]
    exc = [m for m in motif_cols 
           if freq_matrix[sf][m] >= 50 and 
           max(freq_matrix[osf][m] for osf in SUBFAMILIAS if osf != sf and osf in freq_matrix) < 15]
    aus = [m for m in motif_cols 
           if freq_matrix[sf][m] < 10 and 
           np.mean([freq_matrix[osf][m] for osf in SUBFAMILIAS if osf != sf and osf in freq_matrix]) >= 50]
    
    informe_lines.append(f"  > Motivos Universales: {', '.join(uni) if uni else 'Ninguno'}")
    informe_lines.append(f"  > Motivos Ausentes:    {', '.join(aus) if aus else 'Ninguno'}")
    informe_lines.append(f"  > Motivos Exclusivos:  {', '.join(exc) if exc else 'Ninguno'}")

informe_lines.append("")
informe_lines.append("-" * 60)
informe_lines.append("4. ANÁLISIS DE EXCLUSIVIDAD (ASOCIACIÓN FUERTE)")
informe_lines.append("-" * 60)
for item in exclusivity_report:
    if item['Tipo'] == 'EXCLUSIVO':
        informe_lines.append(f"  ★ {item['Motivo']} ({item['Nombre_Motivo'][:40]}): EXCLUSIVO de {item['Subfamilia']}")

informe_lines.append("")
informe_lines.append("=" * 80)
informe_lines.append("CONCLUSIÓN: Este análisis utiliza un MEME unificado sobre el FASTA combinado")
informe_lines.append("(GFF3 + Exonerate), garantizando que los 12 motivos sean directamente")
informe_lines.append("comparables entre todas las secuencias.")
informe_lines.append("=" * 80)

informe_path = os.path.join(OUT, 'INFORME_HUELLA_MOTIVOS_v2.txt')
with open(informe_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(informe_lines))
print(f"[SAVED] INFORME_HUELLA_MOTIVOS_v2.txt")

# =============================================================================
# PASO 6: VISUALIZACIONES
# =============================================================================
banner("PASO 6: VISUALIZACIONES")

# ── HEATMAP: Frecuencia de motivos por subfamilia ──
print("[PLOT 1/3] Heatmap de frecuencia de motivos...")

heatmap_data = pd.DataFrame(freq_matrix).T[motif_cols]
# Añadir nombres de motivos como segundo nivel de columna
col_labels = [f"{m}\n{motif_names.get(m, '?')[:15]}" for m in motif_cols]

fig, ax = plt.subplots(figsize=(16, 6))
sns.heatmap(
    heatmap_data.values, annot=True, fmt='.0f', cmap='RdYlGn',
    xticklabels=col_labels, yticklabels=heatmap_data.index,
    vmin=0, vmax=100, linewidths=1, linecolor='white',
    cbar_kws={'label': 'Frecuencia (%)'},
    annot_kws={'size': 11, 'fontweight': 'bold'},
    ax=ax
)
ax.set_title('FRECUENCIA DE MOTIVOS MEME POR SUBFAMILIA (%)\n(MEME unificado: GFF3 + Exonerate)',
             fontsize=14, fontweight='bold', pad=15)
ax.set_ylabel('')
ax.set_xlabel('')
plt.xticks(rotation=45, ha='right', fontsize=9)
plt.yticks(fontsize=12)
plt.tight_layout()
fig.savefig(os.path.join(OUT, 'HEATMAP_Motivos_Subfamilia.png'), bbox_inches='tight')
print("  -> HEATMAP_Motivos_Subfamilia.png")
plt.close()

# ── HEATMAP 2: Presencia/Ausencia binario con color por subfamilia ──
print("[PLOT 2/3] Heatmap binario presencia/ausencia...")

fig2, ax2 = plt.subplots(figsize=(16, 6))

# Crear una versión categorizada: 0=ausente(<10%), 1=bajo(10-50%), 2=medio(50-80%), 3=alto(>80%)
cat_data = heatmap_data.copy()
for c in cat_data.columns:
    cat_data[c] = pd.cut(cat_data[c], bins=[-1, 10, 50, 80, 101], labels=[0, 1, 2, 3]).astype(int)

cmap_custom = matplotlib.colors.ListedColormap(['#2C3E50', '#E67E22', '#F1C40F', '#27AE60'])

sns.heatmap(
    cat_data.values, annot=heatmap_data.values.astype(int).astype(str),
    fmt='s', cmap=cmap_custom,
    xticklabels=col_labels, yticklabels=cat_data.index,
    linewidths=2, linecolor='white',
    cbar=False, annot_kws={'size': 12, 'fontweight': 'bold', 'color': 'white'},
    ax=ax2
)

# Leyenda manual
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='#2C3E50', label='Ausente (<10%)'),
    Patch(facecolor='#E67E22', label='Bajo (10-50%)'),
    Patch(facecolor='#F1C40F', label='Medio (50-80%)'),
    Patch(facecolor='#27AE60', label='Alto (>80%)')
]
ax2.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.01, 1),
           title='Categoría', fontsize=10, title_fontsize=11)

ax2.set_title('HUELLA DACTILAR DE MOTIVOS POR SUBFAMILIA\n(Presencia categorizada)',
              fontsize=14, fontweight='bold', pad=15)
ax2.set_ylabel('')
plt.xticks(rotation=45, ha='right', fontsize=9)
plt.yticks(fontsize=12)
plt.tight_layout()
fig2.savefig(os.path.join(OUT, 'HEATMAP_Huella_Dactilar.png'), bbox_inches='tight')
print("  -> HEATMAP_Huella_Dactilar.png")
plt.close()

# ── BARPLOT: Motivos ausentes por subfamilia ──
print("[PLOT 3/3] Barplot de motivos ausentes...")

ausentes_count = {}
for sf in SUBFAMILIAS:
    if sf not in freq_matrix:
        continue
    count = sum(1 for m in motif_cols 
                if freq_matrix[sf][m] < 10 and 
                np.mean([freq_matrix[osf][m] for osf in SUBFAMILIAS if osf != sf and osf in freq_matrix]) >= 50)
    ausentes_count[sf] = count

fig3, ax3 = plt.subplots(figsize=(10, 6))
bars = ax3.bar(ausentes_count.keys(), ausentes_count.values(),
               color=[COLOR_MAP[sf] for sf in ausentes_count.keys()],
               edgecolor='black', linewidth=1.5)

for bar, count in zip(bars, ausentes_count.values()):
    ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
             str(count), ha='center', va='bottom', fontsize=14, fontweight='bold')

ax3.set_ylabel('Número de motivos ausentes', fontsize=12)
ax3.set_title('MOTIVOS AUSENTES O EVITADOS POR SUBFAMILIA\n(Presentes en otras subfamilias pero <10% en esta)',
              fontsize=13, fontweight='bold', pad=15)
ax3.set_xlabel('')
ax3.grid(axis='y', alpha=0.3)
plt.tight_layout()
fig3.savefig(os.path.join(OUT, 'BARPLOT_Motivos_Ausentes.png'), bbox_inches='tight')
print("  -> BARPLOT_Motivos_Ausentes.png")
plt.close()

# =============================================================================
# RESUMEN FINAL
# =============================================================================
banner("ANÁLISIS COMPLETADO")
print(f"Resultados en: {OUT}")
print(f"\nArchivos generados:")
for f in sorted(os.listdir(OUT)):
    size = os.path.getsize(os.path.join(OUT, f))
    print(f"  {f} ({size:,} bytes)")
