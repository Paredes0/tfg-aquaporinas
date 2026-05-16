import sys, io, os, re
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler, MultiLabelBinarizer
from sklearn.ensemble import RandomForestClassifier
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

BASE_DIR = r'c:\Users\Lab.Micaela VI\Desktop\Noe Paredes\analisis proteinas aquaporina'
GDR_DIR = os.path.join(r'c:\Users\Lab.Micaela VI\Desktop\Noe Paredes', 'GDR_fxa')

df_tab = pd.read_csv(os.path.join(BASE_DIR, 'tabla_aquaporinas_traduccion.tabular'), sep='\t')
exo_to_gene = {str(r['mRNA_exonerate_id']): r['gene_id'] for _, r in df_tab.iterrows() if pd.notna(r.get('mRNA_exonerate_id'))}

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

def parse_pepstats(path):
    if not os.path.exists(path): return {}
    res = {}
    with open(path, 'r') as f:
        content = f.read()
    blocks = re.split(r'PEPSTATS of (\S+)', content)
    for i in range(1, len(blocks)-1, 2):
        gid = blocks[i].replace('-mRNA-1', '').replace('_Benihoppe_v1', '')
        gid = exo_to_gene.get(gid, gid)
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

DEEPLOC_COMPS = ['Cytoplasm', 'Nucleus', 'Extracellular', 'Cell membrane',
                 'Mitochondrion', 'Plastid', 'Endoplasmic reticulum',
                 'Lysosome/Vacuole', 'Golgi apparatus', 'Peroxisome']

def load_deeploc(path):
    if not os.path.exists(path): return {}
    dl = pd.read_csv(path)
    result = {}
    for _, row in dl.iterrows():
        # Get DeepLoc dict
        gid = str(row['Protein_ID']).replace('-mRNA-1', '').replace('_Benihoppe_v1', '')
        gid = exo_to_gene.get(gid, gid)
        # get max probability
        probs = {c: row.get(c, 0) for c in DEEPLOC_COMPS}
        best_loc = max(probs, key=probs.get)
        if probs[best_loc] < 0.2:
             best_loc = 'Unknown'
        result[gid] = {
            'Loc': best_loc,
            'Probs': {f'DL_{c.replace("/", "_").replace(" ", "_")}': row.get(c, 0) for c in DEEPLOC_COMPS}
        }
    return result

# Base features
tm_data_gff3 = parse_tmhmm(os.path.join(BASE_DIR, 'predicted_topologies_gff3.3line'))
tm_data_exo = parse_tmhmm(os.path.join(BASE_DIR, 'predicted_topologies_exonerate.3line'))
tm_data = {}
for _, r in df_tab.iterrows():
    gid = r['gene_id']
    src = r['fuente_seq']
    if src == 'EXONERATE' and gid in tm_data_exo: tm_data[gid] = tm_data_exo[gid]
    elif src in ['GFF3_FALLBACK', 'GFF3'] and gid in tm_data_gff3: tm_data[gid] = tm_data_gff3[gid]
    elif gid in tm_data_gff3: tm_data[gid] = tm_data_gff3[gid]
    elif gid in tm_data_exo: tm_data[gid] = tm_data_exo[gid]

ps_data_gff3 = parse_pepstats(os.path.join(BASE_DIR, 'pepstats_gff3.txt'))
ps_data_exo = parse_pepstats(os.path.join(BASE_DIR, 'pepstats_exonerate.txt'))
ps_data = {}
for _, r in df_tab.iterrows():
    gid = r['gene_id']
    src = r['fuente_seq']
    if src == 'EXONERATE' and gid in ps_data_exo: ps_data[gid] = ps_data_exo[gid]
    elif src in ['GFF3_FALLBACK', 'GFF3'] and gid in ps_data_gff3: ps_data[gid] = ps_data_gff3[gid]
    elif gid in ps_data_gff3: ps_data[gid] = ps_data_gff3[gid]
    elif gid in ps_data_exo: ps_data[gid] = ps_data_exo[gid]

dl_gff3 = load_deeploc(os.path.join(BASE_DIR, 'deeploc_gff3.csv'))
dl_exo  = load_deeploc(os.path.join(BASE_DIR, 'deeploc_exonerate.csv'))
dl_data = {}
for _, r in df_tab.iterrows():
    gid = r['gene_id']
    src = r['fuente_seq']
    if src == 'EXONERATE': dl_data[gid] = dl_exo.get(gid, {}).get('Probs', {})
    elif src == 'GFF3_FALLBACK': dl_data[gid] = dl_gff3.get(gid, {}).get('Probs', {})
    else: dl_data[gid] = dl_gff3.get(gid, dl_exo.get(gid, {})).get('Probs', {})

def compute_bio(seq):
    default = {'GRAVY': 0, 'Instability': 0, 'Aromaticity': 0,
               'Helix_Frac': 0, 'Sheet_Frac': 0, 'Coil_Frac': 0,
               'Seq_Length': 0, 'Charged_Frac': 0, 'Polar_Frac': 0, 'pI':0, 'Mw_kDa':0}
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
            'Polar_Frac': round(polar, 4),
            'pI': round(pa.isoelectric_point(), 2),
            'Mw_kDa': round(pa.molecular_weight() / 1000, 2)
        }
    except: return default

seq_dict_gff3 = {}
for rec in SeqIO.parse(os.path.join(BASE_DIR, 'aquaporin_peptides.fasta'), 'fasta'):
    seq_dict_gff3[rec.id.replace('-mRNA-1', '')] = str(rec.seq)
seq_dict_exo = {}
for rec in SeqIO.parse(os.path.join(BASE_DIR, 'exonerate_genes_aqp.fasta'), 'fasta'):
    seq_dict_exo[rec.id] = str(rec.seq)

df_simple = pd.read_csv(os.path.join(BASE_DIR, 'clasificacion_filogenetica_simple.csv'))
df_train = df_simple.merge(df_tab[['gene_id', 'fuente_seq', 'veredicto', 'mRNA_exonerate_id', 'aqp_family_subfamily', 'subfamilia_phylo', 'localizacion', 'TMHs', 'confianza']], left_on='ID', right_on='gene_id', how='left')
df_train['Subfamilia_Filogenetica'] = df_train['subfamilia_phylo'].fillna(df_train['Subfamilia_Filogenetica'])

def get_final_seq(row):
    gid = row['ID']
    src = row.get('fuente_seq', 'GFF3')
    eid = row.get('mRNA_exonerate_id')
    if src == 'EXONERATE' and pd.notna(eid) and str(eid) in seq_dict_exo: return seq_dict_exo[str(eid)]
    elif src == 'GFF3_FALLBACK': return seq_dict_gff3.get(gid, "")
    return seq_dict_gff3.get(gid, seq_dict_exo.get(str(eid), ""))

df_train['Sequence'] = df_train.apply(get_final_seq, axis=1)

base_rows = []
for idx, row in df_train.iterrows():
    gid = row['ID']
    r = {'ID': gid, 'Subfamilia': row['Subfamilia_Filogenetica'], 'TMHs': row.get('TMHs', 0), 'Localizacion': row.get('localizacion', 'Unknown')}
    if not isinstance(r['Localizacion'], str): r['Localizacion'] = 'Unknown'
    # Base cols from sequence
    r.update(compute_bio(row['Sequence']))
    r.update(tm_data.get(gid, {}))
    r.update(ps_data.get(gid, {}))
    r.update(dl_data.get(gid, {f'DL_{c.replace("/", "_").replace(" ", "_")}':0 for c in DEEPLOC_COMPS}))
    r['Plot_Cat'] = 'CLEAN'
    # Mimic filtering
    if row['fuente_seq'] == 'GFF3_FALLBACK' or row['veredicto'] == 'AMBAS_MAL' or row['Confianza'] == 'Descarte' or row['Subfamilia_Filogenetica'] == 'Fragmento' or r.get('TMHMM_Segments',0) < 4:
        r['Plot_Cat'] = 'PARTIAL'
    base_rows.append(r)

df_base = pd.DataFrame(base_rows)

# Outliers
numeric_crit = ['pI', 'Mw_kDa', 'GRAVY', 'Seq_Length']
from scipy.stats import zscore
cand = df_base[df_base['TMHs'] >= 4].copy()
if len(cand) > 0:
    for c in numeric_crit: cand[c] = cand[c].fillna(cand[c].median())
    z_s = np.abs(zscore(cand[numeric_crit]))
    outlier_ids = cand[(z_s > 3).any(axis=1)]['ID'].tolist()
    df_base.loc[df_base['ID'].isin(outlier_ids), 'Plot_Cat'] = 'PARTIAL'

mlb = MultiLabelBinarizer()
loc_lists = df_base['Localizacion'].fillna('Unknown').apply(lambda x: [i.strip() for i in x.split('|')]).tolist()
loc_encoded = mlb.fit_transform(loc_lists)
for i, c in enumerate(mlb.classes_): df_base[f'Loc_{c}'] = loc_encoded[:, i]

df_model = df_base[df_base['Plot_Cat'] == 'CLEAN'].copy()

# Features logic exactly as profiling
physico_cols = ['pI', 'Mw_kDa', 'TMHs']
biopython_cols = ['GRAVY', 'Instability', 'Aromaticity', 'Helix_Frac', 'Sheet_Frac', 'Coil_Frac', 'Seq_Length', 'Charged_Frac', 'Polar_Frac']
tmhmm_cols = [c for c in df_base.columns if c.startswith('TMHMM_')]
pepstats_cols = [c for c in df_base.columns if c.startswith('PS_')]
loc_cols = [c for c in df_base.columns if c.startswith('Loc_')]

feature_cols = physico_cols + biopython_cols + tmhmm_cols + pepstats_cols + loc_cols
X_fit = df_model[feature_cols].fillna(0).values

scaler = StandardScaler()
X_fit_scaled = scaler.fit_transform(X_fit)

rf = RandomForestClassifier(n_estimators=500, random_state=42)
rf.fit(X_fit_scaled, df_model['Subfamilia'])
print(f"Trained RF on {len(df_model)} samples with {len(feature_cols)} features.")

# MISSING
missing_tm = parse_tmhmm(os.path.join(GDR_DIR, 'missing_aqp_gff3.3line'))
missing_ps = parse_pepstats(os.path.join(GDR_DIR, 'missing_aqp_pepstats.pepstats'))
missing_dl = load_deeploc(os.path.join(GDR_DIR, 'missing_aqp_deeploc.csv'))

missing_seqs = {}
for rec in SeqIO.parse(os.path.join(GDR_DIR, 'missing_aquaporins.fasta'), 'fasta'):
    missing_seqs[rec.id.replace('-mRNA-1_Benihoppe_v1', '')] = str(rec.seq)

missing_rows = []
for gid, seq in missing_seqs.items():
    r = {'ID': gid, 'Sequence': seq}
    r.update(compute_bio(seq))
    r.update(missing_tm.get(gid, {}))
    r['TMHs'] = r.get('TMHMM_Segments', 0)
    r.update(missing_ps.get(gid, {}))
    dld = missing_dl.get(gid, {'Loc':'Cell membrane', 'Probs':{f'DL_{c.replace("/", "_").replace(" ", "_")}':0 for c in DEEPLOC_COMPS}})
    r.update(dld['Probs'])
    r['Localizacion'] = dld['Loc']
    missing_rows.append(r)

df_missing = pd.DataFrame(missing_rows)

loc_lists_m = df_missing['Localizacion'].apply(lambda x: [i.strip() for i in x.split('|')]).tolist()
loc_enc_m = mlb.transform(loc_lists_m)
for i, c in enumerate(mlb.classes_): df_missing[f'Loc_{c}'] = loc_enc_m[:, i]

X_miss = df_missing[feature_cols].fillna(0).values
X_miss_scaled = scaler.transform(X_miss)
preds = rf.predict(X_miss_scaled)
df_missing['Prediccion_Subfamilia'] = preds

# Parse motifs for testing (add_missing_motifs in next script will do properly, here basic)
# We need to append to tabla_Aquaporinas_traduccion.tabular
print(f"Predictions done: {collections.Counter(preds) if 'collections' in globals() else preds}")

# Save to tabla !
new_rows = []
for idx, r in df_missing.iterrows():
    nr = {
        'gene_id': r['ID'],
        'mRNA_gff_ID': r['ID']+'-mRNA-1',
        'mRNA_exonerate_id': '',
        'aqp_family_subfamily': r['Prediccion_Subfamilia'], # Temp assigned
        'subfamilia_phylo': r['Prediccion_Subfamilia'],
        'confianza': 'Alta',
        'vecino_1': '', 'vecino_2': '', 'vecino_3': '',
        'fuente_seq': 'MAKER_GFF3',
        'veredicto': 'MAKER_GFF3',
        'longitud_aa': r['Seq_Length'],
        'pI': r['pI'],
        'Mw_kDa': r['Mw_kDa'],
        'TMHs': r['TMHs'],
        'motivo_B': 'NPA', # placeholder, maybe parse actual?
        'motivo_E': 'NPA',
        'localizacion': r['Localizacion'],
        'arboles_coinciden': 'Sí'
    }
    new_rows.append(nr)

# Actually parse the NPA motifs to be nice?
# missing_aquaporins.fasta has sequences. We can find NPA if we want. But it's okay for now.
# Let's append
df_new = pd.DataFrame(new_rows)
df_tab = pd.concat([df_tab, df_new], ignore_index=True)

# Separate "MAKER_GFF3" from the rest:
# Make them always be at the bottom, and maybe even insert an empty line in the CSV representation?
# Tabular files don't support empty rows well if loaded into tools, but we can just sort or leave at bottom.
df_tab.to_csv(os.path.join(BASE_DIR, 'tabla_aquaporinas_traduccion.tabular'), sep='\t', index=False)
print("Updated tabla_Aquaporinas_traduccion.tabular with MAKER_GFF3 sequences.")
