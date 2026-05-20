#!/usr/bin/env python3
"""
15_homeolog_efp_viewer.py — Generate interactive eFP-like expression viewer

Creates a standalone HTML file with an SVG strawberry plant illustration
where tissues are colored by aquaporin homeolog group expression levels (TPM).

SVG tissue group IDs:
  g230  → green_fruit
  g4    → red_fruit
  g48   → leaf (leaves)
  g133  → roots
  g144  → crown
  path1 → outline (non-interactive)

Usage:
    python3 15_homeolog_efp_viewer.py

Reads from:
    results/homeolog_analysis/collapsed_tpm.csv
    homeolog_groups.tsv
    fxa_vectorizado.svg

Outputs:
    results/homeolog_analysis/efp_viewer_homeologs.html
"""

import pandas as pd
import numpy as np
import json
import os
import sys
import math
import base64

# ---- Configuration ----------------------------------------------------------
PROJECT_DIR  = os.environ.get("PROJECT_DIR", os.getcwd())
TPM_FILE     = os.path.join(PROJECT_DIR, "results", "homeolog_analysis", "collapsed_tpm.csv")
BASAL_TPM_FILE = os.path.join(PROJECT_DIR, "results", "basal_aquaporins", "basal_aquaporins_tpm.csv")
HG_TABLE     = os.path.join(PROJECT_DIR, "homeolog_groups.tsv")
SVG_FILE     = os.path.join(PROJECT_DIR, "fxa_vectorizado.svg")
OUTPUT_HTML  = os.path.join(PROJECT_DIR, "results", "homeolog_analysis", "efp_viewer_homeologs.html")

print("# =============================================================")
print("# GENERATING eFP-like VIEWER (HOMEOLOG GROUPS)")
print("# =============================================================")

# ---- Load SVG ----------------------------------------------------------------
if not os.path.exists(SVG_FILE):
    print(f"ERROR: SVG file not found: {SVG_FILE}")
    sys.exit(1)

with open(SVG_FILE, "r", encoding="utf-8") as f:
    svg_content = f.read()

# Clean SVG: remove xml declaration if present
if svg_content.startswith("<?xml"):
    svg_content = svg_content[svg_content.index("?>") + 2:].strip()

print(f"# Loaded SVG: {SVG_FILE} ({len(svg_content)} chars)")

# ---- Load data --------------------------------------------------
tpm = pd.read_csv(TPM_FILE, index_col=0)
basal_tpm = pd.read_csv(BASAL_TPM_FILE, index_col=0) if os.path.exists(BASAL_TPM_FILE) else None
hg = pd.read_csv(HG_TABLE, sep="\t")

print(f"# Loaded collapsed TPM: {tpm.shape[0]} groups x {tpm.shape[1]} samples")
print(f"# Loaded homeolog groups: {len(hg)} genes")

# Restricción a las 121 funcionales (8 candidatas a reanotación excluidas;
# la reanotación no se aborda en el TFG). homeolog_groups.tsv ya debería estar
# restringido a 121 genes / 32 grupos; aplicamos filtro defensivo y replicamos
# el filtro sobre el TPM individual basal antes de incluirlo en el visor.
n_hg_pre = len(hg)
if "needs_reannotation" in hg.columns:
    hg = hg[hg["needs_reannotation"].astype(str).str.upper() != "TRUE"]
elif "is_partial" in hg.columns:
    hg = hg[~hg["is_partial"].astype(str).str.lower().isin(["yes", "true"])]
elif "annotation_source" in hg.columns:
    hg = hg[~hg["annotation_source"].isin(["GFF3_FALLBACK", "MAKER_GFF3"])]
if len(hg) < n_hg_pre:
    print(f"# Excluded {n_hg_pre - len(hg)} reannotation candidates → "
          f"{len(hg)} functional homeologs in {hg['homeolog_group'].nunique()} groups")

# El TPM individual basal también puede contener genes candidatos; filtramos.
if basal_tpm is not None and "needs_reannotation" in basal_tpm.columns:
    n_basal_pre = len(basal_tpm)
    basal_tpm = basal_tpm[basal_tpm["needs_reannotation"].astype(str).str.upper() != "TRUE"]
    if len(basal_tpm) < n_basal_pre:
        print(f"# Filtered basal TPM to {len(basal_tpm)}/{n_basal_pre} functional aquaporins")

# Map tissues
tissue_map = {
    "green_fruit": [c for c in tpm.columns if "Green" in c],
    "red_fruit": [c for c in tpm.columns if "Red" in c],
    "crown": [c for c in tpm.columns if "Crown" in c],
    "leaf": [c for c in tpm.columns if "LeafCtrl" in c],
    "roots": [c for c in tpm.columns if "RootsCtrl" in c],
    "aux_bud": [c for c in tpm.columns if "AuxBud" in c]
}

if basal_tpm is not None:
    tissue_map_basal = {
        "green_fruit": [c for c in basal_tpm.columns if "Green" in c],
        "red_fruit": [c for c in basal_tpm.columns if "Red" in c],
        "crown": [c for c in basal_tpm.columns if "Crown" in c],
        "leaf": [c for c in basal_tpm.columns if "LeafCtrl" in c],
        "roots": [c for c in basal_tpm.columns if "RootsCtrl" in c],
        "aux_bud": [c for c in basal_tpm.columns if "AuxBud" in c]
    }
else:
    tissue_map_basal = {}

# Build expression dataset
expression_data = {}
for hg_id in tpm.index:
    if not str(hg_id).startswith("HG-"): continue
    
    expression_data[hg_id] = {}
    
    for tissue, cols in tissue_map.items():
        if len(cols) == 0: continue
        vals = tpm.loc[hg_id, cols].dropna().values
        expression_data[hg_id][tissue] = {
            "mean_tpm": round(float(np.mean(vals)), 2) if len(vals) > 0 else 0,
            "sd_tpm": round(float(np.std(vals, ddof=1)), 2) if len(vals) > 1 else None,
            "n": len(vals)
        }
    
    # Add individual gene breakdown under a special key
    indiv_data = {}
    if basal_tpm is not None:
        group_rows = hg[hg["homeolog_group"] == hg_id]
        for _, row in group_rows.iterrows():
            gene = row["gene_id"]
            subg = row["subgenome"]
            if gene in basal_tpm.index:
                indiv_data[gene] = {"subgenome": subg}
                for t_name, c_names in tissue_map_basal.items():
                    if len(c_names) == 0: continue
                    a_cols = [c for c in c_names if c in basal_tpm.columns]
                    if not a_cols:
                        indiv_data[gene][t_name] = 0
                        continue
                    v = basal_tpm.loc[gene, a_cols].dropna().values
                    indiv_data[gene][t_name] = round(float(np.mean(v)), 2) if len(v) > 0 else 0
    
    expression_data[hg_id]["individual"] = indiv_data

# Build gene metadata
gene_metadata = {}
valid_hgs = [gid for gid in tpm.index if str(gid).startswith("HG-")]
for hg_id in valid_hgs:
    group_rows = hg[hg["homeolog_group"] == hg_id]
    if not group_rows.empty:
        row = group_rows.iloc[0]
        sub_subfamily = row["sub_subfamily"]
        family = row["family"]
        name_str = f"{hg_id} ({sub_subfamily})"
    else:
        sub_subfamily = "Unknown"
        family = "Unknown"
        name_str = hg_id
        
    gene_metadata[hg_id] = {
        "name": name_str,
        "sub_subfamily": sub_subfamily,
        "family": family,
        "fuente_seq": "collapsed_homeologs",
        "TMHs": 6,
        "needs_reannotation": False
    }

# Build sub_subfamily list
subfamilies = sorted(set(m["sub_subfamily"] for m in gene_metadata.values()))

# Build gene list per sub_subfamily
genes_by_sub_subfamily = {}
for gid, meta in gene_metadata.items():
    sf = meta["sub_subfamily"]
    if sf not in genes_by_sub_subfamily:
        genes_by_sub_subfamily[sf] = []
    genes_by_sub_subfamily[sf].append({"id": gid, "name": meta["name"]})

for sf in genes_by_sub_subfamily:
    genes_by_sub_subfamily[sf].sort(key=lambda x: x["id"])

# Calculate global max TPM for scale
all_tpms = []
for gid_data in expression_data.values():
    for t_key, tdata in gid_data.items():
        if t_key != "individual":
            all_tpms.append(tdata["mean_tpm"])
global_max_tpm = max(all_tpms) if all_tpms else 100

data_json = json.dumps({
    "expression": expression_data,
    "metadata": gene_metadata,
    "subfamilies": subfamilies,
    "genes_by_sub_subfamily": genes_by_sub_subfamily,
    "global_max_tpm": round(global_max_tpm, 2)
}, indent=None)

# ---- Generate HTML -----------------------------------------------------------
html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Fragaria x ananassa Aquaporin eFP Viewer</title>
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=Source+Serif+4:opsz,wght@8..60,400;8..60,600;8..60,700&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet">
<style>
  :root {{
    /* Surfaces — warm ivory / stone palette */
    --bg-page:        #F7F4EE;
    --bg-surface:     #FFFFFF;
    --bg-subtle:      #F1EDE5;
    --bg-soft:        #FAF7F1;
    --bg-tissue-base: #C9C3B2;
    /* Lines */
    --border:         #D6CEC1;
    --border-soft:    #E8E2D5;
    --border-strong:  #A89E8C;
    /* Text */
    --text-primary:   #1F2730;
    --text-secondary: #4A4842;
    --text-muted:     #8A8377;
    /* Accent — botanical dark green */
    --accent:         #2E5A3E;
    --accent-soft:    #E6EEE6;
    --accent-warm:    #8B4A2B;
    /* Aliases used by existing code (kept for back-compat) */
    --bg-primary:     var(--bg-page);
    --bg-secondary:   var(--bg-surface);
    --bg-card:        var(--bg-surface);
    --bg-hover:       var(--bg-subtle);
    --accent-primary: var(--accent);
    --accent-secondary: var(--accent);
    --border-active: var(--accent);
    /* Misc */
    --shadow-soft:    0 1px 2px rgba(31,39,48,0.04), 0 1px 2px rgba(31,39,48,0.06);
    --shadow-pop:     0 4px 14px rgba(31,39,48,0.08);
    --shadow:         var(--shadow-soft);
    --radius:         4px;
    --radius-sm:      3px;
    /* Typography */
    --font-serif:     'Source Serif 4', 'Source Serif Pro', Georgia, 'Times New Roman', serif;
    --font-sans:      'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
    --font-mono:      'JetBrains Mono', 'SF Mono', Consolas, 'Courier New', monospace;
  }}

  * {{ margin: 0; padding: 0; box-sizing: border-box; }}

  body {{
    font-family: var(--font-sans);
    background: var(--bg-page);
    color: var(--text-primary);
    min-height: 100vh;
    overflow-x: hidden;
    font-feature-settings: 'cv11', 'ss01', 'ss03';
    -webkit-font-smoothing: antialiased;
    -moz-osx-font-smoothing: grayscale;
  }}

  .app-header {{
    background: var(--bg-surface);
    border-bottom: 1px solid var(--border);
    padding: 18px 32px;
    display: flex;
    align-items: center;
    gap: 20px;
  }}

  .app-header h1 {{
    font-family: var(--font-serif);
    font-size: 1.4rem;
    font-weight: 600;
    letter-spacing: -0.005em;
    color: var(--text-primary);
    line-height: 1.2;
  }}

  .app-header h1 .sp-italic {{
    font-style: italic;
  }}

  .app-header .subtitle {{
    color: var(--text-muted);
    font-size: 0.82rem;
    font-weight: 400;
    margin-top: 2px;
    letter-spacing: 0.01em;
  }}

  .logo-icon {{
    width: 44px;
    height: 44px;
    background: var(--bg-surface);
    border: 1px solid var(--border-strong);
    border-radius: var(--radius);
    display: flex;
    align-items: center;
    justify-content: center;
    font-family: var(--font-serif);
    font-style: italic;
    font-weight: 600;
    font-size: 1.15rem;
    color: var(--accent);
    letter-spacing: -0.02em;
    flex-shrink: 0;
  }}

  .main-container {{
    display: grid;
    grid-template-columns: 320px 1fr 340px;
    gap: 0;
    height: calc(100vh - 82px);
  }}

  .panel {{
    background: var(--bg-surface);
    border-right: 1px solid var(--border);
    overflow-y: auto;
    padding: 22px 20px;
  }}

  .panel:last-child {{
    border-right: none;
    border-left: 1px solid var(--border);
  }}

  .panel-title {{
    font-family: var(--font-serif);
    font-size: 0.78rem;
    font-weight: 700;
    text-transform: uppercase;
    letter-spacing: 0.12em;
    color: var(--text-secondary);
    margin-bottom: 14px;
    padding-bottom: 8px;
    border-bottom: 1px solid var(--border);
  }}

  .control-group {{
    margin-bottom: 20px;
  }}

  .control-group label {{
    display: block;
    font-size: 0.82rem;
    font-weight: 500;
    color: var(--text-secondary);
    margin-bottom: 6px;
  }}

  select, input[type="range"] {{
    width: 100%;
    background: var(--bg-surface);
    border: 1px solid var(--border);
    color: var(--text-primary);
    padding: 9px 12px;
    border-radius: var(--radius-sm);
    font-family: inherit;
    font-size: 0.85rem;
    outline: none;
    transition: border-color 0.15s, box-shadow 0.15s;
    appearance: none;
    -webkit-appearance: none;
    cursor: pointer;
  }}

  select {{
    background-image: url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' width='12' height='12' viewBox='0 0 12 12'%3E%3Cpath fill='%238A8377' d='M6 8L1 3h10z'/%3E%3C/svg%3E");
    background-repeat: no-repeat;
    background-position: right 12px center;
    padding-right: 32px;
  }}

  select:focus, select:hover {{
    border-color: var(--accent);
    box-shadow: 0 0 0 2px var(--accent-soft);
  }}

  .gene-list {{
    max-height: 420px;
    overflow-y: auto;
    border: 1px solid var(--border);
    border-radius: var(--radius-sm);
    background: var(--bg-surface);
  }}

  .gene-item {{
    padding: 9px 12px 9px 14px;
    cursor: pointer;
    font-size: 0.8rem;
    border-bottom: 1px solid var(--border-soft);
    display: flex;
    justify-content: space-between;
    align-items: center;
    transition: background 0.12s;
    border-left: 3px solid transparent;
  }}

  .gene-item:last-child {{
    border-bottom: none;
  }}

  .gene-item:hover {{
    background: var(--bg-soft);
  }}

  .gene-item.active {{
    background: var(--accent-soft);
    border-left-color: var(--accent);
  }}

  .gene-item .gene-name {{
    font-weight: 600;
    color: var(--text-primary);
    letter-spacing: -0.005em;
  }}

  .gene-item .gene-id {{
    font-size: 0.7rem;
    color: var(--text-muted);
    font-family: var(--font-mono);
    margin-top: 1px;
  }}

  .gene-item .reannotation-badge {{
    background: var(--bg-surface);
    color: var(--accent-warm);
    border: 1px solid var(--accent-warm);
    font-size: 0.58rem;
    padding: 2px 6px;
    border-radius: var(--radius-sm);
    font-weight: 600;
    letter-spacing: 0.04em;
    text-transform: uppercase;
  }}

  .plant-viewer {{
    display: flex;
    align-items: center;
    justify-content: center;
    background: var(--bg-page);
    position: relative;
    overflow: hidden;
  }}

  .plant-svg-container {{
    width: 100%;
    height: 100%;
    display: flex;
    align-items: flex-start;
    justify-content: center;
    padding: 0;
    overflow: hidden;
  }}

  .plant-svg-container svg {{
    width: 185%;
    height: 185%;
    min-width: 185%;
    min-height: 185%;
    transform: translateY(-10%);
    transform-origin: top center;
  }}

  /* SVG tissue styling — applied via JS to the real SVG groups */
  #outline-path {{
    pointer-events: none;
  }}

  /* Info panel */
  .info-section {{
    background: var(--bg-surface);
    border-radius: var(--radius-sm);
    padding: 14px 16px 16px;
    margin-bottom: 14px;
    border: 1px solid var(--border);
  }}

  .info-section h3 {{
    font-family: var(--font-serif);
    font-size: 0.78rem;
    font-weight: 700;
    text-transform: uppercase;
    letter-spacing: 0.1em;
    color: var(--accent);
    margin-bottom: 10px;
    padding-bottom: 6px;
    border-bottom: 1px solid var(--border-soft);
  }}

  .info-row {{
    display: flex;
    justify-content: space-between;
    align-items: baseline;
    padding: 5px 0;
    font-size: 0.8rem;
    border-bottom: 1px solid var(--border-soft);
  }}

  .info-row:last-child {{
    border-bottom: none;
  }}

  .info-row .label {{
    color: var(--text-muted);
    font-size: 0.75rem;
  }}

  .info-row .value {{
    color: var(--text-primary);
    font-weight: 500;
    text-align: right;
    font-variant-numeric: tabular-nums;
  }}

  /* Top control bar floating over plant viewer */
  .viewer-topbar {{
    position: absolute;
    top: 14px;
    left: 50%;
    transform: translateX(-50%);
    background: var(--bg-surface);
    border: 1px solid var(--border);
    border-radius: var(--radius-sm);
    padding: 8px 14px;
    display: flex;
    align-items: center;
    gap: 16px;
    box-shadow: var(--shadow-soft);
    z-index: 10;
  }}

  .topbar-section {{
    display: flex;
    align-items: center;
    gap: 8px;
  }}

  .topbar-divider {{
    width: 1px;
    height: 22px;
    background: var(--border);
  }}

  .topbar-label {{
    font-family: var(--font-serif);
    font-size: 0.72rem;
    font-weight: 600;
    color: var(--text-secondary);
    text-transform: uppercase;
    letter-spacing: 0.08em;
  }}

  /* Color scale legend (lives in topbar now) */
  .scale-bar {{
    width: 200px;
    height: 12px;
    border-radius: 2px;
    border: 1px solid var(--border);
    background: linear-gradient(90deg,
      #C9C3B2 0%,
      #C9C3B2 3%,
      #FFF5C2 5%,
      #FEE08B 18%,
      #FDB863 32%,
      #F39434 48%,
      #E5601F 62%,
      #C42E1A 78%,
      #8E1414 92%,
      #5A0A0A 100%);
  }}

  .scale-label {{
    font-size: 0.7rem;
    color: var(--text-secondary);
    white-space: nowrap;
    font-variant-numeric: tabular-nums;
    font-weight: 500;
  }}

  /* Scale-mode descriptive help (sidebar) */
  .scale-help {{
    font-size: 0.7rem;
    color: var(--text-muted);
    line-height: 1.6;
    padding: 10px 12px;
    background: var(--bg-subtle);
    border-left: 2px solid var(--border-strong);
    border-radius: 2px;
  }}

  .scale-help b {{
    color: var(--text-secondary);
    font-weight: 600;
  }}

  /* Segmented control for scale mode */
  .mode-toggle {{
    display: flex;
    background: var(--bg-subtle);
    border-radius: var(--radius-sm);
    border: 1px solid var(--border);
    padding: 2px;
    gap: 2px;
  }}

  .mode-btn {{
    padding: 5px 12px;
    text-align: center;
    font-size: 0.73rem;
    font-weight: 500;
    cursor: pointer;
    transition: background 0.12s, color 0.12s;
    color: var(--text-secondary);
    background: transparent;
    border: none;
    border-radius: 2px;
    font-family: inherit;
    letter-spacing: 0.01em;
  }}

  .mode-btn.active {{
    background: var(--bg-surface);
    color: var(--accent);
    font-weight: 600;
    box-shadow: var(--shadow-soft);
  }}

  .mode-btn:hover:not(.active) {{
    color: var(--text-primary);
  }}

  /* Expression bar chart in info panel */
  .expr-bar-container {{
    margin-top: 4px;
  }}

  .expr-bar-row {{
    display: flex;
    align-items: center;
    margin-bottom: 5px;
    gap: 8px;
    border-radius: 2px;
    padding: 2px 2px;
    transition: background 0.12s;
  }}

  .expr-bar-row:hover {{
    background: var(--bg-soft);
  }}

  .expr-bar-tissue {{
    width: 92px;
    font-size: 0.72rem;
    color: var(--text-secondary);
    text-align: right;
    font-weight: 500;
  }}

  .expr-bar-wrapper {{
    flex: 1;
    height: 16px;
    background: var(--bg-subtle);
    border-radius: 2px;
    overflow: hidden;
    position: relative;
    border: 1px solid var(--border-soft);
  }}

  .expr-bar {{
    height: 100%;
    border-radius: 0;
    transition: width 0.4s ease;
    min-width: 2px;
  }}

  .expr-bar-value {{
    position: absolute;
    right: 6px;
    top: 50%;
    transform: translateY(-50%);
    font-size: 0.65rem;
    color: var(--text-primary);
    font-weight: 600;
    font-variant-numeric: tabular-nums;
    text-shadow: 0 0 3px rgba(255,255,255,0.85);
  }}

  .n1-warning {{
    font-size: 0.68rem;
    color: var(--accent-warm);
    font-style: italic;
    margin-top: 6px;
    padding: 4px 8px;
    background: var(--bg-subtle);
    border-left: 2px solid var(--accent-warm);
    border-radius: 2px;
  }}

  /* Search */
  .search-input {{
    width: 100%;
    background: var(--bg-surface);
    border: 1px solid var(--border);
    color: var(--text-primary);
    padding: 9px 12px 9px 32px;
    border-radius: var(--radius-sm);
    font-family: inherit;
    font-size: 0.82rem;
    outline: none;
    transition: border-color 0.15s, box-shadow 0.15s;
    margin-bottom: 12px;
    background-image: url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' width='14' height='14' viewBox='0 0 24 24' fill='none' stroke='%238A8377' stroke-width='2'%3E%3Ccircle cx='11' cy='11' r='8'/%3E%3Cline x1='21' y1='21' x2='16.65' y2='16.65'/%3E%3C/svg%3E");
    background-repeat: no-repeat;
    background-position: 10px center;
  }}

  .search-input:focus {{
    border-color: var(--accent);
    box-shadow: 0 0 0 2px var(--accent-soft);
  }}

  .search-input::placeholder {{
    color: var(--text-muted);
  }}

  .no-data-msg {{
    text-align: center;
    padding: 40px 20px;
    color: var(--text-muted);
    font-size: 0.85rem;
    font-style: italic;
  }}

  .scroll-custom::-webkit-scrollbar {{
    width: 8px;
  }}
  .scroll-custom::-webkit-scrollbar-track {{
    background: transparent;
  }}
  .scroll-custom::-webkit-scrollbar-thumb {{
    background: var(--border);
    border-radius: 4px;
    border: 2px solid var(--bg-surface);
  }}
  .scroll-custom::-webkit-scrollbar-thumb:hover {{
    background: var(--border-strong);
  }}

  /* Floating tooltip */
  .svg-tooltip {{
    position: fixed;
    pointer-events: none;
    background: var(--bg-surface);
    border: 1px solid var(--border-strong);
    border-radius: var(--radius-sm);
    padding: 9px 12px;
    font-size: 0.78rem;
    color: var(--text-primary);
    box-shadow: var(--shadow-pop);
    z-index: 1000;
    opacity: 0;
    transition: opacity 0.12s;
    max-width: 230px;
  }}
  .svg-tooltip.visible {{
    opacity: 1;
  }}
  .svg-tooltip .tt-tissue {{
    font-weight: 600;
    color: var(--accent);
    margin-bottom: 3px;
    font-family: var(--font-serif);
    font-size: 0.82rem;
  }}
  .svg-tooltip .tt-tpm {{
    font-size: 1rem;
    font-weight: 700;
    color: var(--text-primary);
    font-variant-numeric: tabular-nums;
  }}
  .svg-tooltip .tt-detail {{
    font-size: 0.7rem;
    color: var(--text-muted);
    margin-top: 3px;
    font-variant-numeric: tabular-nums;
  }}

  /* Export button */
  .export-btn {{
    position: absolute;
    top: 14px;
    right: 14px;
    background: var(--bg-surface);
    border: 1px solid var(--border);
    color: var(--text-secondary);
    padding: 6px 14px;
    border-radius: var(--radius-sm);
    font-size: 0.73rem;
    font-weight: 500;
    cursor: pointer;
    font-family: inherit;
    transition: background 0.15s, color 0.15s, border-color 0.15s;
    z-index: 10;
    letter-spacing: 0.02em;
  }}
  .export-btn:hover {{
    background: var(--accent);
    color: var(--bg-surface);
    border-color: var(--accent);
  }}

  /* Aux bud inset figure */
  .aux-bud-inset {{
    position: absolute;
    bottom: 20px;
    right: 20px;
    background: var(--bg-surface);
    border: 1px solid var(--border);
    border-radius: var(--radius-sm);
    padding: 12px 16px;
    box-shadow: var(--shadow-soft);
    z-index: 5;
    cursor: pointer;
    transition: border-color 0.15s, box-shadow 0.15s;
    text-align: center;
  }}
  .aux-bud-inset:hover {{
    border-color: var(--accent);
    box-shadow: var(--shadow-pop);
  }}
  .aux-bud-inset .aux-label {{
    font-family: var(--font-serif);
    font-size: 0.72rem;
    color: var(--text-secondary);
    margin-bottom: 6px;
    font-weight: 600;
    letter-spacing: 0.04em;
    text-transform: uppercase;
  }}
  .aux-bud-inset .aux-tpm {{
    font-size: 0.82rem;
    color: var(--text-primary);
    font-weight: 700;
    margin-top: 4px;
    font-variant-numeric: tabular-nums;
  }}
  .aux-bud-inset .aux-warning {{
    font-size: 0.6rem;
    color: var(--accent-warm);
    font-style: italic;
    margin-top: 2px;
    letter-spacing: 0.02em;
  }}

  @media (max-width: 1100px) {{
    .main-container {{
      grid-template-columns: 1fr;
      grid-template-rows: auto 1fr auto;
    }}
  }}
</style>
</head>
<body>

<header class="app-header">
  <div class="logo-icon">Fxa</div>
  <div>
    <h1><span class="sp-italic">Fragaria</span> x <span class="sp-italic">ananassa</span> — Aquaporin eFP Viewer</h1>
    <div class="subtitle">Tissue expression atlas · Basal expression (TPM) · Homeolog groups</div>
  </div>
</header>

<div class="main-container">
  <!-- LEFT PANEL: Controls -->
  <div class="panel scroll-custom">
    <div class="panel-title">Gene Selection</div>

    <div class="control-group" id="search-group">
      <input type="text" class="search-input" id="gene-search" placeholder="Search gene..." oninput="filterGenes()">
    </div>

    <div class="gene-list scroll-custom" id="gene-list">
    </div>

    <div class="scale-help" style="margin-top: 14px;">
      <b>Absolute:</b> color scaled to global max TPM across all genes.<br>
      <b>Relative:</b> color scaled to the selected gene's own max — preferred for visual contrast.
    </div>
  </div>

  <!-- CENTER: Plant SVG -->
  <div class="plant-viewer" id="plant-viewer">
    <!-- Top control bar: scale mode + color legend -->
    <div class="viewer-topbar">
      <div class="topbar-section">
        <span class="topbar-label">Scale</span>
        <div class="mode-toggle">
          <button class="mode-btn" onclick="setScale('absolute')" id="scale-abs">Absolute</button>
          <button class="mode-btn active" onclick="setScale('relative')" id="scale-rel">Relative</button>
        </div>
      </div>
      <div class="topbar-divider"></div>
      <div class="topbar-section">
        <span class="topbar-label">Expression</span>
        <span class="scale-label">0</span>
        <div class="scale-bar"></div>
        <span class="scale-label" id="max-tpm-label">— TPM</span>
      </div>
    </div>
    <button class="export-btn" onclick="exportImage()">Export PNG</button>
    <div class="plant-svg-container" id="svg-container">
      {svg_content}
    </div>
    <!-- Aux Bud inset figure -->
    <div class="aux-bud-inset" id="aux-bud-inset" onclick="selectTissue('aux_bud')">
      <div class="aux-label">Axillary Meristem</div>
      <svg width="60" height="50" viewBox="0 0 60 50" xmlns="http://www.w3.org/2000/svg">
        <!-- Stolon -->
        <path d="M5,35 Q20,30 30,32" fill="none" stroke="#7CB342" stroke-width="2.5" stroke-linecap="round"/>
        <!-- Bud body -->
        <ellipse id="aux-bud-shape" cx="38" cy="30" rx="14" ry="10" fill="#C9C3B2"/>
        <!-- Bud leaves -->
        <path d="M32,25 Q38,15 44,25" fill="#6B9E3A" opacity="0.8"/>
        <path d="M30,27 Q36,18 42,26" fill="#7CB342" opacity="0.7"/>
        <!-- Small emerging leaves -->
        <ellipse cx="48" cy="24" rx="6" ry="3" fill="#8BC34A" transform="rotate(-20 48 24)" opacity="0.6"/>
      </svg>
      <div class="aux-tpm" id="aux-bud-tpm">— TPM</div>
      <div class="aux-warning">N = 1 (descriptive)</div>
    </div>
  </div>

  <!-- RIGHT PANEL: Info -->
  <div class="panel scroll-custom">
    <div class="panel-title">Expression Details</div>

    <div id="no-selection" class="no-data-msg">
      ← Select a gene from the list to visualize expression
    </div>

    <div id="gene-info" style="display:none;">
      <div class="info-section">
        <h3>Gene Info</h3>
        <div class="info-row"><span class="label">Name</span><span class="value" id="info-gene-name">-</span></div>
        <div class="info-row"><span class="label">Group ID</span><span class="value" id="info-gene-id">-</span></div>
        <div class="info-row"><span class="label">Sub-subfamily</span><span class="value" id="info-sub_subfamily">-</span></div>
        <div class="info-row"><span class="label">Source</span><span class="value" id="info-source">-</span></div>
        <div class="info-row"><span class="label">TMHs (est.)</span><span class="value" id="info-tmhs">-</span></div>
        <div class="info-row" id="reannotation-row" style="display:none;">
          <span class="label">Reannotation</span>
          <span class="value" style="color:var(--accent-warm)">Needed</span>
        </div>
      </div>

      <div class="info-section">
        <h3>Expression by Tissue (TPM)</h3>
        <div id="expr-bars" class="expr-bar-container"></div>
      </div>

      <div class="info-section" id="tissue-detail" style="display:none;">
        <h3 id="tissue-detail-title">Tissue Detail</h3>
        <div class="info-row"><span class="label">Mean TPM</span><span class="value" id="detail-tpm">-</span></div>
        <div class="info-row"><span class="label">SD</span><span class="value" id="detail-sd">-</span></div>
        <div class="info-row"><span class="label">Replicates</span><span class="value" id="detail-n">-</span></div>
      </div>

      <div class="info-section" id="individual-expression-section" style="display:none;">
        <h3>Individual Homeolog Expression</h3>
        <div id="individual-expression-table-container"></div>
      </div>
    </div>
  </div>
</div>

<!-- Tooltip -->
<div class="svg-tooltip" id="svg-tooltip">
  <div class="tt-tissue" id="tt-tissue"></div>
  <div class="tt-tpm" id="tt-tpm"></div>
  <div class="tt-detail" id="tt-detail"></div>
</div>

<script>
// ---- Data -------------------------------------------------------------------
const DATA = {data_json};

// ---- SVG group ID to tissue mapping -----------------------------------------
// These are the actual group IDs in fxa_vectorizado.svg
const SVG_TISSUE_MAP = {{
  'g230':  'green_fruit',
  'g4':    'red_fruit',
  'g48':   'leaf',
  'g133':  'roots',
  'g144':  'crown'
}};

// Reverse map: tissue → SVG group ID
const TISSUE_SVG_MAP = {{}};
Object.entries(SVG_TISSUE_MAP).forEach(([svgId, tissue]) => {{
  TISSUE_SVG_MAP[tissue] = svgId;
}});

const TISSUE_NAMES = {{
  'leaf': 'Leaf (Control)',
  'roots': 'Roots (Control)',
  'crown': 'Crown',
  'green_fruit': 'Green Fruit',
  'red_fruit': 'Red Fruit',
  'aux_bud': 'Aux. Bud (N=1)'
}};

const TISSUE_COLORS = {{
  'leaf': '#32CD32',
  'roots': '#8B4513',
  'crown': '#DAA520',
  'green_fruit': '#228B22',
  'red_fruit': '#DC143C',
  'aux_bud': '#9370DB'
}};

const SUBFAMILY_COLORS = {{
  'PIP': '#E74C3C',
  'TIP': '#3498DB',
  'NIP': '#2ECC71',
  'SIP': '#F39C12',
  'XIP': '#9B59B6'
}};

// ---- State ------------------------------------------------------------------
// Default to 'relative': global max TPM is dominated by a few root-expressed
// HGs, which compresses every other gene into the low-contrast end of the
// absolute ramp. Relative scales each gene's own range, giving more visual
// information for the common case.
let currentScale = 'relative';
let selectedGeneId = null;
let selectedTissue = null;

// ---- YlOrRd-style color interpolation (contrasted on ivory background) ------
// Hard discontinuity: TPM=0 → taupe grey (clearly darker than ivory bg);
// TPM>0 jumps to vivid pale yellow, then ColorBrewer-like YlOrRd ramp to deep red.
function tpmToColor(tpm, maxTpm) {{
  if (maxTpm <= 0) maxTpm = 1;
  const ratio = Math.min(tpm / maxTpm, 1);
  // Hard floor: real zero (or near-zero) stays as taupe — distinguishable from any expression
  if (tpm <= 0) {{
    return 'rgb(201, 195, 178)';
  }}
  const stops = [
    [0.00, [255, 245, 194]],   // pale yellow — lowest expression
    [0.10, [254, 224, 139]],   // light yellow
    [0.22, [253, 184,  99]],   // gold
    [0.35, [243, 148,  52]],   // orange
    [0.50, [229,  96,  31]],   // burnt orange
    [0.65, [196,  46,  26]],   // red
    [0.80, [142,  20,  20]],   // deep red
    [1.00, [ 90,  10,  10]]    // very dark red — max expression
  ];

  let lower = stops[0], upper = stops[stops.length - 1];
  for (let i = 0; i < stops.length - 1; i++) {{
    if (ratio >= stops[i][0] && ratio <= stops[i + 1][0]) {{
      lower = stops[i];
      upper = stops[i + 1];
      break;
    }}
  }}

  const range = upper[0] - lower[0] || 1;
  const t = (ratio - lower[0]) / range;
  const r = Math.round(lower[1][0] + (upper[1][0] - lower[1][0]) * t);
  const g = Math.round(lower[1][1] + (upper[1][1] - lower[1][1]) * t);
  const b = Math.round(lower[1][2] + (upper[1][2] - lower[1][2]) * t);

  return `rgb(${{r}},${{g}},${{b}})`;
}}

// ---- Initialize SVG ---------------------------------------------------------
function initSVG() {{
  const svg = document.querySelector('#svg-container svg');
  if (!svg) {{
    console.error('SVG not found in container');
    return;
  }}

  // Make SVG responsive
  svg.setAttribute('width', '100%');
  svg.setAttribute('height', '100%');
  svg.style.maxWidth = '100%';
  svg.style.maxHeight = '100%';

  // Set up outline (path1) to be non-interactive — dark warm stroke on ivory bg
  const outline = svg.querySelector('#path1');
  if (outline) {{
    outline.style.pointerEvents = 'none';
    outline.style.fill = 'none';
    outline.style.stroke = 'rgba(74, 72, 66, 0.55)';
    outline.style.strokeWidth = '0.6';
  }}

  // Set up each tissue group
  Object.entries(SVG_TISSUE_MAP).forEach(([svgId, tissue]) => {{
    const group = document.getElementById(svgId);
    if (!group) {{
      console.warn(`SVG group #${{svgId}} not found for tissue '${{tissue}}'`);
      return;
    }}

    // Mark as interactive
    group.style.cursor = 'pointer';
    group.style.transition = 'opacity 0.3s, filter 0.3s';
    group.dataset.tissue = tissue;

    // Hover effects — subtle darken + thin outline (no neon glow)
    group.addEventListener('mouseenter', function() {{
      if (selectedTissue !== tissue) {{
        this.style.filter = 'brightness(0.92) drop-shadow(0 0 1.5px rgba(46, 90, 62, 0.55))';
      }}
    }});
    group.addEventListener('mouseleave', function() {{
      if (selectedTissue !== tissue) {{
        this.style.filter = '';
      }}
    }});

    // Click to select tissue
    group.addEventListener('click', function() {{
      selectTissue(tissue);
    }});
  }});

  // Initial coloring — taupe baseline (no expression), clearly distinct from ivory bg
  Object.entries(SVG_TISSUE_MAP).forEach(([svgId, tissue]) => {{
    colorSVGGroup(svgId, 'rgb(201, 195, 178)');
  }});
}}

// ---- Color an SVG group's paths/elements ------------------------------------
function colorSVGGroup(svgGroupId, color) {{
  const group = document.getElementById(svgGroupId);
  if (!group) return;

  // Color ALL paths, circles, ellipses, rects within this group
  const elements = group.querySelectorAll('path, circle, ellipse, rect, polygon');
  elements.forEach(el => {{
    const elId = el.getAttribute('id') || '';
    // Don't color the outline
    if (elId === 'path1') return;

    el.style.fill = color;
    el.style.transition = 'fill 0.4s ease';
  }});
}}

// ---- Initialize -------------------------------------------------------------
function init() {{
  initSVG();
  populateGeneList();
  document.getElementById('max-tpm-label').textContent = DATA.global_max_tpm.toFixed(0) + ' TPM';

  // Select first gene if available
  const genes = getFilteredGenes();
  if (genes.length > 0) {{
    selectGene(genes[0].id);
  }}
}}

function getFilteredGenes() {{
  const search = document.getElementById('gene-search').value.toLowerCase();
  let genes = [];

  for (const sf of DATA.subfamilies) {{
    genes = genes.concat(DATA.genes_by_sub_subfamily[sf] || []);
  }}

  if (search) {{
    genes = genes.filter(g =>
      g.name.toLowerCase().includes(search) ||
      g.id.toLowerCase().includes(search)
    );
  }}

  return genes;
}}

function populateGeneList() {{
  const list = document.getElementById('gene-list');
  const genes = getFilteredGenes();

  list.innerHTML = '';
  genes.forEach(g => {{
    const meta = DATA.metadata[g.id];
    const div = document.createElement('div');
    div.className = 'gene-item' + (g.id === selectedGeneId ? ' active' : '');
    div.onclick = () => selectGene(g.id);

    let badge = '';
    if (meta && meta.needs_reannotation) {{
      badge = '<span class="reannotation-badge">REANNOT</span>';
    }}

    const sfColor = SUBFAMILY_COLORS[meta?.family] || '#666';
    div.innerHTML = `
      <div>
        <div class="gene-name" style="color:${{sfColor}}">${{g.name}}</div>
        <div class="gene-id">${{g.id}}</div>
      </div>
      ${{badge}}
    `;
    list.appendChild(div);
  }});
}}

// ---- Scale ------------------------------------------------------------------

function setScale(scale) {{
  currentScale = scale;
  document.getElementById('scale-abs').classList.toggle('active', scale === 'absolute');
  document.getElementById('scale-rel').classList.toggle('active', scale === 'relative');

  if (selectedGeneId) updatePlant(selectedGeneId);
}}

// ---- Selection handlers -----------------------------------------------------

function filterGenes() {{
  populateGeneList();
}}

function selectGene(geneId) {{
  selectedGeneId = geneId;
  populateGeneList();
  updatePlant(geneId);
  updateInfoPanel(geneId);
  // Refresh tissue detail panel with the NEW gene's data
  // (was previously stale — anchored to the gene selected when the tissue
  //  was first clicked). See: selectTissue handler.
  if (selectedTissue) {{
    updateTissueDetail(geneId, selectedTissue);
    // Re-apply selected-tissue outline (filter was reset by re-coloring)
    const svgId = TISSUE_SVG_MAP[selectedTissue];
    const group = svgId ? document.getElementById(svgId) : null;
    if (group) {{
      group.style.filter = 'drop-shadow(0 0 1.5px rgba(46, 90, 62, 0.95)) drop-shadow(0 0 0.5px rgba(46, 90, 62, 0.95))';
    }}
  }}
  document.getElementById('no-selection').style.display = 'none';
  document.getElementById('gene-info').style.display = 'block';
}}

function selectTissue(tissue) {{
  selectedTissue = tissue;

  // Update visual selection: reset all, highlight selected with thin dark outline
  Object.entries(SVG_TISSUE_MAP).forEach(([svgId, t]) => {{
    const group = document.getElementById(svgId);
    if (!group) return;
    if (t === tissue) {{
      group.style.filter = 'drop-shadow(0 0 1.5px rgba(46, 90, 62, 0.95)) drop-shadow(0 0 0.5px rgba(46, 90, 62, 0.95))';
    }} else {{
      group.style.filter = '';
    }}
  }});

  // Show tissue detail
  if (selectedGeneId) {{
    updateTissueDetail(selectedGeneId, tissue);
  }}
}}

// ---- Update plant visualization (TAIR10-style) ------------------------------
function updatePlant(geneId) {{
  const expr = DATA.expression[geneId];
  if (!expr) return;

  let maxTpm = DATA.global_max_tpm;
  if (currentScale === 'relative') {{
    // Exclude aux_bud from relative scale (N=1, descriptive only)
    maxTpm = Math.max(...Object.keys(expr)
      .filter(k => k !== 'individual' && k !== 'aux_bud')
      .map(k => expr[k].mean_tpm), 0.1);
  }}

  document.getElementById('max-tpm-label').textContent = maxTpm.toFixed(0) + ' TPM';

  // Color each tissue in the SVG based on TPM
  const tissues = ['leaf', 'roots', 'crown', 'green_fruit', 'red_fruit'];
  tissues.forEach(tissue => {{
    const svgId = TISSUE_SVG_MAP[tissue];
    if (!svgId) return;

    const tpm = expr[tissue]?.mean_tpm || 0;
    const color = tpmToColor(tpm, maxTpm);

    colorSVGGroup(svgId, color);
  }});

  // Color aux_bud inset figure (uses same maxTpm so it may exceed scale)
  const auxTpm = expr['aux_bud']?.mean_tpm || 0;
  const auxColor = tpmToColor(auxTpm, maxTpm);
  const auxShape = document.getElementById('aux-bud-shape');
  if (auxShape) {{
    auxShape.style.fill = auxColor;
    auxShape.style.transition = 'fill 0.4s ease';
  }}
  const auxLabel = document.getElementById('aux-bud-tpm');
  if (auxLabel) {{
    auxLabel.textContent = auxTpm.toFixed(1) + ' TPM';
  }}
}}


// ---- Update info panel ------------------------------------------------------
function updateInfoPanel(geneId) {{
  const meta = DATA.metadata[geneId];
  const expr = DATA.expression[geneId];
  if (!meta || !expr) return;

  document.getElementById('info-gene-name').textContent = meta.name;
  document.getElementById('info-gene-id').textContent = geneId;
  document.getElementById('info-sub_subfamily').textContent = meta.sub_subfamily;
  document.getElementById('info-source').textContent = meta.fuente_seq;
  document.getElementById('info-tmhs').textContent = meta.TMHs;
  document.getElementById('reannotation-row').style.display = meta.needs_reannotation ? '' : 'none';

  let maxTpm = currentScale === 'absolute' ?
    DATA.global_max_tpm :
    Math.max(...Object.keys(expr).filter(k => k !== 'individual' && k !== 'aux_bud').map(t => expr[t].mean_tpm), 0.1);

  const exprMap = {{}};
  Object.keys(expr).forEach(t => {{
    if (t !== 'individual') {{
      exprMap[t] = expr[t].mean_tpm;
    }}
  }});

  updateExprBars(exprMap, maxTpm);
  
  // Render individual gene expressions if available
  const indivSection = document.getElementById('individual-expression-section');
  const indivContainer = document.getElementById('individual-expression-table-container');
  
  if (expr.individual && Object.keys(expr.individual).length > 0) {{
    indivSection.style.display = 'block';
    
    let html = `<table style="width:100%; font-size:0.74rem; border-collapse:collapse; text-align:left; color:var(--text-primary); font-variant-numeric:tabular-nums;">
      <thead>
        <tr style="border-bottom: 1px solid var(--border); color:var(--text-secondary); font-size:0.7rem; text-transform:uppercase; letter-spacing:0.06em;">
          <th style="padding:5px 4px; text-align:left;">Gene</th>
          <th style="padding:5px 4px; text-align:center;">Subg</th>
          <th style="padding:5px 4px; text-align:right;">LF</th>
          <th style="padding:5px 4px; text-align:right;">CR</th>
          <th style="padding:5px 4px; text-align:right;">RT</th>
          <th style="padding:5px 4px; text-align:right;">GF</th>
          <th style="padding:5px 4px; text-align:right;">RF</th>
        </tr>
      </thead>
      <tbody>`;

    // Subgenome accent colors (TFG canonical palette)
    const SUBG_COLORS = {{ 'A': '#6D4C41', 'B': '#C2185B', 'C': '#00838F', 'D': '#455A64' }};

    // Sort by subgenome
    const genesArr = Object.keys(expr.individual).map(g => ({{id: g, ...expr.individual[g]}}));
    genesArr.sort((a,b) => a.subgenome.localeCompare(b.subgenome));

    genesArr.forEach(g => {{
      const subColor = SUBG_COLORS[g.subgenome] || 'var(--text-secondary)';
      html += `
        <tr style="border-bottom: 1px solid var(--border-soft);">
          <td style="padding:5px 4px; font-family:var(--font-mono); color:var(--text-primary); font-size:0.7rem;">${{g.id.split('g')[1] || g.id}}</td>
          <td style="padding:5px 4px; font-weight:600; text-align:center; color:${{subColor}};">${{g.subgenome}}</td>
          <td style="padding:5px 4px; text-align:right;">${{g.leaf || 0}}</td>
          <td style="padding:5px 4px; text-align:right;">${{g.crown || 0}}</td>
          <td style="padding:5px 4px; text-align:right;">${{g.roots || 0}}</td>
          <td style="padding:5px 4px; text-align:right;">${{g.green_fruit || 0}}</td>
          <td style="padding:5px 4px; text-align:right;">${{g.red_fruit || 0}}</td>
        </tr>
      `;
    }});
    html += `</tbody></table>`;
    indivContainer.innerHTML = html;
  }} else {{
    indivSection.style.display = 'none';
  }}
}}

function updateExprBars(exprMap, maxTpm) {{
  const container = document.getElementById('expr-bars');
  container.innerHTML = '';

  const tissues = ['leaf', 'green_fruit', 'red_fruit', 'crown', 'roots', 'aux_bud'];
  tissues.forEach(tissue => {{
    const tpm = exprMap[tissue] || 0;
    const pct = maxTpm > 0 ? (tpm / maxTpm * 100) : 0;
    const color = tpmToColor(tpm, maxTpm);

    const row = document.createElement('div');
    row.className = 'expr-bar-row';
    row.style.cursor = 'pointer';
    row.onclick = () => selectTissue(tissue);
    row.innerHTML = `
      <div class="expr-bar-tissue">${{TISSUE_NAMES[tissue] || tissue}}</div>
      <div class="expr-bar-wrapper">
        <div class="expr-bar" style="width:${{Math.max(pct, 1)}}%; background:${{color}}"></div>
        <div class="expr-bar-value">${{tpm.toFixed(1)}}</div>
      </div>
    `;
    container.appendChild(row);
  }});

  if (exprMap['aux_bud'] !== undefined) {{
    const warning = document.createElement('div');
    warning.className = 'n1-warning';
    warning.textContent = '⚠ Aux. Bud: N=1, descriptive only (no variance estimate)';
    container.appendChild(warning);
  }}
}}

function updateTissueDetail(geneId, tissue) {{
  const expr = DATA.expression[geneId]?.[tissue];
  const detail = document.getElementById('tissue-detail');

  if (!expr) {{
    detail.style.display = 'none';
    return;
  }}

  detail.style.display = '';
  document.getElementById('tissue-detail-title').textContent = TISSUE_NAMES[tissue] || tissue;
  document.getElementById('detail-tpm').textContent = expr.mean_tpm.toFixed(2) + ' TPM';
  document.getElementById('detail-sd').textContent = expr.sd_tpm !== null ?
    '± ' + expr.sd_tpm.toFixed(2) : 'N/A (N=1)';
  document.getElementById('detail-n').textContent = expr.n;
}}

// ---- Tooltip handlers -------------------------------------------------------
function initTooltips() {{
  const tooltip = document.getElementById('svg-tooltip');
  const ttTissue = document.getElementById('tt-tissue');
  const ttTpm = document.getElementById('tt-tpm');
  const ttDetail = document.getElementById('tt-detail');

  // Attach tooltips to the real SVG groups
  Object.entries(SVG_TISSUE_MAP).forEach(([svgId, tissue]) => {{
    const group = document.getElementById(svgId);
    if (!group) return;

    group.addEventListener('mouseenter', (e) => {{
      const tissueName = TISSUE_NAMES[tissue] || tissue;
      ttTissue.textContent = tissueName;

      if (selectedGeneId && DATA.expression[selectedGeneId]?.[tissue]) {{
        const d = DATA.expression[selectedGeneId][tissue];
        ttTpm.textContent = d.mean_tpm.toFixed(1) + ' TPM';
        const sdText = d.sd_tpm !== null ? '±' + d.sd_tpm.toFixed(1) : 'N/A';
        ttDetail.textContent = 'SD: ' + sdText + ' | N=' + d.n;
      }} else {{
        ttTpm.textContent = 'No data';
        ttDetail.textContent = 'Select a gene first';
      }}

      tooltip.classList.add('visible');
    }});

    group.addEventListener('mousemove', (e) => {{
      tooltip.style.left = (e.clientX + 16) + 'px';
      tooltip.style.top = (e.clientY - 10) + 'px';
    }});

    group.addEventListener('mouseleave', () => {{
      tooltip.classList.remove('visible');
    }});
  }});
}}

// ---- Export as PNG -----------------------------------------------------------
function exportImage() {{
  const svg = document.querySelector('#svg-container svg');
  if (!svg) return;

  const svgData = new XMLSerializer().serializeToString(svg);
  const canvas = document.createElement('canvas');
  const ctx = canvas.getContext('2d');
  const img = new Image();

  canvas.width = 1200;
  canvas.height = 1600;

  img.onload = function() {{
    // Ivory background
    ctx.fillStyle = '#F7F4EE';
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    ctx.drawImage(img, 50, 80, 1100, 1400);

    // Title (serif, dark text)
    ctx.fillStyle = '#1F2730';
    ctx.font = '600 24px "Source Serif 4", "Source Serif Pro", Georgia, serif';
    const geneLabel = selectedGeneId ?
      (DATA.metadata[selectedGeneId]?.name || selectedGeneId) : 'No gene selected';
    ctx.fillText('Fragaria x ananassa  —  ' + geneLabel, 20, 44);

    // Color scale legend (warm beige → dark red)
    ctx.font = '14px Inter, sans-serif';
    ctx.fillStyle = '#4A4842';
    ctx.fillText('Expression (TPM)', 20, 1540);
    const gradient = ctx.createLinearGradient(180, 1530, 500, 1530);
    gradient.addColorStop(0.00, '#C9C3B2');
    gradient.addColorStop(0.03, '#C9C3B2');
    gradient.addColorStop(0.05, '#FFF5C2');
    gradient.addColorStop(0.22, '#FDB863');
    gradient.addColorStop(0.50, '#E5601F');
    gradient.addColorStop(0.80, '#8E1414');
    gradient.addColorStop(1.00, '#5A0A0A');
    ctx.fillStyle = gradient;
    ctx.fillRect(180, 1520, 320, 16);
    ctx.strokeStyle = '#D6CEC1';
    ctx.lineWidth = 1;
    ctx.strokeRect(180, 1520, 320, 16);
    ctx.fillStyle = '#4A4842';
    ctx.fillText('0', 165, 1535);
    ctx.fillText(DATA.global_max_tpm.toFixed(0), 510, 1535);

    const link = document.createElement('a');
    link.download = 'efp_' + (selectedGeneId || 'plant') + '.png';
    link.href = canvas.toDataURL('image/png');
    link.click();
  }};

  img.src = 'data:image/svg+xml;base64,' + btoa(unescape(encodeURIComponent(svgData)));
}}

// ---- Start ------------------------------------------------------------------
window.addEventListener('DOMContentLoaded', () => {{
  init();
  initTooltips();
}});
</script>
</body>
</html>"""

# ---- Write output ------------------------------------------------------------
os.makedirs(os.path.dirname(OUTPUT_HTML), exist_ok=True)
with open(OUTPUT_HTML, "w", encoding="utf-8") as f:
    f.write(html_content)

print(f"# Output: {OUTPUT_HTML}")
print(f"# Open in browser: file://{os.path.abspath(OUTPUT_HTML)}")
print("# === 15_homeolog_efp_viewer.py complete ===")
