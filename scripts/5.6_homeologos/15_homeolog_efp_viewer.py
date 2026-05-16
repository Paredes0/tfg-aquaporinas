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
<title>Fragaria × ananassa Aquaporin eFP Viewer</title>
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap" rel="stylesheet">
<style>
  :root {{
    --bg-primary: #0f1117;
    --bg-secondary: #1a1d27;
    --bg-card: #222536;
    --bg-hover: #2a2e42;
    --text-primary: #e8eaed;
    --text-secondary: #9aa0b0;
    --text-muted: #6b7185;
    --accent-primary: #6c63ff;
    --accent-secondary: #4ecdc4;
    --accent-warm: #ff6b6b;
    --border: #2d3147;
    --border-active: #6c63ff;
    --shadow: 0 4px 24px rgba(0,0,0,0.3);
    --radius: 12px;
    --radius-sm: 8px;
  }}

  * {{ margin: 0; padding: 0; box-sizing: border-box; }}

  body {{
    font-family: 'Inter', -apple-system, sans-serif;
    background: var(--bg-primary);
    color: var(--text-primary);
    min-height: 100vh;
    overflow-x: hidden;
  }}

  .app-header {{
    background: linear-gradient(135deg, var(--bg-secondary), var(--bg-card));
    border-bottom: 1px solid var(--border);
    padding: 20px 32px;
    display: flex;
    align-items: center;
    gap: 20px;
    box-shadow: var(--shadow);
  }}

  .app-header h1 {{
    font-size: 1.4rem;
    font-weight: 600;
    background: linear-gradient(135deg, var(--accent-secondary), var(--accent-primary));
    -webkit-background-clip: text;
    -webkit-text-fill-color: transparent;
    background-clip: text;
  }}

  .app-header .subtitle {{
    color: var(--text-secondary);
    font-size: 0.85rem;
    font-weight: 300;
  }}

  .logo-icon {{
    width: 40px; height: 40px;
    background: linear-gradient(135deg, #ff6b6b, #ff8e53);
    border-radius: 10px;
    display: flex; align-items: center; justify-content: center;
    font-size: 1.5rem;
  }}

  .main-container {{
    display: grid;
    grid-template-columns: 320px 1fr 340px;
    gap: 0;
    height: calc(100vh - 80px);
  }}

  .panel {{
    background: var(--bg-secondary);
    border-right: 1px solid var(--border);
    overflow-y: auto;
    padding: 20px;
  }}

  .panel:last-child {{
    border-right: none;
    border-left: 1px solid var(--border);
  }}

  .panel-title {{
    font-size: 0.75rem;
    font-weight: 600;
    text-transform: uppercase;
    letter-spacing: 1.2px;
    color: var(--text-muted);
    margin-bottom: 16px;
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
    background: var(--bg-card);
    border: 1px solid var(--border);
    color: var(--text-primary);
    padding: 10px 12px;
    border-radius: var(--radius-sm);
    font-family: inherit;
    font-size: 0.85rem;
    outline: none;
    transition: border-color 0.2s;
    appearance: none;
    -webkit-appearance: none;
    cursor: pointer;
  }}

  select {{
    background-image: url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' width='12' height='12' viewBox='0 0 12 12'%3E%3Cpath fill='%239aa0b0' d='M6 8L1 3h10z'/%3E%3C/svg%3E");
    background-repeat: no-repeat;
    background-position: right 12px center;
    padding-right: 32px;
  }}

  select:focus, select:hover {{
    border-color: var(--accent-primary);
  }}

  .gene-list {{
    max-height: 400px;
    overflow-y: auto;
    border: 1px solid var(--border);
    border-radius: var(--radius-sm);
  }}

  .gene-item {{
    padding: 8px 12px;
    cursor: pointer;
    font-size: 0.8rem;
    border-bottom: 1px solid var(--border);
    display: flex;
    justify-content: space-between;
    align-items: center;
    transition: background 0.15s;
  }}

  .gene-item:hover {{
    background: var(--bg-hover);
  }}

  .gene-item.active {{
    background: rgba(108, 99, 255, 0.15);
    border-left: 3px solid var(--accent-primary);
  }}

  .gene-item .gene-name {{
    font-weight: 500;
    color: var(--text-primary);
  }}

  .gene-item .gene-id {{
    font-size: 0.7rem;
    color: var(--text-muted);
    font-family: 'Courier New', monospace;
  }}

  .gene-item .reannotation-badge {{
    background: var(--accent-warm);
    color: white;
    font-size: 0.6rem;
    padding: 2px 6px;
    border-radius: 4px;
    font-weight: 600;
  }}

  .plant-viewer {{
    display: flex;
    align-items: center;
    justify-content: center;
    background: radial-gradient(ellipse at center, #1a2030 0%, var(--bg-primary) 100%);
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
    background: var(--bg-card);
    border-radius: var(--radius-sm);
    padding: 16px;
    margin-bottom: 16px;
    border: 1px solid var(--border);
  }}

  .info-section h3 {{
    font-size: 0.85rem;
    font-weight: 600;
    color: var(--accent-secondary);
    margin-bottom: 10px;
  }}

  .info-row {{
    display: flex;
    justify-content: space-between;
    padding: 5px 0;
    font-size: 0.8rem;
    border-bottom: 1px solid rgba(255,255,255,0.05);
  }}

  .info-row:last-child {{
    border-bottom: none;
  }}

  .info-row .label {{
    color: var(--text-muted);
  }}

  .info-row .value {{
    color: var(--text-primary);
    font-weight: 500;
    text-align: right;
  }}

  /* Color scale legend — TAIR10 style */
  .color-scale {{
    display: flex;
    align-items: center;
    gap: 8px;
    padding: 12px 16px;
    background: var(--bg-card);
    border-radius: var(--radius-sm);
    margin-bottom: 16px;
    border: 1px solid var(--border);
  }}

  .color-scale .scale-bar {{
    flex: 1;
    height: 14px;
    border-radius: 7px;
    background: linear-gradient(90deg, #1a1a2e, #4a3a20, #8a7a30, #ccaa20, #ffcc00, #ff8800, #dd2200, #aa0000);
  }}

  .color-scale .scale-label {{
    font-size: 0.7rem;
    color: var(--text-muted);
    white-space: nowrap;
  }}

  .mode-toggle {{
    display: flex;
    background: var(--bg-card);
    border-radius: var(--radius-sm);
    overflow: hidden;
    margin-bottom: 16px;
    border: 1px solid var(--border);
  }}

  .mode-btn {{
    flex: 1;
    padding: 10px;
    text-align: center;
    font-size: 0.78rem;
    font-weight: 500;
    cursor: pointer;
    transition: all 0.2s;
    color: var(--text-muted);
    background: transparent;
    border: none;
    font-family: inherit;
  }}

  .mode-btn.active {{
    background: var(--accent-primary);
    color: white;
  }}

  .mode-btn:hover:not(.active) {{
    background: var(--bg-hover);
    color: var(--text-primary);
  }}

  /* Expression bar chart in info panel */
  .expr-bar-container {{
    margin-top: 8px;
  }}

  .expr-bar-row {{
    display: flex;
    align-items: center;
    margin-bottom: 6px;
    gap: 8px;
  }}

  .expr-bar-tissue {{
    width: 80px;
    font-size: 0.72rem;
    color: var(--text-secondary);
    text-align: right;
  }}

  .expr-bar-wrapper {{
    flex: 1;
    height: 18px;
    background: var(--bg-primary);
    border-radius: 4px;
    overflow: hidden;
    position: relative;
  }}

  .expr-bar {{
    height: 100%;
    border-radius: 4px;
    transition: width 0.5s ease;
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
    text-shadow: 0 0 4px rgba(0,0,0,0.8);
  }}

  .n1-warning {{
    font-size: 0.68rem;
    color: var(--accent-warm);
    font-style: italic;
    margin-top: 4px;
  }}

  /* Search */
  .search-input {{
    width: 100%;
    background: var(--bg-card);
    border: 1px solid var(--border);
    color: var(--text-primary);
    padding: 10px 12px 10px 34px;
    border-radius: var(--radius-sm);
    font-family: inherit;
    font-size: 0.82rem;
    outline: none;
    transition: border-color 0.2s;
    margin-bottom: 12px;
    background-image: url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' width='16' height='16' viewBox='0 0 24 24' fill='none' stroke='%236b7185' stroke-width='2'%3E%3Ccircle cx='11' cy='11' r='8'/%3E%3Cline x1='21' y1='21' x2='16.65' y2='16.65'/%3E%3C/svg%3E");
    background-repeat: no-repeat;
    background-position: 10px center;
  }}

  .search-input:focus {{
    border-color: var(--accent-primary);
  }}

  .no-data-msg {{
    text-align: center;
    padding: 40px 20px;
    color: var(--text-muted);
    font-size: 0.85rem;
  }}

  .scroll-custom::-webkit-scrollbar {{
    width: 6px;
  }}
  .scroll-custom::-webkit-scrollbar-track {{
    background: transparent;
  }}
  .scroll-custom::-webkit-scrollbar-thumb {{
    background: var(--border);
    border-radius: 3px;
  }}

  /* Floating tooltip */
  .svg-tooltip {{
    position: fixed;
    pointer-events: none;
    background: rgba(20, 22, 35, 0.95);
    border: 1px solid var(--accent-secondary);
    border-radius: 8px;
    padding: 10px 14px;
    font-size: 0.78rem;
    color: var(--text-primary);
    box-shadow: 0 8px 32px rgba(0,0,0,0.5);
    z-index: 1000;
    opacity: 0;
    transition: opacity 0.15s;
    max-width: 220px;
    backdrop-filter: blur(8px);
  }}
  .svg-tooltip.visible {{
    opacity: 1;
  }}
  .svg-tooltip .tt-tissue {{
    font-weight: 600;
    color: var(--accent-secondary);
    margin-bottom: 4px;
  }}
  .svg-tooltip .tt-tpm {{
    font-size: 1.1rem;
    font-weight: 700;
    color: #ffd93d;
  }}
  .svg-tooltip .tt-detail {{
    font-size: 0.7rem;
    color: var(--text-muted);
    margin-top: 2px;
  }}

  /* Export button */
  .export-btn {{
    position: absolute;
    top: 12px;
    right: 12px;
    background: var(--bg-card);
    border: 1px solid var(--border);
    color: var(--text-secondary);
    padding: 6px 14px;
    border-radius: 6px;
    font-size: 0.75rem;
    cursor: pointer;
    font-family: inherit;
    transition: all 0.2s;
    z-index: 10;
  }}
  .export-btn:hover {{
    background: var(--accent-primary);
    color: white;
    border-color: var(--accent-primary);
  }}

  /* Aux bud inset figure */
  .aux-bud-inset {{
    position: absolute;
    bottom: 20px;
    right: 20px;
    background: rgba(20, 22, 35, 0.85);
    border: 1px solid var(--border);
    border-radius: var(--radius-sm);
    padding: 12px 16px;
    backdrop-filter: blur(6px);
    z-index: 5;
    cursor: pointer;
    transition: border-color 0.2s, box-shadow 0.2s;
    text-align: center;
  }}
  .aux-bud-inset:hover {{
    border-color: var(--accent-secondary);
    box-shadow: 0 0 12px rgba(78, 205, 196, 0.2);
  }}
  .aux-bud-inset .aux-label {{
    font-size: 0.72rem;
    color: var(--text-secondary);
    margin-bottom: 6px;
    font-weight: 600;
  }}
  .aux-bud-inset .aux-tpm {{
    font-size: 0.8rem;
    color: var(--text-primary);
    font-weight: 700;
    margin-top: 4px;
  }}
  .aux-bud-inset .aux-warning {{
    font-size: 0.6rem;
    color: var(--accent-warm);
    font-style: italic;
    margin-top: 2px;
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
  <div class="logo-icon">🍓</div>
  <div>
    <h1>Fragaria × ananassa — Aquaporin eFP Viewer</h1>
    <div class="subtitle">Interactive tissue expression atlas • Basal expression (TPM) • TAIR10-style coloring</div>
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

    <div class="color-scale" style="margin-top: 20px;">
      <span class="scale-label">0</span>
      <div class="scale-bar"></div>
      <span class="scale-label" id="max-tpm-label">TPM</span>
    </div>

    <div class="info-section" style="margin-top: 8px;">
      <h3>📊 Scale Mode</h3>
      <div class="mode-toggle">
        <button class="mode-btn active" onclick="setScale('absolute')" id="scale-abs">Absolute</button>
        <button class="mode-btn" onclick="setScale('relative')" id="scale-rel">Relative</button>
      </div>
      <div style="font-size:0.72rem; color:var(--text-muted); margin-top:8px;">
        <b>Absolute:</b> Color scaled to global max TPM<br>
        <b>Relative:</b> Color scaled to selected gene's max
      </div>
    </div>
  </div>

  <!-- CENTER: Plant SVG -->
  <div class="plant-viewer" id="plant-viewer">
    <button class="export-btn" onclick="exportImage()">📷 Export PNG</button>
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
        <ellipse id="aux-bud-shape" cx="38" cy="30" rx="14" ry="10" fill="#555"/>
        <!-- Bud leaves -->
        <path d="M32,25 Q38,15 44,25" fill="#6B9E3A" opacity="0.8"/>
        <path d="M30,27 Q36,18 42,26" fill="#7CB342" opacity="0.7"/>
        <!-- Small emerging leaves -->
        <ellipse cx="48" cy="24" rx="6" ry="3" fill="#8BC34A" transform="rotate(-20 48 24)" opacity="0.6"/>
      </svg>
      <div class="aux-tpm" id="aux-bud-tpm">— TPM</div>
      <div class="aux-warning">⚠ N=1 (descriptive)</div>
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
        <h3>🧬 Gene Info</h3>
        <div class="info-row"><span class="label">Name</span><span class="value" id="info-gene-name">-</span></div>
        <div class="info-row"><span class="label">Group ID</span><span class="value" id="info-gene-id">-</span></div>
        <div class="info-row"><span class="label">Sub-subfamily</span><span class="value" id="info-sub_subfamily">-</span></div>
        <div class="info-row"><span class="label">Source</span><span class="value" id="info-source">-</span></div>
        <div class="info-row"><span class="label">TMHs (est.)</span><span class="value" id="info-tmhs">-</span></div>
        <div class="info-row" id="reannotation-row" style="display:none;">
          <span class="label">⚠️ Reannotation</span>
          <span class="value" style="color:var(--accent-warm)">Needed</span>
        </div>
      </div>

      <div class="info-section">
        <h3>📊 Expression by Tissue (TPM)</h3>
        <div id="expr-bars" class="expr-bar-container"></div>
      </div>

      <div class="info-section" id="tissue-detail" style="display:none;">
        <h3 id="tissue-detail-title">Tissue Detail</h3>
        <div class="info-row"><span class="label">Mean TPM</span><span class="value" id="detail-tpm">-</span></div>
        <div class="info-row"><span class="label">SD</span><span class="value" id="detail-sd">-</span></div>
        <div class="info-row"><span class="label">Replicates</span><span class="value" id="detail-n">-</span></div>
      </div>

      <div class="info-section" id="individual-expression-section" style="display:none;">
        <h3>🔬 Individual Homeolog Expression</h3>
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
let currentScale = 'absolute';
let selectedGeneId = null;
let selectedTissue = null;

// ---- TAIR10-style color interpolation ---------------------------------------
// TAIR eFP Browser uses a yellow (low) → red (high) scale
// 0 expression = light yellow (#FFFF00)
// max expression = dark red (#CC0000)
function tpmToColor(tpm, maxTpm) {{
  if (maxTpm <= 0) maxTpm = 1;
  const ratio = Math.min(tpm / maxTpm, 1);
  // Progressive: dark bg → faint yellow → yellow → orange → red (like opacity ramp)
  const stops = [
    [0.0,  [ 26,  26,  46]],   // background — no expression
    [0.05, [ 50,  42,  30]],   // barely visible warm hint
    [0.15, [ 90,  75,  25]],   // faint dark yellow
    [0.3,  [160, 140,  20]],   // muted yellow
    [0.45, [210, 190,  15]],   // yellow
    [0.6,  [255, 200,   0]],   // bright yellow
    [0.75, [255, 140,   0]],   // orange
    [0.9,  [220,  40,   0]],   // red
    [1.0,  [170,   0,   0]]    // dark red — max expression
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

  // Set up outline (path1) to be non-interactive
  const outline = svg.querySelector('#path1');
  if (outline) {{
    outline.style.pointerEvents = 'none';
    outline.style.fill = 'none';
    outline.style.stroke = 'rgba(255,255,255,0.25)';
    outline.style.strokeWidth = '0.5';
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

    // Hover effects
    group.addEventListener('mouseenter', function() {{
      this.style.filter = 'brightness(1.3) drop-shadow(0 0 8px rgba(255,255,255,0.3))';
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

  // Initial coloring — matches dark background (no expression)
  Object.entries(SVG_TISSUE_MAP).forEach(([svgId, tissue]) => {{
    colorSVGGroup(svgId, 'rgb(26, 26, 46)');
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
  document.getElementById('no-selection').style.display = 'none';
  document.getElementById('gene-info').style.display = 'block';
}}

function selectTissue(tissue) {{
  selectedTissue = tissue;

  // Update visual selection: reset all, highlight selected
  Object.entries(SVG_TISSUE_MAP).forEach(([svgId, t]) => {{
    const group = document.getElementById(svgId);
    if (!group) return;
    if (t === tissue) {{
      group.style.filter = 'brightness(1.3) drop-shadow(0 0 12px rgba(78,205,196,0.5))';
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
    
    let html = `<table style="width:100%; font-size:0.75rem; border-collapse:collapse; text-align:left; color:var(--text-primary);">
      <thead>
        <tr style="border-bottom: 1px solid var(--border); color:var(--text-muted);">
          <th style="padding:4px;">Gene</th>
          <th style="padding:4px;">Subg</th>
          <th style="padding:4px;">LF</th>
          <th style="padding:4px;">CR</th>
          <th style="padding:4px;">RT</th>
          <th style="padding:4px;">GF</th>
          <th style="padding:4px;">RF</th>
        </tr>
      </thead>
      <tbody>`;
      
    // Sort by subgenome
    const genesArr = Object.keys(expr.individual).map(g => ({{id: g, ...expr.individual[g]}}));
    genesArr.sort((a,b) => a.subgenome.localeCompare(b.subgenome));
    
    genesArr.forEach(g => {{
      html += `
        <tr style="border-bottom: 1px dashed #333;">
          <td style="padding:4px; font-family:monospace; color:var(--accent-secondary);">${{g.id.split('g')[1] || g.id}}</td>
          <td style="padding:4px; font-weight:bold;">${{g.subgenome}}</td>
          <td style="padding:4px;">${{g.leaf || 0}}</td>
          <td style="padding:4px;">${{g.crown || 0}}</td>
          <td style="padding:4px;">${{g.roots || 0}}</td>
          <td style="padding:4px;">${{g.green_fruit || 0}}</td>
          <td style="padding:4px;">${{g.red_fruit || 0}}</td>
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
    // Dark background
    ctx.fillStyle = '#0f1117';
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    ctx.drawImage(img, 50, 80, 1100, 1400);

    // Add title
    ctx.fillStyle = '#e8eaed';
    ctx.font = 'bold 24px Inter, sans-serif';
    const geneLabel = selectedGeneId ?
      (DATA.metadata[selectedGeneId]?.name || selectedGeneId) : 'No gene selected';
    ctx.fillText('Fragaria × ananassa — ' + geneLabel, 20, 40);

    // Add color scale legend
    ctx.font = '14px Inter, sans-serif';
    ctx.fillStyle = '#9aa0b0';
    ctx.fillText('Expression (TPM)', 20, 1540);
    const gradient = ctx.createLinearGradient(180, 1530, 500, 1530);
    gradient.addColorStop(0, '#FFFF00');
    gradient.addColorStop(0.5, '#FF6600');
    gradient.addColorStop(1, '#CC0000');
    ctx.fillStyle = gradient;
    ctx.fillRect(180, 1520, 320, 16);
    ctx.fillStyle = '#9aa0b0';
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
