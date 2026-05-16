#!/usr/bin/env python3
"""
10_generate_efp_viewer.py — Generate interactive eFP-like expression viewer

Creates a standalone HTML file with an SVG strawberry plant illustration
where tissues are colored by aquaporin expression levels (TPM).

Usage:
    python3 10_generate_efp_viewer.py

Reads from:
    results/basal_aquaporins/basal_aquaporins_summary.csv
    tabla_Aquaporinas_traduccion.tabular

Outputs:
    results/basal_aquaporins/efp_viewer.html
"""

import pandas as pd
import json
import os
import sys
import math

# ---- Configuration ----------------------------------------------------------
PROJECT_DIR  = os.environ.get("PROJECT_DIR", os.getcwd())
SUMMARY_FILE = os.path.join(PROJECT_DIR, "results", "basal_aquaporins", "basal_aquaporins_summary.csv")
AQP_TABLE    = os.path.join(PROJECT_DIR, "tabla_Aquaporinas_traduccion.tabular")
OUTPUT_HTML  = os.path.join(PROJECT_DIR, "results", "basal_aquaporins", "efp_viewer.html")

print("# =============================================================")
print("# GENERATING eFP-like EXPRESSION VIEWER")
print("# =============================================================")

# ---- Load data or use demo --------------------------------------------------
if os.path.exists(SUMMARY_FILE):
    summary = pd.read_csv(SUMMARY_FILE)
    print(f"# Loaded expression summary: {len(summary)} entries")
else:
    print(f"# WARNING: {SUMMARY_FILE} not found. Using demo data.")
    # Generate minimal demo data
    tissues = ["green_fruit", "red_fruit", "crown", "leaf", "roots", "aux_bud"]
    genes = [f"DemoGene_{i}" for i in range(1, 11)]
    subfamilies = ["PIP", "TIP", "NIP", "SIP", "XIP"]
    rows = []
    import random
    random.seed(42)
    for g in genes:
        sf = random.choice(subfamilies)
        for t in tissues:
            rows.append({
                "gene_id": g,
                "tissue": t,
                "mean_tpm": random.uniform(0, 200),
                "sd_tpm": random.uniform(0, 30),
                "n": 3 if t != "aux_bud" else 1,
                "aqp_family_subfamily": f"Fa{sf}{random.randint(1,5)}",
                "subfamilia_phylo": sf,
                "fuente_seq": "GFF3",
                "TMHs": 6,
                "needs_reannotation": False,
                "note": "N=1_descriptive_only" if t == "aux_bud" else ""
            })
    summary = pd.DataFrame(rows)

if os.path.exists(AQP_TABLE):
    aqp_info = pd.read_csv(AQP_TABLE, sep="\t")
    print(f"# Loaded aquaporin table: {len(aqp_info)} genes")
else:
    aqp_info = None
    print("# WARNING: Aquaporin table not found. Using summary data only.")

# ---- Prepare JSON data for embedding ----------------------------------------
# Build expression matrix: {gene_id: {tissue: {mean_tpm, sd_tpm, n}}}
expression_data = {}
gene_metadata = {}

for _, row in summary.iterrows():
    gid = row["gene_id"]
    tissue = row["tissue"]

    if gid not in expression_data:
        expression_data[gid] = {}
        gene_metadata[gid] = {
            "name": row.get("aqp_family_subfamily", gid),
            "subfamily": row.get("subfamilia_phylo", "Unknown"),
            "fuente_seq": row.get("fuente_seq", "N/A"),
            "TMHs": int(row["TMHs"]) if pd.notna(row.get("TMHs")) else "N/A",
            "needs_reannotation": bool(row.get("needs_reannotation", False))
        }

    expression_data[gid][tissue] = {
        "mean_tpm": round(float(row["mean_tpm"]), 2) if pd.notna(row["mean_tpm"]) else 0,
        "sd_tpm": round(float(row["sd_tpm"]), 2) if pd.notna(row.get("sd_tpm")) else None,
        "n": int(row["n"]) if pd.notna(row.get("n")) else 1
    }

# Build subfamily list
subfamilies = sorted(set(m["subfamily"] for m in gene_metadata.values()))
# Build gene list per subfamily
genes_by_subfamily = {}
for gid, meta in gene_metadata.items():
    sf = meta["subfamily"]
    if sf not in genes_by_subfamily:
        genes_by_subfamily[sf] = []
    genes_by_subfamily[sf].append({"id": gid, "name": meta["name"]})

# Sort genes within each subfamily
for sf in genes_by_subfamily:
    genes_by_subfamily[sf].sort(key=lambda x: x["name"])

# Calculate global max TPM for scale
all_tpms = []
for gid_data in expression_data.values():
    for tdata in gid_data.values():
        all_tpms.append(tdata["mean_tpm"])
global_max_tpm = max(all_tpms) if all_tpms else 100

data_json = json.dumps({
    "expression": expression_data,
    "metadata": gene_metadata,
    "subfamilies": subfamilies,
    "genes_by_subfamily": genes_by_subfamily,
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
    align-items: center;
    justify-content: center;
  }}

  .plant-svg-container svg {{
    max-width: 90%;
    max-height: 90%;
  }}

  /* Tissue SVG parts */
  .tissue-part {{
    cursor: pointer;
    transition: opacity 0.3s, filter 0.3s;
    stroke: rgba(255,255,255,0.15);
    stroke-width: 1;
  }}

  .tissue-part:hover {{
    filter: brightness(1.3) drop-shadow(0 0 8px rgba(255,255,255,0.3));
    stroke: rgba(255,255,255,0.5);
    stroke-width: 2;
  }}

  .tissue-part.selected {{
    stroke: var(--accent-secondary);
    stroke-width: 2.5;
    filter: drop-shadow(0 0 12px rgba(78,205,196,0.4));
  }}

  .tissue-label {{
    font-family: 'Inter', sans-serif;
    font-size: 11px;
    fill: var(--text-secondary);
    pointer-events: none;
    font-weight: 500;
  }}

  .tissue-tpm-label {{
    font-family: 'Inter', sans-serif;
    font-size: 10px;
    fill: white;
    pointer-events: none;
    font-weight: 600;
    text-anchor: middle;
    paint-order: stroke;
    stroke: rgba(0,0,0,0.7);
    stroke-width: 3px;
    stroke-linecap: round;
    stroke-linejoin: round;
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

  /* Color scale legend */
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
    background: linear-gradient(90deg, #1a1a2e, #16213e, #0f3460, #e94560, #ff6b6b, #ffd93d);
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

  .soil-line {{
    stroke: #5a4a3a;
    stroke-width: 2;
    stroke-dasharray: 8 4;
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
    <div class="subtitle">Interactive tissue expression atlas • Basal expression (TPM)</div>
  </div>
</header>

<div class="main-container">
  <!-- LEFT PANEL: Controls -->
  <div class="panel scroll-custom">
    <div class="panel-title">Gene Selection</div>

    <div class="mode-toggle">
      <button class="mode-btn active" onclick="setMode('single')" id="mode-single">Single Gene</button>
      <button class="mode-btn" onclick="setMode('subfamily')" id="mode-subfamily">Subfamily Mean</button>
    </div>

    <div class="control-group">
      <label>Subfamily</label>
      <select id="subfamily-select" onchange="onSubfamilyChange()">
        <option value="ALL">All subfamilies</option>
      </select>
    </div>

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
    <div class="plant-svg-container">
      <svg viewBox="0 0 600 700" xmlns="http://www.w3.org/2000/svg" id="plant-svg">
        <!-- Background gradient -->
        <defs>
          <linearGradient id="sky" x1="0" y1="0" x2="0" y2="1">
            <stop offset="0%" stop-color="#1a2030"/>
            <stop offset="60%" stop-color="#1a2530"/>
            <stop offset="100%" stop-color="#2a2520"/>
          </linearGradient>
          <linearGradient id="soilGrad" x1="0" y1="0" x2="0" y2="1">
            <stop offset="0%" stop-color="#3a3025"/>
            <stop offset="100%" stop-color="#2a2015"/>
          </linearGradient>
          <!-- Leaf shape -->
          <filter id="glow">
            <feGaussianBlur stdDeviation="3" result="coloredBlur"/>
            <feMerge>
              <feMergeNode in="coloredBlur"/>
              <feMergeNode in="SourceGraphic"/>
            </feMerge>
          </filter>
        </defs>

        <!-- Soil background -->
        <rect x="0" y="430" width="600" height="270" fill="url(#soilGrad)" opacity="0.5"/>
        <line x1="30" y1="430" x2="570" y2="430" class="soil-line" opacity="0.6"/>

        <!-- ============ ROOTS ============ -->
        <g id="tissue-roots" class="tissue-part" onclick="selectTissue('roots')" data-tissue="roots">
          <!-- Main roots -->
          <path d="M300,410 Q280,470 250,530 Q240,560 220,600 Q210,620 195,650" fill="none" stroke="#8B7355" stroke-width="6" stroke-linecap="round"/>
          <path d="M300,410 Q310,460 300,520 Q295,560 290,610 Q288,630 285,660" fill="none" stroke="#8B7355" stroke-width="5" stroke-linecap="round"/>
          <path d="M300,410 Q330,470 360,530 Q375,560 390,600 Q398,625 405,650" fill="none" stroke="#8B7355" stroke-width="6" stroke-linecap="round"/>
          <!-- Branch roots -->
          <path d="M270,490 Q245,510 225,540" fill="none" stroke="#8B7355" stroke-width="3" stroke-linecap="round"/>
          <path d="M260,530 Q235,545 215,570" fill="none" stroke="#8B7355" stroke-width="2.5" stroke-linecap="round"/>
          <path d="M340,500 Q365,520 385,545" fill="none" stroke="#8B7355" stroke-width="3" stroke-linecap="round"/>
          <path d="M350,540 Q375,555 395,575" fill="none" stroke="#8B7355" stroke-width="2.5" stroke-linecap="round"/>
          <path d="M295,550 Q275,575 260,600" fill="none" stroke="#8B7355" stroke-width="2" stroke-linecap="round"/>
          <path d="M305,570 Q325,590 345,615" fill="none" stroke="#8B7355" stroke-width="2" stroke-linecap="round"/>
          <!-- Fine roots -->
          <path d="M230,560 Q215,575 200,595" fill="none" stroke="#8B7355" stroke-width="1.5" stroke-linecap="round"/>
          <path d="M370,560 Q385,575 400,590" fill="none" stroke="#8B7355" stroke-width="1.5" stroke-linecap="round"/>
          <path d="M290,620 Q275,635 260,655" fill="none" stroke="#8B7355" stroke-width="1.5" stroke-linecap="round"/>
          <path d="M310,610 Q325,630 340,650" fill="none" stroke="#8B7355" stroke-width="1.5" stroke-linecap="round"/>
        </g>
        <text x="300" y="680" class="tissue-label" text-anchor="middle">Roots</text>

        <!-- ============ CROWN ============ -->
        <g id="tissue-crown" class="tissue-part" onclick="selectTissue('crown')" data-tissue="crown">
          <ellipse cx="300" cy="395" rx="50" ry="30" fill="#B8956A"/>
          <ellipse cx="300" cy="390" rx="45" ry="22" fill="#C9A87C"/>
          <path d="M260,390 Q270,370 280,385" fill="#A08060" stroke="none"/>
          <path d="M320,390 Q330,370 340,385" fill="#A08060" stroke="none"/>
        </g>
        <text x="230" y="410" class="tissue-label" text-anchor="middle">Crown</text>

        <!-- ============ LEAVES (left cluster) ============ -->
        <g id="tissue-leaf" class="tissue-part" onclick="selectTissue('leaf')" data-tissue="leaf">
          <!-- Stem left -->
          <path d="M285,375 Q240,310 180,260" fill="none" stroke="#5A8C3A" stroke-width="4" stroke-linecap="round"/>
          <!-- Trifoliate leaf left -->
          <ellipse cx="155" cy="240" rx="45" ry="28" fill="#3D8B37" transform="rotate(-25 155 240)"/>
          <ellipse cx="185" cy="215" rx="40" ry="25" fill="#4A9E43" transform="rotate(-5 185 215)"/>
          <ellipse cx="145" cy="215" rx="38" ry="24" fill="#3D8B37" transform="rotate(-40 145 215)"/>
          <!-- Leaf veins -->
          <path d="M155,240 L130,225" fill="none" stroke="#2D6B27" stroke-width="1" opacity="0.5"/>
          <path d="M185,215 L205,200" fill="none" stroke="#2D6B27" stroke-width="1" opacity="0.5"/>

          <!-- Stem center-left -->
          <path d="M295,370 Q270,280 260,210" fill="none" stroke="#5A8C3A" stroke-width="4" stroke-linecap="round"/>
          <!-- Trifoliate leaf top -->
          <ellipse cx="258" cy="185" rx="42" ry="26" fill="#4A9E43" transform="rotate(-10 258 185)"/>
          <ellipse cx="280" cy="162" rx="38" ry="24" fill="#3D8B37" transform="rotate(10 280 162)"/>
          <ellipse cx="238" cy="168" rx="36" ry="22" fill="#4A9E43" transform="rotate(-30 238 168)"/>

          <!-- Stem right -->
          <path d="M310,372 Q350,300 390,255" fill="none" stroke="#5A8C3A" stroke-width="4" stroke-linecap="round"/>
          <!-- Trifoliate leaf right -->
          <ellipse cx="410" cy="240" rx="43" ry="27" fill="#3D8B37" transform="rotate(20 410 240)"/>
          <ellipse cx="385" cy="218" rx="38" ry="24" fill="#4A9E43" transform="rotate(5 385 218)"/>
          <ellipse cx="420" cy="218" rx="36" ry="22" fill="#3D8B37" transform="rotate(35 420 218)"/>
        </g>
        <text x="260" y="145" class="tissue-label" text-anchor="middle">Leaves</text>

        <!-- ============ GREEN FRUIT ============ -->
        <g id="tissue-green_fruit" class="tissue-part" onclick="selectTissue('green_fruit')" data-tissue="green_fruit">
          <!-- Stem -->
          <path d="M265,380 Q230,350 200,340" fill="none" stroke="#5A8C3A" stroke-width="3" stroke-linecap="round"/>
          <!-- Fruit -->
          <ellipse cx="185" cy="348" rx="22" ry="26" fill="#7CB342" transform="rotate(10 185 348)"/>
          <!-- Seeds -->
          <circle cx="178" cy="338" r="1.5" fill="#9CCC65"/>
          <circle cx="190" cy="342" r="1.5" fill="#9CCC65"/>
          <circle cx="182" cy="355" r="1.5" fill="#9CCC65"/>
          <circle cx="193" cy="355" r="1.5" fill="#9CCC65"/>
          <!-- Calyx -->
          <path d="M183,323 L178,315 L185,322" fill="#4A7C25"/>
          <path d="M188,322 L193,313 L190,322" fill="#4A7C25"/>
          <path d="M180,324 L172,318 L180,326" fill="#3D6B20"/>
        </g>
        <text x="155" y="310" class="tissue-label" text-anchor="middle">Green Fruit</text>

        <!-- ============ RED FRUIT ============ -->
        <g id="tissue-red_fruit" class="tissue-part" onclick="selectTissue('red_fruit')" data-tissue="red_fruit">
          <!-- Stem -->
          <path d="M335,375 Q380,340 420,325" fill="none" stroke="#5A8C3A" stroke-width="3" stroke-linecap="round"/>
          <!-- Fruit -->
          <path d="M420,320 Q440,330 445,355 Q448,375 435,390 Q420,405 405,395 Q390,385 392,360 Q393,335 420,320Z" fill="#DC143C"/>
          <!-- Seeds -->
          <circle cx="415" cy="345" r="1.8" fill="#FFD700"/>
          <circle cx="428" cy="350" r="1.8" fill="#FFD700"/>
          <circle cx="410" cy="362" r="1.8" fill="#FFD700"/>
          <circle cx="425" cy="368" r="1.8" fill="#FFD700"/>
          <circle cx="435" cy="358" r="1.8" fill="#FFD700"/>
          <circle cx="415" cy="378" r="1.8" fill="#FFD700"/>
          <circle cx="430" cy="380" r="1.8" fill="#FFD700"/>
          <!-- Calyx -->
          <path d="M418,320 L412,308 L420,318" fill="#4A7C25"/>
          <path d="M423,318 L430,306 L425,317" fill="#4A7C25"/>
          <path d="M415,322 L405,312 L416,322" fill="#3D6B20"/>
        </g>
        <text x="460" y="310" class="tissue-label" text-anchor="middle">Red Fruit</text>

        <!-- ============ AXILLARY BUD ============ -->
        <g id="tissue-aux_bud" class="tissue-part" onclick="selectTissue('aux_bud')" data-tissue="aux_bud">
          <!-- Small runner/stolon emerging -->
          <path d="M340,390 Q365,388 385,392" fill="none" stroke="#7CB342" stroke-width="3" stroke-linecap="round"/>
          <!-- Bud -->
          <ellipse cx="395" cy="392" rx="12" ry="8" fill="#8BC34A"/>
          <path d="M390,387 Q395,378 400,387" fill="#6B9E3A"/>
          <path d="M387,389 Q393,380 398,388" fill="#7CB342"/>
        </g>
        <text x="420" y="405" class="tissue-label" font-size="9">Aux. Bud</text>
        <text x="420" y="416" class="tissue-label" font-size="7" fill="#ff6b6b">(N=1)</text>

        <!-- TPM value labels (positioned over each tissue) -->
        <text id="tpm-label-leaf" class="tissue-tpm-label" x="260" y="195" opacity="0"></text>
        <text id="tpm-label-green_fruit" class="tissue-tpm-label" x="185" y="350" opacity="0"></text>
        <text id="tpm-label-red_fruit" class="tissue-tpm-label" x="420" y="360" opacity="0"></text>
        <text id="tpm-label-crown" class="tissue-tpm-label" x="300" y="398" opacity="0"></text>
        <text id="tpm-label-roots" class="tissue-tpm-label" x="300" y="560" opacity="0"></text>
        <text id="tpm-label-aux_bud" class="tissue-tpm-label" x="395" y="395" opacity="0" font-size="8"></text>
      </svg>
    </div>
  </div>

  <!-- RIGHT PANEL: Info -->
  <div class="panel scroll-custom" id="info-panel">
    <div class="panel-title">Expression Details</div>

    <div id="no-selection" class="no-data-msg">
      <div style="font-size: 2rem; margin-bottom: 12px;">🧬</div>
      <p>Select a gene from the left panel<br>or click a tissue on the plant</p>
    </div>

    <div id="gene-info" style="display:none;">
      <div class="info-section">
        <h3 id="info-gene-name">—</h3>
        <div class="info-row">
          <span class="label">Gene ID</span>
          <span class="value" id="info-gene-id">—</span>
        </div>
        <div class="info-row">
          <span class="label">Subfamily</span>
          <span class="value" id="info-subfamily">—</span>
        </div>
        <div class="info-row">
          <span class="label">Annotation</span>
          <span class="value" id="info-source">—</span>
        </div>
        <div class="info-row">
          <span class="label">TMHs</span>
          <span class="value" id="info-tmhs">—</span>
        </div>
        <div class="info-row" id="reannotation-row" style="display:none;">
          <span class="label">Status</span>
          <span class="value" style="color:var(--accent-warm)">⚠ Needs reannotation</span>
        </div>
      </div>

      <div class="info-section">
        <h3>Expression by Tissue (TPM)</h3>
        <div class="expr-bar-container" id="expr-bars"></div>
      </div>

      <div class="info-section" id="tissue-detail" style="display:none;">
        <h3 id="tissue-detail-title">Selected Tissue</h3>
        <div class="info-row">
          <span class="label">Mean TPM</span>
          <span class="value" id="detail-tpm">—</span>
        </div>
        <div class="info-row">
          <span class="label">SD</span>
          <span class="value" id="detail-sd">—</span>
        </div>
        <div class="info-row">
          <span class="label">Replicates</span>
          <span class="value" id="detail-n">—</span>
        </div>
      </div>
    </div>
  </div>
</div>

<!-- Floating tooltip -->
<div class="svg-tooltip" id="svg-tooltip">
  <div class="tt-tissue" id="tt-tissue">—</div>
  <div class="tt-tpm" id="tt-tpm">—</div>
  <div class="tt-detail" id="tt-detail">—</div>
</div>

<script>
// ---- Embedded data ----------------------------------------------------------
const DATA = {data_json};

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
let currentMode = 'single';
let currentScale = 'absolute';
let selectedGeneId = null;
let selectedTissue = null;
let selectedSubfamily = 'ALL';

// ---- Color interpolation ----------------------------------------------------
function tpmToColor(tpm, maxTpm) {{
  if (maxTpm <= 0) maxTpm = 1;
  const ratio = Math.min(tpm / maxTpm, 1);
  // Deep blue → cyan → yellow → red palette
  const stops = [
    [0.0, [26, 26, 46]],
    [0.15, [22, 33, 62]],
    [0.3, [15, 52, 96]],
    [0.5, [78, 160, 180]],
    [0.7, [233, 69, 96]],
    [0.85, [255, 107, 107]],
    [1.0, [255, 217, 61]]
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

// ---- Initialize -------------------------------------------------------------
function init() {{
  populateSubfamilySelect();
  populateGeneList();
  document.getElementById('max-tpm-label').textContent = DATA.global_max_tpm.toFixed(0) + ' TPM';

  // Select first gene if available
  const genes = getFilteredGenes();
  if (genes.length > 0) {{
    selectGene(genes[0].id);
  }}
}}

function populateSubfamilySelect() {{
  const sel = document.getElementById('subfamily-select');
  DATA.subfamilies.forEach(sf => {{
    const opt = document.createElement('option');
    opt.value = sf;
    opt.textContent = sf + ' (' + (DATA.genes_by_subfamily[sf] || []).length + ' genes)';
    opt.style.color = SUBFAMILY_COLORS[sf] || 'inherit';
    sel.appendChild(opt);
  }});
}}

function getFilteredGenes() {{
  const search = document.getElementById('gene-search').value.toLowerCase();
  let genes = [];

  if (selectedSubfamily === 'ALL') {{
    for (const sf of DATA.subfamilies) {{
      genes = genes.concat(DATA.genes_by_subfamily[sf] || []);
    }}
  }} else {{
    genes = DATA.genes_by_subfamily[selectedSubfamily] || [];
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

    const sfColor = SUBFAMILY_COLORS[meta?.subfamily] || '#666';
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

// ---- Mode & Scale -----------------------------------------------------------
function setMode(mode) {{
  currentMode = mode;
  document.getElementById('mode-single').classList.toggle('active', mode === 'single');
  document.getElementById('mode-subfamily').classList.toggle('active', mode === 'subfamily');

  if (mode === 'subfamily' && selectedSubfamily !== 'ALL') {{
    updatePlantSubfamily(selectedSubfamily);
  }} else if (selectedGeneId) {{
    updatePlant(selectedGeneId);
  }}
}}

function setScale(scale) {{
  currentScale = scale;
  document.getElementById('scale-abs').classList.toggle('active', scale === 'absolute');
  document.getElementById('scale-rel').classList.toggle('active', scale === 'relative');

  if (selectedGeneId) updatePlant(selectedGeneId);
}}

// ---- Selection handlers -----------------------------------------------------
function onSubfamilyChange() {{
  selectedSubfamily = document.getElementById('subfamily-select').value;
  populateGeneList();

  if (currentMode === 'subfamily' && selectedSubfamily !== 'ALL') {{
    updatePlantSubfamily(selectedSubfamily);
  }}
}}

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

  // Update tissue visual selection
  document.querySelectorAll('.tissue-part').forEach(el => {{
    el.classList.toggle('selected', el.dataset.tissue === tissue);
  }});

  // Show tissue detail
  if (selectedGeneId) {{
    updateTissueDetail(selectedGeneId, tissue);
  }}
}}

// ---- Update plant visualization ---------------------------------------------
function updatePlant(geneId) {{
  const expr = DATA.expression[geneId];
  if (!expr) return;

  let maxTpm = DATA.global_max_tpm;
  if (currentScale === 'relative') {{
    maxTpm = Math.max(...Object.values(expr).map(e => e.mean_tpm), 0.1);
  }}

  document.getElementById('max-tpm-label').textContent = maxTpm.toFixed(0) + ' TPM';

  const tissues = ['leaf', 'roots', 'crown', 'green_fruit', 'red_fruit', 'aux_bud'];
  tissues.forEach(tissue => {{
    const group = document.getElementById('tissue-' + tissue);
    const label = document.getElementById('tpm-label-' + tissue);
    if (!group) return;

    const tpm = expr[tissue]?.mean_tpm || 0;
    const color = tpmToColor(tpm, maxTpm);

    // Color all paths and shapes in the tissue group
    group.querySelectorAll('path, ellipse, circle, rect').forEach(el => {{
      el.style.fill = color;
    }});
    // Keep root paths as strokes
    if (tissue === 'roots') {{
      group.querySelectorAll('path').forEach(el => {{
        el.style.fill = 'none';
        el.style.stroke = color;
      }});
    }}

    // Update TPM label
    if (label) {{
      label.textContent = tpm.toFixed(1);
      label.setAttribute('opacity', tpm > 0 ? '1' : '0.4');
    }}
  }});
}}

function updatePlantSubfamily(subfamily) {{
  const genes = DATA.genes_by_subfamily[subfamily] || [];
  if (genes.length === 0) return;

  // Calculate mean TPM per tissue across all genes in subfamily
  const tissues = ['leaf', 'roots', 'crown', 'green_fruit', 'red_fruit', 'aux_bud'];
  const meanExpr = {{}};

  tissues.forEach(tissue => {{
    let sum = 0, count = 0;
    genes.forEach(g => {{
      const expr = DATA.expression[g.id];
      if (expr && expr[tissue]) {{
        sum += expr[tissue].mean_tpm;
        count++;
      }}
    }});
    meanExpr[tissue] = count > 0 ? sum / count : 0;
  }});

  let maxTpm = currentScale === 'absolute' ?
    DATA.global_max_tpm :
    Math.max(...Object.values(meanExpr), 0.1);

  document.getElementById('max-tpm-label').textContent = maxTpm.toFixed(0) + ' TPM';

  tissues.forEach(tissue => {{
    const group = document.getElementById('tissue-' + tissue);
    const label = document.getElementById('tpm-label-' + tissue);
    if (!group) return;

    const tpm = meanExpr[tissue];
    const color = tpmToColor(tpm, maxTpm);

    group.querySelectorAll('path, ellipse, circle, rect').forEach(el => {{
      el.style.fill = color;
    }});
    if (tissue === 'roots') {{
      group.querySelectorAll('path').forEach(el => {{
        el.style.fill = 'none';
        el.style.stroke = color;
      }});
    }}

    if (label) {{
      label.textContent = tpm.toFixed(1);
      label.setAttribute('opacity', tpm > 0 ? '1' : '0.4');
    }}
  }});

  // Update info panel for subfamily mode
  document.getElementById('no-selection').style.display = 'none';
  document.getElementById('gene-info').style.display = 'block';
  document.getElementById('info-gene-name').textContent = subfamily + ' subfamily (mean)';
  document.getElementById('info-gene-id').textContent = genes.length + ' genes';
  document.getElementById('info-subfamily').textContent = subfamily;
  document.getElementById('info-source').textContent = '—';
  document.getElementById('info-tmhs').textContent = '—';
  document.getElementById('reannotation-row').style.display = 'none';

  updateExprBars(meanExpr, maxTpm);
}}

// ---- Update info panel ------------------------------------------------------
function updateInfoPanel(geneId) {{
  const meta = DATA.metadata[geneId];
  const expr = DATA.expression[geneId];
  if (!meta || !expr) return;

  document.getElementById('info-gene-name').textContent = meta.name;
  document.getElementById('info-gene-id').textContent = geneId;
  document.getElementById('info-subfamily').textContent = meta.subfamily;
  document.getElementById('info-source').textContent = meta.fuente_seq;
  document.getElementById('info-tmhs').textContent = meta.TMHs;
  document.getElementById('reannotation-row').style.display = meta.needs_reannotation ? '' : 'none';

  let maxTpm = currentScale === 'absolute' ?
    DATA.global_max_tpm :
    Math.max(...Object.values(expr).map(e => e.mean_tpm), 0.1);

  const exprMap = {{}};
  Object.keys(expr).forEach(t => {{
    exprMap[t] = expr[t].mean_tpm;
  }});

  updateExprBars(exprMap, maxTpm);
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

  document.querySelectorAll('.tissue-part').forEach(el => {{
    el.addEventListener('mouseenter', (e) => {{
      const tissue = el.dataset.tissue;
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

    el.addEventListener('mousemove', (e) => {{
      tooltip.style.left = (e.clientX + 16) + 'px';
      tooltip.style.top = (e.clientY - 10) + 'px';
    }});

    el.addEventListener('mouseleave', () => {{
      tooltip.classList.remove('visible');
    }});
  }});
}}

// ---- Export as PNG -----------------------------------------------------------
function exportImage() {{
  const svg = document.getElementById('plant-svg');
  const svgData = new XMLSerializer().serializeToString(svg);
  const canvas = document.createElement('canvas');
  const ctx = canvas.getContext('2d');
  const img = new Image();

  canvas.width = 1200;
  canvas.height = 1400;

  img.onload = function() {{
    // Dark background
    ctx.fillStyle = '#0f1117';
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    ctx.drawImage(img, 0, 0, 1200, 1400);

    // Add title
    ctx.fillStyle = '#e8eaed';
    ctx.font = 'bold 20px Inter, sans-serif';
    const geneLabel = selectedGeneId ?
      (DATA.metadata[selectedGeneId]?.name || selectedGeneId) : 'No gene selected';
    ctx.fillText('Fragaria × ananassa — ' + geneLabel, 20, 30);

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
print("# === 10_generate_efp_viewer.py complete ===")
