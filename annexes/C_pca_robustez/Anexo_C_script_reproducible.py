import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.lines import Line2D
import numpy as np
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.common import config

OUT_DIR = Path(__file__).resolve().parent
PCA_COORDS = config.RESULTS_DIR / "profiling_aqp_motifs_final" / "PCA_Coordenadas_Finales.csv"
TABLA = config.CURADO_DIR / "tabla_aquaporinas_traduccion.tabular"

pca = {}
with open(PCA_COORDS, encoding="utf-8") as f:
    for row in csv.DictReader(f):
        pca[row["ID"]] = {"PC1":float(row["PC1"]),"PC2":float(row["PC2"]),"subfam":row["Subfamilia"]}

tabla = {}
with open(TABLA, encoding="utf-8") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        tabla[row["gene_id"]] = row

CLEAN = {"IGUALES","EXONERATE","PEPTIDE_EMP","MANUAL_CURATED"}
data = []
for gid, p in pca.items():
    t = tabla.get(gid, {})
    sf = t.get("subfamilia_phylo", p["subfam"])
    v = t.get("veredicto", "?")
    cat = "CLEAN" if (v in CLEAN and t.get("TMHs")=="6" and gid != "Fxa5Ag03930") else "PARTIAL"
    data.append({"id":gid,"pc1":p["PC1"],"pc2":p["PC2"],"subfam":sf,"cat":cat})

colors = {"PIP":"#E74C3C","TIP":"#3498DB","NIP":"#2ECC71","SIP":"#F39C12","XIP":"#9B59B6"}

def add_ellipse(ax, xs, ys, color, alpha=0.12, lw=2.5, ls="-", fill=True):
    if len(xs) < 3: return None, None, None
    cov = np.cov(xs, ys)
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = eigvals.argsort()[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]
    angle = np.degrees(np.arctan2(eigvecs[1,0], eigvecs[0,0]))
    w, h = 2 * np.sqrt(eigvals * 5.991)
    cx, cy = np.mean(xs), np.mean(ys)
    e = Ellipse((cx, cy), w, h, angle=angle,
                facecolor=color if fill else "none",
                edgecolor=color, alpha=alpha if fill else 1.0,
                lw=lw, linestyle=ls)
    ax.add_patch(e)
    return cx, cy, e

# Plot 1: dos paneles lado a lado
fig, axes = plt.subplots(1, 2, figsize=(20, 9), sharex=True, sharey=True)
for ax, use_all in zip(axes, [False, True]):
    for sf in ["NIP","PIP","SIP","TIP","XIP"]:
        cleans = [(d["pc1"],d["pc2"]) for d in data if d["subfam"]==sf and d["cat"]=="CLEAN"]
        parts  = [(d["pc1"],d["pc2"]) for d in data if d["subfam"]==sf and d["cat"]=="PARTIAL"]
        if cleans:
            xs, ys = zip(*cleans)
            ax.scatter(xs, ys, s=40, c=colors[sf], alpha=0.85, edgecolor="white", linewidth=0.5, zorder=3)
        if parts:
            xs, ys = zip(*parts)
            ax.scatter(xs, ys, s=80, c=colors[sf], alpha=0.5, edgecolor="black", linewidth=1.2, marker="X", zorder=4)
        pts = (cleans + parts) if use_all else cleans
        if len(pts) >= 3:
            xs, ys = zip(*pts)
            cx, cy, _ = add_ellipse(ax, xs, ys, colors[sf], alpha=0.15, lw=2.5)
            ax.text(cx, cy, sf, fontsize=14, fontweight="bold", color=colors[sf],
                    ha="center", va="center",
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor=colors[sf], alpha=0.85))
    for pid in ["Fxa6Bg00715","Fxa6Cg01391","Fxa5Bg03706","Fxa3Ag00841","Fxa5Ag03930"]:
        d = next((x for x in data if x["id"]==pid), None)
        if d:
            ax.annotate(pid, (d["pc1"], d["pc2"]), fontsize=8, fontweight="bold",
                       xytext=(8,8), textcoords="offset points",
                       bbox=dict(boxstyle="round,pad=0.2", facecolor="yellow", alpha=0.9))
    title = ("VERSION A: elipses solo con CLEAN (n=121)\nLas parciales NO afectan al calculo"
             if not use_all else
             "VERSION B: elipses con TODAS las candidatas (n=144)\nLas parciales SI participan en el calculo")
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xlabel("PC1", fontsize=11)
    ax.set_ylabel("PC2", fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color="gray", linewidth=0.5)
    ax.axvline(0, color="gray", linewidth=0.5)

legend_items = [
    Line2D([0],[0], marker="o", color="w", markerfacecolor="gray", markersize=10, label="CLEAN (121 funcionales)"),
    Line2D([0],[0], marker="X", color="w", markerfacecolor="gray", markeredgecolor="black", markersize=12, label="PARTIAL (23 descartadas)"),
]
fig.legend(handles=legend_items, loc="lower center", ncol=2, fontsize=11, bbox_to_anchor=(0.5, -0.02))
plt.suptitle("Comparacion PCA: elipses al 95% segun se entrenen con CLEAN o con TODAS",
             fontsize=14, fontweight="bold", y=1.00)
plt.tight_layout()
plt.savefig(OUT_DIR / "Anexo_C_figura_paneles_AvsB.png", dpi=180, bbox_inches="tight")
plt.close()
print("Guardado:", OUT_DIR / "Anexo_C_figura_paneles_AvsB.png")

# Plot 2: superpuesto
fig, ax = plt.subplots(figsize=(13, 9))
for sf in ["NIP","PIP","SIP","TIP","XIP"]:
    cleans = [(d["pc1"],d["pc2"]) for d in data if d["subfam"]==sf and d["cat"]=="CLEAN"]
    parts  = [(d["pc1"],d["pc2"]) for d in data if d["subfam"]==sf and d["cat"]=="PARTIAL"]
    if cleans:
        xs, ys = zip(*cleans)
        ax.scatter(xs, ys, s=35, c=colors[sf], alpha=0.85, edgecolor="white", linewidth=0.4, zorder=3)
        add_ellipse(ax, xs, ys, colors[sf], alpha=1.0, lw=2.5, ls="-", fill=False)
    if parts:
        xs, ys = zip(*parts)
        ax.scatter(xs, ys, s=70, c=colors[sf], alpha=0.5, edgecolor="black", linewidth=1.0, marker="X", zorder=4)
    pts_all = cleans + parts
    if len(pts_all) >= 3:
        xs, ys = zip(*pts_all)
        add_ellipse(ax, xs, ys, colors[sf], alpha=1.0, lw=2.0, ls="--", fill=False)

for pid in ["Fxa6Bg00715","Fxa6Cg01391","Fxa5Bg03706","Fxa3Ag00841","Fxa5Ag03930"]:
    d = next((x for x in data if x["id"]==pid), None)
    if d:
        ax.annotate(pid, (d["pc1"], d["pc2"]), fontsize=8, fontweight="bold",
                   xytext=(8,8), textcoords="offset points",
                   bbox=dict(boxstyle="round,pad=0.2", facecolor="yellow", alpha=0.9))

ax.set_title("PCA - elipses al 95%: linea solida = solo CLEAN (n=121); linea discontinua = todas (n=144)",
             fontsize=12, fontweight="bold")
ax.set_xlabel("PC1", fontsize=11)
ax.set_ylabel("PC2", fontsize=11)
ax.grid(True, alpha=0.3)
ax.axhline(0, color="gray", linewidth=0.5)
ax.axvline(0, color="gray", linewidth=0.5)
plt.tight_layout()
plt.savefig(OUT_DIR / "Anexo_C_figura_elipses_AvsB.png", dpi=180, bbox_inches="tight")
plt.close()
print("Guardado:", OUT_DIR / "Anexo_C_figura_elipses_AvsB.png")
