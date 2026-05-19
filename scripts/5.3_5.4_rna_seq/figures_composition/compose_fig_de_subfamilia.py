import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from scipy import stats
from pathlib import Path

base = Path("/home/noe/work/RNA-seq_test/results")
sf_colors = {"PIP": "#E74C3C", "TIP": "#3498DB", "NIP": "#2ECC71", "SIP": "#F39C12", "XIP": "#9B59B6"}
subfams = ["PIP", "TIP", "NIP", "SIP", "XIP"]

leaf_df = pd.read_csv(base / "de_leaf" / "de_aquaporins_leaf.csv")
roots_df = pd.read_csv(base / "de_roots" / "de_aquaporins_roots.csv")
for d in (leaf_df, roots_df):
    d["log2FoldChange"] = pd.to_numeric(d["log2FoldChange"], errors="coerce")
    d["PAdj"] = pd.to_numeric(d["PAdj"], errors="coerce")
leaf_df = leaf_df.dropna(subset=["log2FoldChange"])
roots_df = roots_df.dropna(subset=["log2FoldChange"])

tissues = [("Hoja (3 ctrl vs 3 PEG)", leaf_df), ("Raiz (2 ctrl vs 3 PEG)", roots_df)]

fig, axes = plt.subplots(1, 2, figsize=(18, 9), dpi=140, sharey=True)
rng = np.random.default_rng(7)

YMIN, YMAX = -7.5, 7
SUB_W = 0.85  # ancho efectivo de cada subfamilia (centro a borde)

for ax_idx, (ax, (label, d)) in enumerate(zip(axes, tissues)):
    # Sombreado direccional
    ax.axhspan(0, YMAX, facecolor="#88c98a", alpha=0.07, zorder=0)
    ax.axhspan(YMIN, 0, facecolor="#e07b7b", alpha=0.07, zorder=0)

    # Pirate plot per subfamily
    for s_idx, sf in enumerate(subfams):
        # Sombra de subfamilia
        rect = Rectangle((s_idx - SUB_W/2, YMIN), SUB_W, YMAX - YMIN,
                         facecolor=sf_colors[sf], alpha=0.06, edgecolor=sf_colors[sf],
                         linewidth=0.5, zorder=0.5)
        ax.add_patch(rect)

        sub = d[d["subfamilia_phylo"] == sf]
        vals = sub["log2FoldChange"].dropna().values
        if len(vals) == 0:
            continue
        # Violin bilateral
        if len(vals) >= 3 and vals.std() > 1e-6:
            kde = stats.gaussian_kde(vals, bw_method=0.4)
            ys = np.linspace(vals.min() - 0.5, vals.max() + 0.5, 80)
            density = kde(ys)
            density = density / density.max() * (SUB_W * 0.40)
            ax.fill_betweenx(ys, s_idx - density, s_idx + density,
                             facecolor=sf_colors[sf], alpha=0.25, edgecolor=sf_colors[sf], linewidth=0.7, zorder=1)
        # Significant mask
        sig_mask = (sub["PAdj"] < 0.05) & (sub["log2FoldChange"].abs() > 1)
        sig_mask = sig_mask.fillna(False).values
        # Scatter
        jitter = rng.uniform(-SUB_W * 0.35, SUB_W * 0.35, size=len(vals))
        ns = ~sig_mask
        ax.scatter(np.full(ns.sum(), s_idx) + jitter[ns], vals[ns],
                   s=22, facecolor="white", edgecolor=sf_colors[sf], linewidth=0.8,
                   alpha=0.75, zorder=2)
        ax.scatter(np.full(sig_mask.sum(), s_idx) + jitter[sig_mask], vals[sig_mask],
                   s=50, color=sf_colors[sf], edgecolor="black", linewidth=0.9,
                   alpha=0.95, zorder=4)
        # Etiquetas en puntos significativos: FaXXXn-FxaXXXXX
        sig_indices = np.where(sig_mask)[0]
        sub_reset = sub.reset_index(drop=True)
        for i_sig in sig_indices:
            x_pt = s_idx + jitter[i_sig]
            y_pt = vals[i_sig]
            gene_id = str(sub_reset.iloc[i_sig]["name"])
            aqp_subfam = str(sub_reset.iloc[i_sig]["aqp_family_subfamily"])
            tag_label = f"{aqp_subfam}-{gene_id}"
            tag_side = "left" if jitter[i_sig] >= 0 else "right"
            tag_dx = 7 if tag_side == "left" else -7
            ax.annotate(tag_label, xy=(x_pt, y_pt), xytext=(tag_dx, 0),
                        textcoords="offset points",
                        fontsize=6, ha=tag_side, va="center",
                        color="#222222", alpha=0.92, zorder=7,
                        bbox=dict(facecolor="white", edgecolor="none",
                                  alpha=0.7, pad=0.6))
        # Mean
        mean_v = float(np.mean(vals))
        ax.hlines(mean_v, s_idx - SUB_W * 0.40, s_idx + SUB_W * 0.40, color="black", linewidth=3.4, zorder=5)
        ax.hlines(mean_v, s_idx - SUB_W * 0.40, s_idx + SUB_W * 0.40, color=sf_colors[sf], linewidth=1.8, zorder=6)
        # Sig count
        n_sig = int(sig_mask.sum())
        if n_sig > 0:
            ax.text(s_idx, YMAX - 0.6, f"{n_sig} sig", ha="center", va="center",
                    fontsize=11, fontweight="bold", color=sf_colors[sf])

    ax.axhline(0, color="black", linewidth=1.0, alpha=0.7, zorder=2)
    ax.axhline(1, color="gray", linestyle="--", linewidth=0.6, alpha=0.6, zorder=1)
    ax.axhline(-1, color="gray", linestyle="--", linewidth=0.6, alpha=0.6, zorder=1)
    ax.set_xticks(range(len(subfams)))
    ax.set_xticklabels(subfams, fontsize=13, fontweight="bold")
    ax.set_title(label, fontsize=13.5, fontweight="bold")
    ax.grid(axis="y", linestyle=":", alpha=0.4, zorder=0)
    ax.set_axisbelow(True)
    ax.set_ylim(YMIN, YMAX)
    ax.set_xlim(-0.6, len(subfams) - 0.4)

    # Direction labels on the right side (only on rightmost subplot)
    if ax_idx == 1:
        ax.text(1.03, 0.75, "INDUCIDO\npor estres\n(estres > control)",
                transform=ax.transAxes, fontsize=11, fontweight="bold", color="#2d8a30",
                ha="left", va="center",
                bbox=dict(facecolor="#e6f4e6", edgecolor="#88c98a", boxstyle="round,pad=0.4"))
        ax.text(1.03, 0.25, "REPRIMIDO\npor estres\n(estres < control)",
                transform=ax.transAxes, fontsize=11, fontweight="bold", color="#a82020",
                ha="left", va="center",
                bbox=dict(facecolor="#fbeaea", edgecolor="#e07b7b", boxstyle="round,pad=0.4"))
        ax.annotate("", xy=(1.18, 0.92), xytext=(1.18, 0.58),
                    xycoords="axes fraction", arrowprops=dict(arrowstyle="->", color="#2d8a30", lw=2))
        ax.annotate("", xy=(1.18, 0.08), xytext=(1.18, 0.42),
                    xycoords="axes fraction", arrowprops=dict(arrowstyle="->", color="#a82020", lw=2))

axes[0].set_ylabel("log2(FoldChange) — Estres / Control", fontsize=12.5)

legend_elems = [
    Line2D([0], [0], marker="o", color="w", markerfacecolor="lightgray", markeredgecolor="gray", markersize=10, label="No significativo"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="gray", markeredgecolor="black", markersize=12, label="padj<0,05 & |LFC|>1"),
    Line2D([0], [0], color="black", linewidth=3.2, label="Media por subfamilia"),
    Line2D([0], [0], color="black", linewidth=1.0, label="LFC = 0 (sin cambio)"),
    Line2D([0], [0], color="gray", linestyle="--", linewidth=0.6, label="|LFC| = 1 (umbral)"),
]
axes[0].legend(handles=legend_elems, loc="lower left", fontsize=9.5, framealpha=0.95)

fig.tight_layout(rect=[0, 0, 0.93, 1.0])

out_pdf = base / "figura_de_stripplot.pdf"
out_png = base / "figura_de_stripplot.png"
fig.savefig(str(out_pdf), dpi=200, bbox_inches="tight")
fig.savefig(str(out_png), dpi=200, bbox_inches="tight")
print("PDF:", out_pdf.stat().st_size)
print("PNG:", out_png.stat().st_size)
