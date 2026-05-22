import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from scipy import stats
from pathlib import Path

base = Path(os.environ.get("TFG_RNA_SEQ_ROOT","/home/noe/work/RNA-seq_test")+"/results/basal_aquaporins")
df = pd.read_csv(base / "basal_aquaporins_summary.csv")
if "needs_reannotation" in df.columns:
    df = df[df["needs_reannotation"].astype(str).str.upper() == "FALSE"].copy()
df = df.rename(columns={"subfamilia_phylo": "subfamily", "mean_tpm": "tpm"})

tissue_order = ["red_fruit", "roots", "green_fruit", "crown", "aux_bud", "leaf"]
tissue_labels = ["Fruto rojo", "Raiz", "Fruto verde", "Corona", "Yema axilar*", "Hoja"]
subfams = ["PIP", "TIP", "NIP", "SIP", "XIP"]
sf_colors = {"PIP": "#E74C3C", "TIP": "#3498DB", "NIP": "#2ECC71", "SIP": "#F39C12", "XIP": "#9B59B6"}

pivot = df.pivot_table(index=["gene_id", "subfamily"], columns="tissue", values="tpm", aggfunc="mean")
pivot["TOTAL"] = pivot[tissue_order].sum(axis=1)
totals_df = pivot.reset_index()[["gene_id", "subfamily", "TOTAL"]].rename(columns={"TOTAL": "tpm"})
totals_df["tissue"] = "TOTAL"

n_sub = len(subfams)
sub_spacing = 1.0 / n_sub
group_gap = 1.05  # menos espacio entre tejidos
fig_w = 22  # mucho mas ancho
fig, ax = plt.subplots(figsize=(fig_w, 8.5), dpi=140)
ax.set_yscale("symlog", linthresh=1)
ax.set_ylim(-0.5, 12000)
rng = np.random.default_rng(42)
xtick_pos = []

def draw_cell(x_center, vals, color):
    # Background subtle
    box_w_total = sub_spacing * 0.92
    ax.add_patch(Rectangle((x_center - box_w_total/2, -0.5), box_w_total, 12000.5,
                           facecolor=color, alpha=0.07, edgecolor=color, linewidth=0.4, zorder=0))
    if len(vals) == 0:
        return
    if len(vals) == 1:
        ax.scatter([x_center], vals, s=40, color=color, edgecolor="black", linewidth=0.5, alpha=0.95, zorder=4)
        mean_v = float(vals[0])
        ax.hlines(mean_v, x_center - sub_spacing * 0.4, x_center + sub_spacing * 0.4,
                  color="black", linewidth=2.6, zorder=5)
        ax.text(x_center, mean_v * 1.45 if mean_v >= 1 else mean_v + 0.5,
                f"{mean_v:.1f}" if mean_v >= 1 else f"{mean_v:.2f}",
                ha="center", va="bottom", fontsize=9, color=color, fontweight="bold")
        return
    # Densidad kernel bilateral
    if vals.std() > 1e-6 and len(vals) >= 3:
        # transformacion para KDE en symlog
        with np.errstate(divide="ignore"):
            log_vals = np.where(vals >= 1, np.log10(vals) + 1, vals)
        if log_vals.std() > 1e-6:
            kde = stats.gaussian_kde(log_vals, bw_method=0.32)
            ys_log = np.linspace(max(log_vals.min() - 0.4, -0.5), log_vals.max() + 0.4, 80)
            density = kde(ys_log)
            density = density / density.max() * (sub_spacing * 0.42)
            ys = np.where(ys_log >= 1, 10 ** (ys_log - 1), ys_log)
            # Violin bilateral (ancho completo de la sub-columna)
            ax.fill_betweenx(ys, x_center - density, x_center + density,
                             facecolor=color, alpha=0.30, edgecolor=color, linewidth=0.7, zorder=1)
    # Scatter con jitter ancho ocupando casi toda la sub-columna
    jitter = rng.uniform(-sub_spacing * 0.38, sub_spacing * 0.38, size=len(vals))
    ax.scatter(np.full(len(vals), x_center) + jitter, vals,
               s=20, color=color, edgecolor="black", linewidth=0.3, alpha=0.92, zorder=3)
    # Linea media cruzando todo el ancho de la sub-columna
    mean_v = float(np.mean(vals))
    ax.hlines(mean_v, x_center - sub_spacing * 0.42, x_center + sub_spacing * 0.42,
              color="black", linewidth=3.2, zorder=5)
    ax.hlines(mean_v, x_center - sub_spacing * 0.42, x_center + sub_spacing * 0.42,
              color=color, linewidth=1.6, zorder=6)
    ax.text(x_center, mean_v * 1.45 if mean_v >= 1 else mean_v + 0.5,
            f"{mean_v:.1f}" if mean_v >= 1 else f"{mean_v:.2f}",
            ha="center", va="bottom", fontsize=9, color=color, fontweight="bold")

def draw_group(group_center, vals_lookup):
    for s_idx, sf in enumerate(subfams):
        x_center = group_center + (s_idx - (n_sub - 1)/2) * sub_spacing
        draw_cell(x_center, vals_lookup(sf), sf_colors[sf])

for t_idx, tissue in enumerate(tissue_order):
    group_center = t_idx * group_gap
    xtick_pos.append(group_center)
    draw_group(group_center, lambda sf, t=tissue: df[(df["subfamily"] == sf) & (df["tissue"] == t)]["tpm"].dropna().values)

total_sep_x = (len(tissue_order) - 0.5) * group_gap + 0.25
ax.axvline(total_sep_x, color="dimgray", linewidth=2.2, alpha=0.85, zorder=0.7)
total_center = len(tissue_order) * group_gap + 0.45
xtick_pos.append(total_center)
draw_group(total_center, lambda sf: totals_df[totals_df["subfamily"] == sf]["tpm"].dropna().values)

for t_idx in range(len(tissue_order) - 1):
    sep_x = (t_idx + 0.5) * group_gap
    ax.axvline(sep_x, color="lightgray", linewidth=0.8, alpha=0.4, zorder=0.5)

ax.set_xticks(xtick_pos)
ax.set_xticklabels(tissue_labels + ["TOTAL\n(suma 6 tejidos)"], fontsize=12)
ax.set_ylabel("TPM medio por gen (escala symlog)", fontsize=12)
ax.axhline(1, color="gray", linestyle="--", linewidth=0.8, alpha=0.7, zorder=1)
ax.grid(axis="y", linestyle=":", alpha=0.4, zorder=0)
ax.set_axisbelow(True)

legend_elems = [Line2D([0], [0], marker="o", color="w", label=sf,
                       markerfacecolor=sf_colors[sf], markeredgecolor="black", markersize=11)
                for sf in subfams]
legend_elems += [
    Line2D([0], [0], color="black", linewidth=3.2, label="Media por gen"),
    Line2D([0], [0], marker="s", color="w", markerfacecolor="gray", markeredgecolor="gray", alpha=0.4, markersize=12, label="Densidad kernel (violin)"),
]
ax.legend(handles=legend_elems, loc="upper right", ncol=2, fontsize=10.5, framealpha=0.95)

fig.tight_layout()
out_pdf = base / "figura6_perfiles_subfamilia.pdf"
out_png = base / "figura6_perfiles_subfamilia.png"
fig.savefig(str(out_pdf), dpi=200, bbox_inches="tight")
fig.savefig(str(out_png), dpi=200, bbox_inches="tight")
print("PDF:", out_pdf.stat().st_size)
print("PNG:", out_png.stat().st_size)
