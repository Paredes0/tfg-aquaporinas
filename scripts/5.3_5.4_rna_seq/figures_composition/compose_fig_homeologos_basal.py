import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from scipy import stats
from pathlib import Path

base = Path("/home/noe/work/RNA-seq_test")
hom = base / "results" / "homeolog_analysis"

hg_sum = pd.read_csv(base / "homeolog_groups_summary.tsv", sep="\t")
ct = pd.read_csv(hom / "collapsed_tpm.csv").rename(columns={"Unnamed: 0": "homeolog_group"})
do = pd.read_csv(hom / "dominance_overall.csv")

design = pd.read_csv(base / "design" / "design_basal.csv")
sample_to_tissue = dict(zip(design["sample"], design["tissue"]))

tissue_order = ["red_fruit", "roots", "green_fruit", "crown", "aux_bud", "leaf"]
tissue_labels = ["Fruto rojo", "Raíz", "Fruto verde", "Corona", "Yema axilar*", "Hoja"]
ct_long = ct.melt(id_vars="homeolog_group", var_name="sample", value_name="tpm")
ct_long["tissue"] = ct_long["sample"].map(sample_to_tissue)
ct_long = ct_long.dropna(subset=["tissue"])
ct_mean = ct_long.groupby(["homeolog_group", "tissue"], as_index=False)["tpm"].mean()
ct_mean = ct_mean.merge(hg_sum[["homeolog_group", "family"]], on="homeolog_group")

pivot = ct_mean.pivot_table(index=["homeolog_group", "family"], columns="tissue", values="tpm")
pivot["TOTAL"] = pivot[tissue_order].sum(axis=1)
totals_df = pivot.reset_index()[["homeolog_group", "family", "TOTAL"]].rename(columns={"TOTAL": "tpm", "family": "subfamily"})
totals_df["tissue"] = "TOTAL"

ct_mean = ct_mean.rename(columns={"family": "subfamily"})

subfams = ["PIP", "TIP", "NIP", "SIP", "XIP"]
sf_colors = {"PIP": "#E74C3C", "TIP": "#3498DB", "NIP": "#2ECC71", "SIP": "#F39C12", "XIP": "#9B59B6"}
# Paleta NUEVA para subgenomas que NO choca con las subfamilias
sg_colors = {"A": "#6D4C41", "B": "#C2185B", "C": "#00838F", "D": "#455A64"}

fig = plt.figure(figsize=(22, 15), dpi=140)
gs = gridspec.GridSpec(2, 1, height_ratios=[1.0, 0.55], hspace=0.35)
ax_a = fig.add_subplot(gs[0])
ax_b = fig.add_subplot(gs[1])

n_sub = len(subfams)
sub_spacing = 1.0 / n_sub
group_gap = 1.05
ax_a.set_yscale("symlog", linthresh=1)
ax_a.set_ylim(-0.5, 12000)

rng = np.random.default_rng(42)
xtick_pos_a = []

def draw_cell(ax, x_center, vals, color):
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
        return
    if len(vals) >= 3 and vals.std() > 1e-6:
        with np.errstate(divide="ignore"):
            log_vals = np.where(vals >= 1, np.log10(vals) + 1, vals)
        if log_vals.std() > 1e-6:
            kde = stats.gaussian_kde(log_vals, bw_method=0.32)
            ys_log = np.linspace(max(log_vals.min() - 0.4, -0.5), log_vals.max() + 0.4, 80)
            density = kde(ys_log)
            density = density / density.max() * (sub_spacing * 0.42)
            ys = np.where(ys_log >= 1, 10 ** (ys_log - 1), ys_log)
            ax.fill_betweenx(ys, x_center - density, x_center + density,
                             facecolor=color, alpha=0.30, edgecolor=color, linewidth=0.7, zorder=1)
    jitter = rng.uniform(-sub_spacing * 0.38, sub_spacing * 0.38, size=len(vals))
    ax.scatter(np.full(len(vals), x_center) + jitter, vals,
               s=22, color=color, edgecolor="black", linewidth=0.3, alpha=0.92, zorder=3)
    mean_v = float(np.mean(vals))
    ax.hlines(mean_v, x_center - sub_spacing * 0.42, x_center + sub_spacing * 0.42,
              color="black", linewidth=3.2, zorder=5)
    ax.hlines(mean_v, x_center - sub_spacing * 0.42, x_center + sub_spacing * 0.42,
              color=color, linewidth=1.6, zorder=6)
    ax.text(x_center, mean_v * 1.45 if mean_v >= 1 else mean_v + 0.5,
            f"{mean_v:.1f}" if mean_v >= 1 else f"{mean_v:.2f}",
            ha="center", va="bottom", fontsize=9, color=color, fontweight="bold")

def draw_group(ax, group_center, vals_lookup):
    for s_idx, sf in enumerate(subfams):
        x_center = group_center + (s_idx - (n_sub - 1)/2) * sub_spacing
        draw_cell(ax, x_center, vals_lookup(sf), sf_colors[sf])

for t_idx, tissue in enumerate(tissue_order):
    group_center = t_idx * group_gap
    xtick_pos_a.append(group_center)
    draw_group(ax_a, group_center, lambda sf, t=tissue: ct_mean[(ct_mean["subfamily"] == sf) & (ct_mean["tissue"] == t)]["tpm"].dropna().values)

total_sep_x = (len(tissue_order) - 0.5) * group_gap + 0.25
ax_a.axvline(total_sep_x, color="dimgray", linewidth=2.2, alpha=0.85, zorder=0.7)
total_center = len(tissue_order) * group_gap + 0.45
xtick_pos_a.append(total_center)
draw_group(ax_a, total_center, lambda sf: totals_df[totals_df["subfamily"] == sf]["tpm"].dropna().values)

for t_idx in range(len(tissue_order) - 1):
    ax_a.axvline((t_idx + 0.5) * group_gap, color="lightgray", linewidth=0.8, alpha=0.4, zorder=0.5)

ax_a.set_xticks(xtick_pos_a)
ax_a.set_xticklabels(tissue_labels + ["TOTAL\n(suma 6 tejidos)"], fontsize=12)
ax_a.set_ylabel("TPM colapsado por grupo HG\n(SUMA de homeologos, escala symlog)", fontsize=12)
ax_a.axhline(1, color="gray", linestyle="--", linewidth=0.8, alpha=0.7, zorder=1)
ax_a.grid(axis="y", linestyle=":", alpha=0.4, zorder=0)
ax_a.set_axisbelow(True)
legend_a = [Line2D([0], [0], marker="o", color="w", label=sf,
                   markerfacecolor=sf_colors[sf], markeredgecolor="black", markersize=11)
            for sf in subfams]
legend_a += [Line2D([0], [0], color="black", linewidth=3.2, label="Media por subfamilia"),
             Line2D([0], [0], marker="s", color="w", markerfacecolor="gray", markeredgecolor="gray", alpha=0.4, markersize=12, label="Densidad kernel")]
# Leyenda FUERA del area del plot (a la derecha del axes)
ax_a.legend(handles=legend_a, loc="upper left", bbox_to_anchor=(1.005, 1.0), fontsize=10.5, framealpha=0.95, borderaxespad=0)

subgenomes = ["A", "B", "C", "D"]
ax_b.set_ylim(-0.05, 1.18)
ax_b.set_xlim(-0.5, len(subgenomes) - 0.5)
SUB_B_W = 0.85
ax_b.axhline(0.25, color="gray", linestyle="--", linewidth=0.7, alpha=0.7)
ax_b.axhline(0.5, color="dimgray", linestyle=":", linewidth=0.7, alpha=0.7)

for s_idx, sg in enumerate(subgenomes):
    rect = Rectangle((s_idx - SUB_B_W/2, -0.05), SUB_B_W, 1.23,
                     facecolor=sg_colors[sg], alpha=0.07, edgecolor=sg_colors[sg], linewidth=0.5, zorder=0)
    ax_b.add_patch(rect)
    vals = do[do["subgenome"] == sg]["mean_proportion"].dropna().values
    color = sg_colors[sg]
    if len(vals) >= 3 and vals.std() > 1e-6:
        kde = stats.gaussian_kde(vals, bw_method=0.32)
        ys = np.linspace(0, 1, 100)
        density = kde(ys)
        density = density / density.max() * (SUB_B_W * 0.42)
        ax_b.fill_betweenx(ys, s_idx - density, s_idx + density,
                           facecolor=color, alpha=0.30, edgecolor=color, linewidth=0.7, zorder=1)
    jitter = rng.uniform(-SUB_B_W * 0.38, SUB_B_W * 0.38, size=len(vals))
    ax_b.scatter(np.full(len(vals), s_idx) + jitter, vals,
                 s=28, color=color, edgecolor="black", linewidth=0.4, alpha=0.92, zorder=3)
    mean_v = float(np.mean(vals))
    ax_b.hlines(mean_v, s_idx - SUB_B_W * 0.42, s_idx + SUB_B_W * 0.42,
                color="black", linewidth=3.4, zorder=5)
    ax_b.hlines(mean_v, s_idx - SUB_B_W * 0.42, s_idx + SUB_B_W * 0.42,
                color=color, linewidth=1.8, zorder=6)
    n_dom = (do.loc[do.groupby("homeolog_group")["mean_proportion"].idxmax()]["subgenome"] == sg).sum()
    ax_b.text(s_idx, 1.10, f"Dominante en {n_dom} grupos", ha="center", va="bottom",
              fontsize=11, fontweight="bold", color=color)

ax_b.set_xticks(range(len(subgenomes)))
ax_b.set_xticklabels(["Subgenoma A", "Subgenoma B", "Subgenoma C", "Subgenoma D"], fontsize=12, fontweight="bold")
ax_b.set_ylabel("Proporción media de TPM\npor grupo homeólogo", fontsize=12)
ax_b.grid(axis="y", linestyle=":", alpha=0.4, zorder=0)
ax_b.set_axisbelow(True)

legend_b = [Line2D([0], [0], marker="o", color="w", label=f"Subgenoma {sg}",
                   markerfacecolor=sg_colors[sg], markeredgecolor="black", markersize=11)
            for sg in subgenomes]
legend_b += [Line2D([0], [0], color="black", linewidth=3.2, label="Media por subgenoma"),
             Line2D([0], [0], color="gray", linestyle="--", linewidth=0.7, label="Equireparto (0,25)")]
# Leyenda FUERA del area del plot
ax_b.legend(handles=legend_b, loc="upper left", bbox_to_anchor=(1.005, 1.0), fontsize=10.5, framealpha=0.95, borderaxespad=0)

# Ajustar para que la leyenda fuera quepa
fig.subplots_adjust(left=0.05, right=0.86, top=0.95, bottom=0.06)

out_pdf = base / "results" / "figura_homeologos_basal.pdf"
out_png = base / "results" / "figura_homeologos_basal.png"
fig.savefig(str(out_pdf), dpi=200, bbox_inches="tight")
fig.savefig(str(out_png), dpi=200, bbox_inches="tight")
print("PDF:", out_pdf.stat().st_size)
print("PNG:", out_png.stat().st_size)
