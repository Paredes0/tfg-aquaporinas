import os
"""
Figura 9 — Duplicaciones tandem pre-poliploidizacion en el clado NIP1.

Panel A: coordenadas genomicas de los pares tandem HG-7/HG-8 (chr_3, NIP1;1)
y HG-18/HG-19 (chr_6, NIP1;2) en los 4 subgenomas A/B/C/D.
Panel B: sub-arbol filogenetico del clado NIP1 extraido del arbol final
del 6.2 (Q.PLANT+R6, log L = -45.149,26), con apoyos SH-aLRT/aBayes/UFboot.

Datos verbatim del subagente A v3 (specs/2026-05-19-rnaseq-6.3.3-tandems-analysis.md).
"""
import io
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from matplotlib.patches import FancyArrow, Rectangle, FancyArrowPatch
from matplotlib.lines import Line2D
from Bio import Phylo

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
from scripts.common import config


# ---------------------------------------------------------------------------
# Paleta canonica del TFG (compose_fig_homeologos_basal.py)
# ---------------------------------------------------------------------------
SG_COLORS = {"A": "#6D4C41", "B": "#C2185B", "C": "#00838F", "D": "#455A64"}
NIP_GREEN = "#2ECC71"  # subfamilia NIP (no usado para genes, solo referencia)

# Paleta para los 4 HG — tonos saturados que no chocan con subfamilias/subgenomas
HG_COLORS = {
    "HG-7":  "#F4A261",  # naranja pastel
    "HG-8":  "#2A9D8F",  # verde-azul
    "HG-18": "#E76F51",  # coral
    "HG-19": "#264653",  # azul oscuro
}

# ---------------------------------------------------------------------------
# Datos genomicos (de tandems-analysis.md secciones A y B)
# ---------------------------------------------------------------------------
GENES_CHR3 = [
    # (gen_id, subgenoma, hg, start, end, strand)
    ("Fxa3Ag00839", "A", "HG-7", 4_074_170, 4_076_451, "-"),
    ("Fxa3Ag00840", "A", "HG-8", 4_080_131, 4_082_053, "-"),
    ("Fxa3Bg00729", "B", "HG-7", 3_630_103, 3_632_372, "-"),
    ("Fxa3Bg00730", "B", "HG-8", 3_634_882, 3_636_775, "-"),
    ("Fxa3Cg00714", "C", "HG-7", 3_588_659, 3_590_637, "-"),
    ("Fxa3Cg00715", "C", "HG-8", 3_593_118, 3_595_056, "-"),
    ("Fxa3Dg00712", "D", "HG-7", 3_531_234, 3_532_975, "-"),
    ("Fxa3Dg00713", "D", "HG-8", 3_536_882, 3_538_801, "-"),
]

GENES_CHR6 = [
    ("Fxa6Ag00747", "A", "HG-18", 4_392_830, 4_395_536, "-"),
    ("Fxa6Ag00748", "A", "HG-19", 4_399_329, 4_401_677, "+"),
    ("Fxa6Bg00716", "B", "HG-19", 3_824_767, 3_827_062, "+"),
    ("Fxa6Dg00623", "D", "HG-18", 3_530_505, 3_533_151, "-"),
    ("Fxa6Dg00624", "D", "HG-19", 3_533_995, 3_536_249, "+"),
]

# Distancias intergenicas (de la tabla de tandems-analysis.md)
INTERGENIC_CHR3 = {
    "A": 3_679,
    "B": 2_509,
    "C": 2_480,
    "D": 3_906,
}
INTERGENIC_CHR6 = {
    "A": 3_792,
    "D": 843,
}

# Subgenomas que faltan
MISSING_CHR6 = {
    "B": ["HG-18"],          # HG-18 no esta en B
    "C": ["HG-18", "HG-19"], # ninguno en C
}

# ---------------------------------------------------------------------------
# Sub-arbol Newick (clado NIP1) — copiado verbatim de tandems-analysis.md
# ---------------------------------------------------------------------------
NEWICK = """(
  (
    (MdNIP1_2:0.05590,MdNIP1_1:0.04220)99.8/1/100:0.17283,
    (
      (
        (
          (FaNIP1_1:0.00000,Fxa6Ag00747:0.00000)0/0.333/100:0.00000,
          Fxa6Dg00623:0.00995
        )100/1/100:0.23419,
        (
          (
            (Fxa3Ag00840:0.00479,Fxa3Dg00713:0.01437)85.8/0.999/94:0.01000,
            Fxa3Bg00730:0.01910
          )44.8/0.934/78:0.02888,
          Fxa3Cg00715:0.02508
        )99/1/100:0.08199
      )97.8/1/100:0.08861,
      (
        (
          (Fxa3Ag00839:0.00974,
            (Fxa3Dg00712:0.00975,Fxa3Cg00714:0.02951)0/0.333/40:0.00000
          )78/0.772/75:0.00467,
          Fxa3Bg00729:0.01248
        )83.2/0.985/83:0.01257,
        (Fxa6Ag00748:0.00998,
          (Fxa6Bg00716:0.00492,Fxa6Dg00624:0.00000)90.1/1/99:0.01443
        )94.9/1/99:0.02531
      )61.5/0.869/75:0.02958
    )95.8/1/98:0.08418
  )96.5/1/99:0.11607,
  (
    (HbNIP1_1:0.07703,HbNIP1_2:0.07995)99.7/1/100:0.19630,
    (
      (AtNIP1_1:0.23326,AtNIP1_2:0.05504)98.4/1/100:0.21827,
      (AtNIP2_1:0.70925,AtNIP3_1:0.77325)92.2/0.981/93:0.14929
    )95.8/1/97:0.16330
  )97.2/1/98:0.13947
);"""

# Mapeo gen -> HG / subgenoma para colorear hojas
LEAF_TO_HG = {}
LEAF_TO_SUBGENOME = {}
for gid, sg, hg, *_ in GENES_CHR3 + GENES_CHR6:
    LEAF_TO_HG[gid] = hg
    LEAF_TO_SUBGENOME[gid] = sg


# ---------------------------------------------------------------------------
# Construccion de la figura
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(20, 14), dpi=140)
gs = gridspec.GridSpec(2, 1, height_ratios=[1.0, 1.1], hspace=0.30)
ax_a = fig.add_subplot(gs[0])
ax_b = fig.add_subplot(gs[1])

# ===========================================================================
# PANEL A — Coordenadas genomicas
# ===========================================================================
# Layout: dos sub-bandas (chr_3 arriba, chr_6 abajo) cada una con 4 sub-pistas
# (A/B/C/D). Eje X = coordenada genomica en kb relativa al inicio del primer
# gen de cada subgenoma (alineamos los pares para visualizar la sintenia).
#
# Estrategia de alineamiento por subgenoma:
#   - Para cada subgenoma, restamos el start del primer gen para que todos
#     empiecen en x = 0 (coordenadas relativas). Asi los pares tandem en
#     los 4 subgenomas quedan visualmente comparables.
#   - El eje X esta en kb (dividimos por 1000).
#
# Ordenamos chr_3 y chr_6 en 8 filas (4 + 4); pero queremos UNA sub-banda
# por cromosoma con etiquetado tipo "chr_3" y "chr_6" a la izquierda. Asi
# que dibujamos 8 lineas, agrupadas con un separador horizontal.

# Coordenadas relativas (kb): start_relativo = start - min_start_subgenoma
def relativize(genes):
    """Devuelve dict subgenoma -> (offset_kb_func, lista_genes)."""
    by_sg = {}
    for gid, sg, hg, start, end, strand in genes:
        by_sg.setdefault(sg, []).append((gid, hg, start, end, strand))
    rel = {}
    for sg, items in by_sg.items():
        items.sort(key=lambda x: x[2])
        base_start = items[0][2]
        rel[sg] = [
            (gid, hg, (s - base_start) / 1000.0, (e - base_start) / 1000.0, strand)
            for (gid, hg, s, e, strand) in items
        ]
    return rel

rel_chr3 = relativize(GENES_CHR3)
rel_chr6 = relativize(GENES_CHR6)

# Configuracion de Y para Panel A
# Filas (de arriba a abajo en pantalla con invert_yaxis):
#   chr_3: A (y=7), B (y=6), C (y=5), D (y=4)
#   chr_6: A (y=2), B (y=1.5 si presente), C (y=1 si presente), D (y=0)
# Usamos espaciado constante.
row_height = 1.0
gap_between_chrs = 1.2  # gap visual entre chr_3 y chr_6

# Posicion Y de cada fila (mayor = mas arriba)
chr3_y = {"A": 7.0, "B": 6.0, "C": 5.0, "D": 4.0}
chr6_y = {"A": 2.0, "B": 1.0, "C": 0.0, "D": -1.0}
all_y = {**{(("chr_3", sg)): y for sg, y in chr3_y.items()},
         **{(("chr_6", sg)): y for sg, y in chr6_y.items()}}

# Escala X — usamos el rango mas amplio (chr_3) como referencia visual; cada
# punto en kb es proporcional. Maximo end (kb) entre todos los pares es ~8.0 kb
# (HG-7 ocupa 0-2.3 kb + intergenico 3.9 kb + HG-8 ~2 kb = ~8 kb).
xmax_kb = 10.0
ax_a.set_xlim(-0.6, xmax_kb)
ax_a.set_ylim(-2.0, 8.4)

# Dibujar lineas de fondo (cromosoma) para cada subgenoma
chr_line_color = "#999999"
chr_line_lw = 1.3
for sg, y in chr3_y.items():
    ax_a.plot([0, xmax_kb - 0.3], [y, y], color=chr_line_color,
              linewidth=chr_line_lw, zorder=1, solid_capstyle="round")
for sg, y in chr6_y.items():
    ax_a.plot([0, xmax_kb - 0.3], [y, y], color=chr_line_color,
              linewidth=chr_line_lw, zorder=1, solid_capstyle="round")

# Etiquetas a la izquierda: cromosoma (centrado) + subgenoma (cada fila)
# Etiqueta chr_3 centrada verticalmente entre las 4 filas chr_3
ax_a.text(-0.55, np.mean(list(chr3_y.values())), "chr_3",
          ha="right", va="center", fontsize=16, fontweight="bold",
          rotation=90, color="#222222")
ax_a.text(-0.55, np.mean(list(chr6_y.values())), "chr_6",
          ha="right", va="center", fontsize=16, fontweight="bold",
          rotation=90, color="#222222")

# Etiqueta Subgenoma X a la izquierda de cada fila (en x = -0.25)
for sg, y in chr3_y.items():
    ax_a.text(-0.25, y, f"Subgenoma {sg}", ha="right", va="center",
              fontsize=11, color=SG_COLORS[sg], fontweight="bold")
for sg, y in chr6_y.items():
    ax_a.text(-0.25, y, f"Subgenoma {sg}", ha="right", va="center",
              fontsize=11, color=SG_COLORS[sg], fontweight="bold")

# Funcion para dibujar un gen con flecha de orientacion
def draw_gene(ax, x_start_kb, x_end_kb, y, hg, strand, gene_id):
    """
    Rectangulo coloreado por HG + flecha de strand.
    x_start_kb < x_end_kb siempre (coordenadas genomicas).
    La flecha apunta a la derecha si strand == '+', a la izquierda si '-'.
    El cuerpo del gen ocupa start..end; la cabeza de flecha sobresale ~0.15 kb.
    """
    color = HG_COLORS[hg]
    height = 0.45  # altura del rectangulo
    width = x_end_kb - x_start_kb
    head_w = 0.18  # cabeza flecha (kb)

    if strand == "+":
        # Cuerpo desde x_start hasta x_end-head_w; cabeza hasta x_end
        body_w = max(width - head_w, width * 0.5)
        arrow = FancyArrow(
            x_start_kb, y, body_w, 0,
            width=height * 0.85,
            head_width=height * 1.05,
            head_length=head_w,
            length_includes_head=False,
            facecolor=color, edgecolor="black", linewidth=0.7,
            zorder=4,
        )
    else:
        # Strand '-': flecha apunta a la izquierda; dibujamos desde x_end hacia x_start
        body_w = max(width - head_w, width * 0.5)
        arrow = FancyArrow(
            x_end_kb, y, -body_w, 0,
            width=height * 0.85,
            head_width=height * 1.05,
            head_length=head_w,
            length_includes_head=False,
            facecolor=color, edgecolor="black", linewidth=0.7,
            zorder=4,
        )
    ax.add_patch(arrow)

    # Etiqueta del gen encima del rectangulo (fontsize pequeno)
    label_y = y + 0.45
    ax.text((x_start_kb + x_end_kb) / 2, label_y, gene_id,
            ha="center", va="bottom", fontsize=7.5,
            color="#222222", zorder=5,
            family="monospace")

# Dibujar genes chr_3
for sg, items in rel_chr3.items():
    y = chr3_y[sg]
    for gid, hg, start_kb, end_kb, strand in items:
        draw_gene(ax_a, start_kb, end_kb, y, hg, strand, gid)
    # Anotacion de distancia intergenica entre HG-7 (item[0]) y HG-8 (item[1])
    if len(items) >= 2:
        g1 = items[0]
        g2 = items[1]
        dist_kb = INTERGENIC_CHR3[sg] / 1000.0
        mid_x = (g1[3] + g2[2]) / 2  # entre end de HG-7 y start de HG-8
        ax_a.text(mid_x, y - 0.42, f"{dist_kb:.2f} kb".replace(".", ","),
                  ha="center", va="top", fontsize=8.5,
                  color="#444444", style="italic", zorder=5)

# Dibujar genes chr_6
for sg, items in rel_chr6.items():
    y = chr6_y[sg]
    for gid, hg, start_kb, end_kb, strand in items:
        draw_gene(ax_a, start_kb, end_kb, y, hg, strand, gid)
    # Anotacion de distancia intergenica solo si hay ambos genes
    if sg in INTERGENIC_CHR6:
        dist_kb = INTERGENIC_CHR6[sg] / 1000.0
        # Calcular el mid entre el HG-18 y el HG-19 (orden creciente)
        items_sorted = sorted(items, key=lambda x: x[2])  # por start_kb
        if len(items_sorted) >= 2:
            mid_x = (items_sorted[0][3] + items_sorted[1][2]) / 2
            ax_a.text(mid_x, y - 0.42, f"{dist_kb:.2f} kb".replace(".", ","),
                      ha="center", va="top", fontsize=8.5,
                      color="#444444", style="italic", zorder=5)

# Marcadores de copias ausentes en chr_6
# Subgenoma B: falta HG-18 (presente HG-19). HG-19 esta en chr_6B en strand +.
# Subgenoma C: faltan ambos.
def draw_missing(ax, x_center_kb, y, hg):
    """Caja vacia con simbolo '□' y etiqueta debajo."""
    color = HG_COLORS[hg]
    box_w = 1.6
    box_h = 0.42
    ax.add_patch(Rectangle(
        (x_center_kb - box_w/2, y - box_h/2),
        box_w, box_h,
        facecolor="white", edgecolor=color, linewidth=1.0,
        linestyle=(0, (3, 2)), zorder=4,
    ))
    ax.text(x_center_kb, y, "ausente", ha="center", va="center",
            fontsize=7.5, color=color, style="italic", zorder=5)
    ax.text(x_center_kb, y + 0.45, hg, ha="center", va="bottom",
            fontsize=7.5, color="#222222", zorder=5,
            family="monospace")

# Subgenoma B chr_6: pintar el HG-19 (Fxa6Bg00716) tal cual y marcar HG-18 ausente
# en el espacio "antes" del HG-19. Como rel_chr6["B"] tiene solo HG-19 alineado
# a x=0, ponemos el "ausente" en x = -? No, mejor: alineamos por la posicion
# del HG-18 del subgenoma A (referencia: HG-18 va antes que HG-19).
# Hacemos esto manualmente: en B el HG-19 esta en start_kb=0 (porque es el unico),
# pero deberia estar a la derecha del intergenico ausente. Re-alineamos:
# Para visualizacion comparativa, mostramos el HG-19 de B desplazado para que
# coincida verticalmente con el HG-19 de A/D.

# Re-alineacion de subgenoma B chr_6:
# HG-19 de A esta en (start_A - start_HG18_A)/1000 ≈ (4399329 - 4392830)/1000 = 6.499 kb
# Asi que ponemos el HG-19 de B en x_start = 6.499 kb (y end = 6.499 + len/1000)

# Borramos el dibujo previo del subgenoma B (lo redibujamos manualmente)
# OJO: ya lo hemos pintado en rel_chr6["B"], donde estaba alineado a x=0.
# Para evitarlo, recalculamos rel_chr6 con un offset fijo por subgenoma referido
# al primer gen del subgenoma A.

# Mejor estrategia: re-relativizar usando como ancla el HG-18 (el "izquierdo")
# de cada subgenoma. Para B y C, donde no hay HG-18, usamos como referencia
# la posicion teorica que tendria si estuviera (extrapolada de A/D).

# Recomputamos rel_chr6 con la nueva regla:
def relativize_chr6():
    """
    Alineamos chr_6 tomando como ancla la posicion del HG-18 (o, si falta,
    la posicion teorica antes del HG-19 igual a la distancia intergenica
    media observada).
    """
    by_sg = {sg: [] for sg in "ABCD"}
    for gid, sg, hg, start, end, strand in GENES_CHR6:
        by_sg[sg].append((gid, hg, start, end, strand))

    rel = {}
    for sg in "ABCD":
        items = sorted(by_sg[sg], key=lambda x: x[2])
        if not items:
            rel[sg] = []
            continue

        # Ancla: si hay HG-18, base_start = start de HG-18; si no, base_start
        # se calcula para que el HG-19 quede alineado con el HG-19 de A.
        hg18_items = [it for it in items if it[1] == "HG-18"]
        if hg18_items:
            base_start = hg18_items[0][2]
        else:
            # Solo HG-19. Asumimos que el HG-19 deberia comenzar en la misma
            # posicion relativa que en A: 6.499 kb desde base_start.
            # Buscamos el HG-19 y desplazamos.
            hg19_items = [it for it in items if it[1] == "HG-19"]
            assert hg19_items, f"Subgenoma {sg} sin HG-18 ni HG-19"
            base_start = hg19_items[0][2] - 6_499  # en bp
        rel[sg] = [
            (gid, hg, (s - base_start) / 1000.0, (e - base_start) / 1000.0, strand)
            for (gid, hg, s, e, strand) in items
        ]
    return rel

rel_chr6 = relativize_chr6()

# Repintamos chr_6: BORRAR previo. Como matplotlib no permite borrar facilmente
# selectivos, hemos pintado los chr_6 una vez ya. Necesitamos limpiar y rehacer.
# Solucion: redibujamos ax_a desde cero para chr_6.
# Por simplicidad, vamos a refactorizar: limpiamos ax_a y rehacemos todo en orden.

# === Refactor: limpiar ax_a y rehacer en el orden correcto ==================
ax_a.clear()
ax_a.set_xlim(-1.0, xmax_kb)
ax_a.set_ylim(-2.0, 8.4)
ax_a.invert_yaxis()  # mas natural: chr_3 arriba, chr_6 abajo
# Invertir Y reordena la lectura visual: mantener chr3 arriba (Y mayor en pantalla)
# Recolocar Y: tras invert_yaxis(), Y mayor queda abajo. Re-coloquemos sin invertir.
ax_a.invert_yaxis()  # vuelve al estado normal

# Reset y dibujar todo en orden:
ax_a.set_xlim(-1.0, xmax_kb)
ax_a.set_ylim(-2.0, 8.4)

# (1) Lineas de cromosoma de fondo
for sg, y in chr3_y.items():
    ax_a.plot([0, xmax_kb - 0.3], [y, y], color=chr_line_color,
              linewidth=chr_line_lw, zorder=1, solid_capstyle="round")
for sg, y in chr6_y.items():
    ax_a.plot([0, xmax_kb - 0.3], [y, y], color=chr_line_color,
              linewidth=chr_line_lw, zorder=1, solid_capstyle="round")

# (2) Etiquetas chr_3 / chr_6 a la izquierda (rotacion vertical)
ax_a.text(-0.95, np.mean(list(chr3_y.values())), "chr_3",
          ha="center", va="center", fontsize=16, fontweight="bold",
          rotation=90, color="#222222")
ax_a.text(-0.95, np.mean(list(chr6_y.values())), "chr_6",
          ha="center", va="center", fontsize=16, fontweight="bold",
          rotation=90, color="#222222")

# (3) Etiqueta de subgenoma por fila
for sg, y in chr3_y.items():
    ax_a.text(-0.55, y, f"Subg. {sg}", ha="center", va="center",
              fontsize=10, color=SG_COLORS[sg], fontweight="bold")
for sg, y in chr6_y.items():
    ax_a.text(-0.55, y, f"Subg. {sg}", ha="center", va="center",
              fontsize=10, color=SG_COLORS[sg], fontweight="bold")

# (4) Linea horizontal separadora entre chr_3 y chr_6
sep_y = (min(chr3_y.values()) + max(chr6_y.values())) / 2
ax_a.axhline(sep_y, color="#cccccc", linewidth=0.8, linestyle="--", zorder=0)

# (5) Dibujar genes chr_3
for sg, items in rel_chr3.items():
    y = chr3_y[sg]
    for gid, hg, start_kb, end_kb, strand in items:
        draw_gene(ax_a, start_kb, end_kb, y, hg, strand, gid)
    if len(items) >= 2:
        g1 = items[0]
        g2 = items[1]
        dist_kb = INTERGENIC_CHR3[sg] / 1000.0
        mid_x = (g1[3] + g2[2]) / 2
        ax_a.text(mid_x, y - 0.42, f"{dist_kb:.2f} kb".replace(".", ","),
                  ha="center", va="top", fontsize=8.5,
                  color="#444444", style="italic", zorder=5)

# (6) Dibujar genes chr_6 (re-relativizado)
for sg, items in rel_chr6.items():
    y = chr6_y[sg]
    items_sorted = sorted(items, key=lambda x: x[2])
    for gid, hg, start_kb, end_kb, strand in items_sorted:
        draw_gene(ax_a, start_kb, end_kb, y, hg, strand, gid)

    # Distancia intergenica
    if sg in INTERGENIC_CHR6 and len(items_sorted) >= 2:
        dist_kb = INTERGENIC_CHR6[sg] / 1000.0
        mid_x = (items_sorted[0][3] + items_sorted[1][2]) / 2
        ax_a.text(mid_x, y - 0.42, f"{dist_kb:.2f} kb".replace(".", ","),
                  ha="center", va="top", fontsize=8.5,
                  color="#444444", style="italic", zorder=5)

    # Marcadores de copias ausentes
    if sg in MISSING_CHR6:
        for hg_missing in MISSING_CHR6[sg]:
            # Posicion teorica del HG-missing: copiamos el rango de A
            ref_items = [it for it in rel_chr6["A"] if it[1] == hg_missing]
            if ref_items:
                gid_ref, _, s_ref, e_ref, _ = ref_items[0]
                x_center = (s_ref + e_ref) / 2
                draw_missing(ax_a, x_center, y, hg_missing)

# (7) Lineas verticales discontinuas conectando los pares tandem entre subgenomas
# Para chr_3: tres lineas verticales finas, una por el centro de cada HG (7, 8)
# atravesando las 4 filas A/B/C/D, alineando visualmente la sintenia.
def draw_vertical_synteny(ax, rel_dict, y_dict, hg_list):
    """Dibuja una linea discontinua vertical en el centro promedio de cada HG."""
    for hg in hg_list:
        xs = []
        ys = []
        for sg, items in rel_dict.items():
            for gid, _hg, s, e, _strand in items:
                if _hg == hg:
                    xs.append((s + e) / 2)
                    ys.append(y_dict[sg])
        if xs and ys:
            mean_x = np.mean(xs)
            y_min = min(ys) - 0.25
            y_max = max(ys) + 0.25
            ax.plot([mean_x, mean_x], [y_min, y_max],
                    color=HG_COLORS[hg], linewidth=0.7,
                    linestyle=(0, (2, 3)), alpha=0.55, zorder=0.5)

draw_vertical_synteny(ax_a, rel_chr3, chr3_y, ["HG-7", "HG-8"])
draw_vertical_synteny(ax_a, rel_chr6, chr6_y, ["HG-18", "HG-19"])

# (8) Eje X
ax_a.set_xlabel("Coordenada genomica relativa (kb, alineada por subgenoma)",
                fontsize=11)
ax_a.set_xticks(np.arange(0, xmax_kb + 0.1, 1.0))
ax_a.set_xticklabels([f"{int(t)}" for t in np.arange(0, xmax_kb + 0.1, 1.0)],
                     fontsize=9)
ax_a.tick_params(axis="x", length=3)
ax_a.set_yticks([])  # ocultamos eje Y (las filas se etiquetan a mano)

# Ocultar bordes superior/derecho/izquierdo
for spine in ["top", "right", "left"]:
    ax_a.spines[spine].set_visible(False)

# (9) Leyenda de Panel A: 4 HG + flechas
legend_a_handles = [
    Line2D([0], [0], marker="s", color="w", markerfacecolor=HG_COLORS["HG-7"],
           markeredgecolor="black", markersize=12, label="HG-7 (NIP1;1)"),
    Line2D([0], [0], marker="s", color="w", markerfacecolor=HG_COLORS["HG-8"],
           markeredgecolor="black", markersize=12, label="HG-8 (NIP1;1)"),
    Line2D([0], [0], marker="s", color="w", markerfacecolor=HG_COLORS["HG-18"],
           markeredgecolor="black", markersize=12, label="HG-18 (NIP1;2)"),
    Line2D([0], [0], marker="s", color="w", markerfacecolor=HG_COLORS["HG-19"],
           markeredgecolor="black", markersize=12, label="HG-19 (NIP1;2)"),
    Line2D([0], [0], marker=">", color="black", markerfacecolor="black",
           markersize=10, linewidth=0, label="Hebra + (sentido $\\rightarrow$)"),
    Line2D([0], [0], marker="<", color="black", markerfacecolor="black",
           markersize=10, linewidth=0, label="Hebra $-$ (sentido $\\leftarrow$)"),
    Line2D([0], [0], marker="s", color="w", markerfacecolor="white",
           markeredgecolor="#888888", markersize=12,
           linestyle=(0, (3, 2)), label="Copia ausente"),
]
ax_a.legend(handles=legend_a_handles, loc="upper left",
            bbox_to_anchor=(1.005, 1.0), fontsize=9.5,
            framealpha=0.95, borderaxespad=0)


# ===========================================================================
# PANEL B — Sub-arbol filogenetico del clado NIP1
# ===========================================================================
# Parseamos el Newick con Bio.Phylo y dibujamos manualmente para controlar
# colores de rama por HG y etiquetado de soportes.

tree = Phylo.read(io.StringIO(NEWICK), "newick")

# Asignar coordenadas X (longitud de rama acumulada) e Y (orden vertical de hojas)
# Bio.Phylo no expone directamente x/y; implementamos un layout simple cladograma
# proporcional a longitudes de rama.

# Recopilar hojas en orden de pre-order DFS
terminals = tree.get_terminals()
# Asignar y_pos a cada terminal en orden
y_positions = {t.name: i for i, t in enumerate(terminals)}

# Calcular x_pos (distancia desde la raiz) para cada nodo
def assign_x(node, x=0.0):
    branch_len = node.branch_length if node.branch_length else 0.0
    node._x = x + branch_len
    for child in node.clades:
        assign_x(child, node._x)

assign_x(tree.root, 0.0)

# Calcular y_pos para nodos internos (media de hojas descendientes)
def assign_y(node):
    if node.is_terminal():
        node._y = y_positions[node.name]
        return node._y
    ys = [assign_y(c) for c in node.clades]
    node._y = (min(ys) + max(ys)) / 2
    return node._y

assign_y(tree.root)

# Funcion para asignar color a una clade segun el HG que contiene (mayoria de hojas Fxa)
def clade_hg(node):
    """Devuelve el HG dominante de las hojas terminales bajo este nodo, o None."""
    terms = node.get_terminals()
    hgs = [LEAF_TO_HG.get(t.name) for t in terms if LEAF_TO_HG.get(t.name)]
    if not hgs:
        return None
    # Si todas las hojas Fxa pertenecen al mismo HG, devolvemos ese
    unique = set(hgs)
    if len(unique) == 1:
        return next(iter(unique))
    return None  # mezcla -> usamos color neutro

# Dibujar el arbol
def draw_tree(ax, node, parent_x=0.0):
    """Dibuja ramas recursivamente."""
    x = getattr(node, "_x", 0.0)
    y = getattr(node, "_y", 0.0)

    # Determinar color de la rama:
    # - Si toda la descendencia es Fxa de un solo HG -> color del HG
    # - Si es interna mezclando HGs Fxa -> color "ancestral" gris
    # - Si es outgroup (At, Hb, Md, FaNIP1_1) -> gris oscuro
    hg = clade_hg(node)
    if hg:
        color = HG_COLORS[hg]
        lw = 2.2
    else:
        # Comprobar si la rama es outgroup vs ancestral interno
        terms = [t.name for t in node.get_terminals()]
        fxa_terms = [t for t in terms if t in LEAF_TO_HG]
        if not fxa_terms:
            color = "#777777"  # outgroup puro
            lw = 1.2
        else:
            color = "#444444"  # interna mezclando Fxa
            lw = 1.5

    # Rama horizontal desde parent_x hasta x
    ax.plot([parent_x, x], [y, y], color=color, linewidth=lw,
            solid_capstyle="round", zorder=3)

    if not node.is_terminal():
        # Rama vertical conectando hijos
        ys_children = [c._y for c in node.clades]
        ax.plot([x, x], [min(ys_children), max(ys_children)],
                color=color, linewidth=lw,
                solid_capstyle="round", zorder=3)
        for c in node.clades:
            draw_tree(ax, c, parent_x=x)

draw_tree(ax_b, tree.root, parent_x=0.0)

# Etiquetas de hojas a la derecha del valor x_max + offset
x_max = max(t._x for t in terminals)
label_offset = 0.005

# Calcular xmax_b para ax_b limits (dejar espacio a las etiquetas)
ax_b_xmax = x_max + 0.55
ax_b.set_xlim(-0.02, ax_b_xmax)
ax_b.set_ylim(-1.5, len(terminals) + 0.5)
ax_b.invert_yaxis()  # hojas de arriba a abajo en orden natural

for t in terminals:
    name = t.name
    sg = LEAF_TO_SUBGENOME.get(name)
    hg = LEAF_TO_HG.get(name)

    if hg:
        # Hoja Fxa: color por subgenoma; sufijo " [HG-X, Subg. Y]"
        sg_color = SG_COLORS[sg]
        label_text = f"{name}  [{hg}, Subg. {sg}]"
        ax_b.text(t._x + label_offset, t._y, label_text,
                  ha="left", va="center", fontsize=8.5,
                  color=sg_color, family="monospace", zorder=4)
    else:
        # Hoja outgroup (FaNIP1_1, MdNIP1_*, HbNIP1_*, AtNIP*)
        label_text = name
        ax_b.text(t._x + label_offset, t._y, label_text,
                  ha="left", va="center", fontsize=8.5,
                  color="#555555", style="italic",
                  family="monospace", zorder=4)

# Anotar apoyos en las ramas internas criticas
# Bio.Phylo guarda los valores como Clade.confidence (cuando es float) o como
# Clade.name si la cadena no es parseable. En nuestro Newick los soportes son
# strings con el formato SH-aLRT/aBayes/UFboot.
SUPPORT_LABELS = {}  # mapear nodo -> string a mostrar

# Recorrer todos los nodos no-terminales y leer confidence/name
def collect_supports(node):
    if not node.is_terminal():
        # En Bio.Phylo, los labels de nodo interno (post ')') van a Clade.confidence
        # si son numericos, o a Clade.name si son strings tipo "0.5/0.8/95".
        # Probamos ambos.
        support_str = None
        if node.confidence is not None:
            support_str = str(node.confidence)
        elif node.name:
            support_str = node.name
        if support_str:
            SUPPORT_LABELS[id(node)] = support_str
        for c in node.clades:
            collect_supports(c)

collect_supports(tree.root)

# Filtrar y formatear: queremos mostrar solo las ramas con apoyos relevantes,
# no todas (para no saturar). Mostramos las del clado NIP1 interno con apoyo > 80
# o las criticas (MRCAs de HG-7/HG-8/HG-18/HG-19/HG-7+HG-19/HG-8+HG-18).

CRITICAL_SUPPORTS = {
    "95.8/1/98",     # MRCA HG-7+HG-19
    "97.8/1/100",    # MRCA HG-8+HG-18+FaNIP1_1
    "100/1/100",     # MRCA HG-18
    "99/1/100",      # MRCA HG-8
    "83.2/0.985/83", # MRCA HG-7
    "94.9/1/99",     # MRCA HG-19
    "61.5/0.869/75", # union HG-7/HG-19
    "96.5/1/99",     # estrella NIP1 Md+Fxa
}

def format_support(s):
    """Convierte '95.8/1/98' a '96/1/98' (sin decimales para SH-aLRT, sin
    cambio en UFboot, aBayes redondeado a 2 decimales si <1)."""
    parts = s.split("/")
    if len(parts) != 3:
        return s
    sh, aB, uf = parts
    try:
        sh_n = float(sh)
        sh_str = f"{sh_n:.0f}"
    except ValueError:
        sh_str = sh
    try:
        aB_n = float(aB)
        aB_str = f"{aB_n:.0f}" if aB_n == 1 else f"{aB_n:.2f}".rstrip("0").rstrip(".")
    except ValueError:
        aB_str = aB
    try:
        uf_n = float(uf)
        uf_str = f"{uf_n:.0f}"
    except ValueError:
        uf_str = uf
    return f"{sh_str}/{aB_str}/{uf_str}"

# Dibujar etiquetas de apoyo
def annotate_supports(node, parent_x=0.0):
    if not node.is_terminal():
        sup = SUPPORT_LABELS.get(id(node))
        if sup and sup in CRITICAL_SUPPORTS:
            # Coordenadas: a la izquierda del nodo (cerca de parent), un poco arriba
            x_label = (parent_x + node._x) / 2
            y_label = node._y - 0.35
            fmt = format_support(sup)
            ax_b.text(x_label, y_label, fmt, ha="center", va="center",
                      fontsize=7.5, color="#222222",
                      bbox=dict(boxstyle="round,pad=0.18", facecolor="white",
                                edgecolor="#aaaaaa", linewidth=0.5, alpha=0.92),
                      zorder=5)
        for c in node.clades:
            annotate_supports(c, parent_x=node._x)

annotate_supports(tree.root, parent_x=0.0)

# Cuadro de anotacion lateral con topologia retirado por dos motivos:
# (i) solapaba con la rama del outgroup MdNIP1_1 del arbol.
# (ii) contenia jerga interna del pipeline (Fase 2 chromosome splitting,
#      outgroup_between, MAX_OUTGROUP_INSERTIONS) que no procede en una
#      figura de TFG. La topologia de pares hermanos (HG-7/HG-19 y
#      HG-8/HG-18) queda visible directamente en el arbol y se discute
#      en el cuerpo del 6.3.3.

# Barra de escala filogenetica (sub/sitio)
scale_x0 = ax_b_xmax - 0.20
scale_x1 = scale_x0 + 0.10
scale_y = len(terminals) - 0.5
ax_b.plot([scale_x0, scale_x1], [scale_y, scale_y], color="black", linewidth=1.5)
ax_b.text((scale_x0 + scale_x1) / 2, scale_y - 0.3,
          "0,10 sub/sitio".replace("0,10", "0,10"),
          ha="center", va="top", fontsize=8.5, color="#222222")

# Ejes Panel B
ax_b.set_xlabel("Distancia patristica (substituciones esperadas/sitio)",
                fontsize=11)
ax_b.set_yticks([])
for spine in ["top", "right", "left"]:
    ax_b.spines[spine].set_visible(False)
ax_b.tick_params(axis="x", length=3)

# Leyenda Panel B: HG colors + outgroup
legend_b_handles = [
    Line2D([0], [0], color=HG_COLORS["HG-7"], linewidth=2.5, label="Rama HG-7"),
    Line2D([0], [0], color=HG_COLORS["HG-8"], linewidth=2.5, label="Rama HG-8"),
    Line2D([0], [0], color=HG_COLORS["HG-18"], linewidth=2.5, label="Rama HG-18"),
    Line2D([0], [0], color=HG_COLORS["HG-19"], linewidth=2.5, label="Rama HG-19"),
    Line2D([0], [0], color="#444444", linewidth=2.0, label="Rama interna ancestral"),
    Line2D([0], [0], color="#777777", linewidth=1.5, label="Rama outgroup"),
]
# Subgenomas
for sg in "ABCD":
    legend_b_handles.append(
        Line2D([0], [0], marker="s", color="w",
               markerfacecolor=SG_COLORS[sg], markeredgecolor="black",
               markersize=10, label=f"Etiqueta Subg. {sg}")
    )
ax_b.legend(handles=legend_b_handles, loc="upper left",
            bbox_to_anchor=(1.005, 1.0), fontsize=9, framealpha=0.95,
            borderaxespad=0)

# Ajustes finales
fig.subplots_adjust(left=0.06, right=0.85, top=0.97, bottom=0.05, hspace=0.30)

# ---------------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------------
out_dir = config.ensure_results() / "figuras_rnaseq"
out_dir.mkdir(parents=True, exist_ok=True)

out_pdf = out_dir / "figura_tandems_schema.pdf"
out_png = out_dir / "figura_tandems_schema.png"

fig.savefig(str(out_pdf), dpi=200, bbox_inches="tight")
fig.savefig(str(out_png), dpi=200, bbox_inches="tight")
print("PDF:", out_pdf, out_pdf.stat().st_size)
print("PNG:", out_png, out_png.stat().st_size)
