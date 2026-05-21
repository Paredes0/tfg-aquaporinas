# Scripts de composición de figuras del apartado 6.3 (RNA-seq) del TFG

Scripts en Python (matplotlib + scipy + PyMuPDF) que componen las figuras finales del apartado 6.3 a partir de los PDFs/CSVs generados por el pipeline principal (`08_basal_expression.R`, `07_de_analysis.R`, `13_homeolog_expression.R`).

Todos los scripts implementan **Pirate plots** (Phillips 2017): violín bilateral de densidad kernel + puntos crudos con jitter + línea horizontal de la media, con la paleta canónica del proyecto documentada en `memory/feedback_tfg_color_palette.md`.

## Inventario

| Script | Genera | Apartado del TFG | Dependencias del pipeline |
|---|---|---|---|
| `compose_fig6_basal_subfamilia.py` | `figura6_perfiles_subfamilia.{pdf,png}` — Pirate plot de TPM por gen × subfamilia × tejido + columna TOTAL | 6.3.1 | `basal_aquaporins_summary.csv` |
| `compose_fig7_validacion_pca.py` | `figura7_validacion.{pdf,png}` — composición lado a lado del PCA + matriz de correlación | 6.3.1 | `pca_aquaporins.pdf` + `correlation_samples.pdf` |
| `compose_fig_de_subfamilia.py` | `figura_de_stripplot.{pdf,png}` — Pirate plot de log2FC por subfamilia × tejido con sombreado direccional (inducido/reprimido) | 6.3.2 | `de_aquaporins_leaf.csv`, `de_aquaporins_roots.csv` |
| `compose_fig_homeologos_basal.py` | `figura_homeologos_basal.{pdf,png}` — 2 paneles: TPM colapsado por grupo HG (Panel A) + dominancia subgenómica (Panel B) | 6.3.3 | `collapsed_tpm.csv`, `dominance_overall.csv`, `homeolog_groups_summary.tsv` |

## Ejecución

Los scripts asumen que el pipeline principal ya se ha ejecutado y que los archivos del directorio `Z:/work/RNA-seq_test/results/` están presentes. Para ejecutar via SSH al PC Linux remoto (entorno conda `rnaseq_aqp`):

```bash
ssh noe@192.168.18.7 "/home/noe/micromamba/envs/rnaseq_aqp/bin/python /tmp/compose_fig6_basal_subfamilia.py"
```

Para ejecutar localmente con un entorno equivalente:

```bash
# Dependencias mínimas
pip install pandas matplotlib scipy PyMuPDF

python compose_fig6_basal_subfamilia.py
```

Los outputs (.pdf y .png) se guardan en el mismo directorio del que provienen los inputs y se copian al anexo `ANEXOS_TFG_DATOS/AnexoA_RNAseq_datos/figuras_compuestas/` para entrega.

## Paleta de colores

Definida en `memory/feedback_tfg_color_palette.md`. En código:

```python
sf_colors = {"PIP": "#E74C3C", "TIP": "#3498DB", "NIP": "#2ECC71",
             "SIP": "#F39C12", "XIP": "#9B59B6"}
sg_colors = {"A": "#6D4C41", "B": "#C2185B", "C": "#00838F", "D": "#455A64"}
```

## Nomenclatura del gráfico

Las figuras 6, DE y homeólogos basal usan **Pirate plot** (Phillips 2017, paquete `yarrr` de R) — variante de strip plot con violín de densidad + puntos crudos + media superpuesta. Aceptado en literatura biomédica como representación rica de distribuciones con N pequeño.
