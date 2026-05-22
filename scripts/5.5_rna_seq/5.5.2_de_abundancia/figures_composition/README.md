# Scripts de composición de figuras del apartado 6.3 (RNA-seq) del TFG

Scripts en Python (matplotlib + scipy + adjustText/PyMuPDF) que componen las figuras
del apartado 6.3 a partir de las matrices generadas por el pipeline de RNA-seq
(`08_basal_expression.R`, `07_de_analysis.R`, `13_homeolog_expression.R`). Leen sus
inputs de `data/rna_seq/` a través de `scripts/common/config.py`.

Todos implementan **Pirate plots** (Phillips, 2017): violín bilateral de densidad
kernel + puntos crudos con jitter + línea de la media, con la paleta canónica del
proyecto. Las figuras carecen de título embebido y las multipanel llevan etiquetas
(A)/(B) (la descripción va en el pie de figura, norma APA/UCAM).

## Inventario

| Script | Genera | Figura del TFG v11 | Inputs (`data/rna_seq/`) |
|---|---|---|---|
| `compose_fig6_basal_subfamilia.py` | Pirate plot de TPM por gen × subfamilia × tejido + columna TOTAL | **Figura 9** | `basal/basal_aquaporins_summary.csv` |
| `compose_fig_de_subfamilia.py` | Pirate plot de log2FC por subfamilia × tejido (inducido/reprimido) | **Figura 10** | `de/de_aquaporins_leaf.csv`, `de/de_aquaporins_roots.csv` |
| `compose_fig_homeologos_basal.py` | 2 paneles: TPM colapsado por grupo HG (A) + dominancia subgenómica (B) | **Figura 11** | `homeologos/collapsed_tpm.csv`, `dominance_overall.csv`, `homeolog_groups_summary.tsv`, `design_basal.csv` |
| `compose_fig_tandems_schema.py` | 2 paneles: localización genómica de los tándems NIP1 (A) + subárbol filogenético (B) | **Figura 12** | datos embebidos (coordenadas + Newick del subárbol) |

## Ejecución

Lo más cómodo es el runner global, que ejecuta estos y el resto de scripts de figuras
y reúne las salidas numeradas en `results/figuras_TFG/`:

```bash
python scripts/regenerar_figuras.py
```

O cada script por separado (lee de `data/` vía `config.py`, escribe en `results/figuras_rnaseq/`):

```bash
pip install pandas matplotlib scipy adjustText PyMuPDF
python scripts/5.5_rna_seq/5.5.2_de_abundancia/figures_composition/compose_fig6_basal_subfamilia.py
```

## Paleta de colores

```python
sf_colors = {"PIP": "#E74C3C", "TIP": "#3498DB", "NIP": "#2ECC71",
             "SIP": "#F39C12", "XIP": "#9B59B6"}
sg_colors = {"A": "#6D4C41", "B": "#C2185B", "C": "#00838F", "D": "#455A64"}
```
