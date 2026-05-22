# Anexo G — Figuras del cuerpo del TFG en alta resolución

Versiones en alta resolución (PNG 300 dpi + PDF vectorial donde aplica) de las
figuras del cuerpo del TFG que se generan por código en este repositorio. Los
nombres se corresponden **uno a uno con la numeración del TFG**.

| Archivo | Figura del TFG | Apartado | Script generador |
|---|---|---|---|
| `Anexo_G_Figura_1.png/.pdf` | Figura 1 — Distribución de longitudes de acuaporinas en Rosaceae (n=419) | 6.1 | `scripts/5.2_5.3_homologia_curacion/generar_visualizaciones_tfg.py` |
| `Anexo_G_Figura_2.png/.pdf` | Figura 2 — Identidad peptídica GFF3 vs Exonerate | 6.1 | `comparacion_identidad_gff3_exonerate.py` |
| `Anexo_G_Figura_3.png` | Figura 3 — PCA de variables fisicoquímicas (54,1 %) | 6.1 | `profiling_final_integrated.py` |
| `Anexo_G_Figura_5.png` | Figura 5 — Frecuencia de los 15 motivos MEME por sub-subfamilia | 6.2.3 | `analisis_motivos_final.py` |
| `Anexo_G_Figura_6.png/.pdf` | Figura 6 — TPM por gen y tejido por subfamilia (basal) | 6.3.1 | `compose_fig6_basal_subfamilia.py` |
| `Anexo_G_Figura_7.png/.pdf` | Figura 7 — Expresión diferencial (hoja/raíz) bajo estrés hídrico | 6.3.2 | `compose_fig_de_subfamilia.py` |
| `Anexo_G_Figura_8.png/.pdf` | Figura 8 — Grupos homeólogos: TPM colapsado y dominancia de subgenoma | 6.3.3 | `compose_fig_homeologos_basal.py` |
| `Anexo_G_Figura_9.png/.pdf` | Figura 9 — Duplicaciones tándem en el clado NIP1 | 6.3.3 | `compose_fig_tandems_schema.py` |

## Figuras del TFG no incluidas aquí (no se generan por código en el repo)

- **Figura 4** (árbol filogenético + mapa de calor): visualizada y anotada en **iTOL** (web), no es salida de un script del repositorio. Los soportes por nodo del árbol están en el **Anexo E**; el árbol y el alineamiento, en `data/filogenia/`.
- **Figura 10** (visor eFP de grupos homeólogos): es un HTML interactivo, no una imagen estática. Está en el **Anexo I** y publicado como página web (ver `docs/GITHUB_PAGES.md`).

Todas las figuras se regeneran ejecutando sus scripts (leen de `data/` vía
`scripts/common/config.py`) y carecen de título embebido; las multipanel llevan
etiquetas (A)/(B), conforme a la norma APA/UCAM (el título va en el pie de figura).
