# Anexo G — Figuras del cuerpo del TFG en alta resolución

Versiones en alta resolución (PNG 300 dpi + PDF vectorial donde aplica) de las
figuras del cuerpo del TFG (**v11**) que se generan por código en este
repositorio. Los nombres se corresponden **uno a uno con la numeración del TFG v11**.

> En el TFG v11 las Figuras 1–3 son ilustraciones tomadas de otras publicaciones
> (introducción/marco teórico) y la Figura 13 es una captura del visor eFP; ninguna
> de ellas se genera por código y por eso no están en este anexo. Las figuras de
> resultados generadas en este trabajo son las **Figuras 4 a 12**.

| Archivo | Figura del TFG v11 | Apartado | Script generador |
|---|---|---|---|
| `Anexo_G_Figura_4.png/.pdf` | Figura 4 — Distribución de longitudes de acuaporinas en Rosaceae (n=419) | 6.1 | `scripts/5.2_5.3_homologia_curacion/generar_visualizaciones_tfg.py` |
| `Anexo_G_Figura_5.png/.pdf` | Figura 5 — Identidad peptídica GFF3 vs Exonerate | 6.1 | `comparacion_identidad_gff3_exonerate.py` |
| `Anexo_G_Figura_6.png` | Figura 6 — PCA de variables fisicoquímicas (54,1 %) | 6.1 | `profiling_final_integrated.py` |
| `Anexo_G_Figura_8.png` | Figura 8 — Frecuencia de los 15 motivos MEME por sub-subfamilia | 6.2.3 | `analisis_motivos_final.py` |
| `Anexo_G_Figura_9.png/.pdf` | Figura 9 — TPM por gen y tejido por subfamilia (basal) | 6.3.1 | `compose_fig6_basal_subfamilia.py` |
| `Anexo_G_Figura_10.png/.pdf` | Figura 10 — Expresión diferencial (hoja/raíz) bajo estrés hídrico | 6.3.2 | `compose_fig_de_subfamilia.py` |
| `Anexo_G_Figura_11.png/.pdf` | Figura 11 — Grupos homeólogos: TPM colapsado y dominancia de subgenoma | 6.3.3 | `compose_fig_homeologos_basal.py` |
| `Anexo_G_Figura_12.png/.pdf` | Figura 12 — Duplicaciones tándem en el clado NIP1 | 6.3.3 | `compose_fig_tandems_schema.py` |

## Figuras del TFG no incluidas aquí (no se generan por código en el repo)

- **Figuras 1–3**: ilustraciones de otras publicaciones (introducción/marco teórico).
- **Figura 7** (árbol filogenético + mapa de calor): visualizada y anotada en **iTOL** (web), no es salida de un script del repositorio. Los soportes por nodo del árbol están en el **Anexo E**; el árbol y el alineamiento, en `data/filogenia/`.
- **Figura 13** (visor eFP de grupos homeólogos): captura de un HTML interactivo. El visor está en el **Anexo I** y publicado como página web (ver `docs/GITHUB_PAGES.md`).

Todas las figuras 4–12 se regeneran ejecutando sus scripts (leen de `data/` vía
`scripts/common/config.py`), carecen de título embebido y las multipanel llevan
etiquetas (A)/(B), conforme a la norma APA/UCAM (el título va en el pie de figura).
