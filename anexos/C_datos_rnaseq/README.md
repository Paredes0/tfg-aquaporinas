# Anexo C — Datos RNA-seq (matrices y control de calidad)

Matrices de expresión derivadas del análisis transcriptómico (apartado 6.3 del TFG)
e informe de control de calidad. Son los datos que sustentan las **Figuras 6, 7, 8 y 9**.

22 muestras *paired-end*, 6 tejidos de *F.* x *ananassa* 'Benihoppe' (la muestra
RootsCtrl_2, SRR30146487, se excluyó por PCA; evidencia en
`Anexo_C_PCA_outlier_RootsCtrl2.pdf`/`.png`). Las matrices están restringidas a las
**121 acuaporinas funcionales** salvo donde se indique.

## Contenido

### Expresión basal por tejido
| Archivo | Contenido |
|---|---|
| `Anexo_C_basal_aquaporins_tpm.csv` | TPM por gen y muestra. |
| `Anexo_C_basal_aquaporins_normalized.csv` | Recuentos normalizados (DESeq2). |
| `Anexo_C_basal_aquaporins_detection.csv` | Detección por tejido (umbral TPM > 1). |
| `Anexo_C_basal_aquaporins_summary.csv` | Resumen por gen × tejido (media, sd, n). |

### Expresión diferencial (control vs estrés hídrico)
| Archivo | Contenido |
|---|---|
| `Anexo_C_de_aquaporins_leaf.csv` | DE en hoja (3 ctrl vs 3 estrés), DESeq2. |
| `Anexo_C_de_aquaporins_roots.csv` | DE en raíz (2 ctrl vs 3 estrés), DESeq2. |

### Grupos homeólogos
| Archivo | Contenido |
|---|---|
| `Anexo_C_homeolog_groups.tsv` / `_groups_summary.tsv` | 32 grupos homeólogos y su resumen. |
| `Anexo_C_homeolog_collapsed_tpm.csv` / `_collapsed_counts.csv` | TPM/recuentos colapsados por grupo. |
| `Anexo_C_homeolog_dominance_overall.csv` / `_dominance_by_tissue.csv` | Dominancia de subgenoma global y por tejido. |
| `Anexo_C_homeolog_dominant_subgenome.csv` / `_summary_statistics.csv` | Subgenoma dominante por grupo y estadísticos. |

### Diseño experimental y control de calidad
| Archivo | Contenido |
|---|---|
| `Anexo_C_design_basal.csv` | Tabla de diseño (muestra → tejido/condición). |
| `Anexo_C_MultiQC_report.html` | Informe MultiQC integrado (QC + asignación featureCounts de las 21 muestras retenidas; RootsCtrl_2 no consta, pues se generó tras su exclusión). |
| `Anexo_C_MultiQC_general_stats.txt` / `_featurecounts.txt` | Tablas de estadísticas que resume el informe. |
| `Anexo_C_PCA_outlier_RootsCtrl2.pdf` / `.png` | PCA genómico (VST, contraste de raíz) que muestra que la réplica control RootsCtrl_2 (SRR30146487) no agrupa con sus réplicas biológicas (RootsCtrl_1 y RootsCtrl_3); sustenta su exclusión descrita en §5.5.2. Generado en modo `include` reincorporando la muestra purgada. |

Las mismas matrices (sin el prefijo `Anexo_C_`) están en `datos/rna_seq/` del repo, desde
donde las leen los scripts de RNA-seq vía `scripts/common/config.py`.
