# Anexo I — Datos RNA-seq (matrices y control de calidad)

Matrices de expresión derivadas del análisis transcriptómico (apartado 6.3 del TFG)
e informe de control de calidad. Son los datos que sustentan las **Figuras 6, 7, 8 y 9**.

22 muestras *paired-end*, 6 tejidos de *F.* x *ananassa* 'Benihoppe' (la muestra
RootsCtrl_2, SRR30146487, se excluyó por PCA). Las matrices están restringidas a las
**121 acuaporinas funcionales** salvo donde se indique.

## Contenido

### Expresión basal por tejido
| Archivo | Contenido |
|---|---|
| `Anexo_I_basal_aquaporins_tpm.csv` | TPM por gen y muestra. |
| `Anexo_I_basal_aquaporins_normalized.csv` | Recuentos normalizados (DESeq2). |
| `Anexo_I_basal_aquaporins_detection.csv` | Detección por tejido (umbral TPM > 1). |
| `Anexo_I_basal_aquaporins_summary.csv` | Resumen por gen × tejido (media, sd, n). |

### Expresión diferencial (control vs estrés hídrico)
| Archivo | Contenido |
|---|---|
| `Anexo_I_de_aquaporins_leaf.csv` | DE en hoja (3 ctrl vs 3 estrés), DESeq2. |
| `Anexo_I_de_aquaporins_roots.csv` | DE en raíz (2 ctrl vs 3 estrés), DESeq2. |

### Grupos homeólogos
| Archivo | Contenido |
|---|---|
| `Anexo_I_homeolog_groups.tsv` / `_groups_summary.tsv` | 32 grupos homeólogos y su resumen. |
| `Anexo_I_homeolog_collapsed_tpm.csv` / `_collapsed_counts.csv` | TPM/recuentos colapsados por grupo. |
| `Anexo_I_homeolog_dominance_overall.csv` / `_dominance_by_tissue.csv` | Dominancia de subgenoma global y por tejido. |
| `Anexo_I_homeolog_dominant_subgenome.csv` / `_summary_statistics.csv` | Subgenoma dominante por grupo y estadísticos. |

### Diseño experimental y control de calidad
| Archivo | Contenido |
|---|---|
| `Anexo_I_design_basal.csv` | Tabla de diseño (muestra → tejido/condición). |
| `Anexo_I_MultiQC_report.html` | Informe MultiQC integrado (QC + asignación featureCounts de las 22 muestras). |
| `Anexo_I_MultiQC_general_stats.txt` / `_featurecounts.txt` | Tablas de estadísticas que resume el informe. |

Las mismas matrices (sin el prefijo `Anexo_I_`) están en `data/rna_seq/` del repo, desde
donde las leen los scripts de RNA-seq vía `scripts/common/config.py`.
