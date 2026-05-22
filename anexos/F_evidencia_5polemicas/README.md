# Anexo F — Evidencia visual del descarte de las secuencias polémicas

Citado en §6.1 del TFG (curado). Reúne la evidencia numérica y visual del descarte de las **5 secuencias "polémicas"**: las que conservan 6 (o casi 6) hélices transmembrana nominales y, por tanto, requieren más justificación que un simple recuento topológico.

## Contenido

| Archivo | Contenido |
|---|---|
| `Anexo_F_memoria_pruebas_descarte.md` | Memoria técnica: tabla resumen, posición de los motivos NPA en el alineamiento, métricas por secuencia y umbrales estadísticos. |
| `Anexo_F_TablaY1_metricas_5_polemicas.csv` | Métricas por secuencia (TMHs, NPA-B/E, longitud, rama terminal, soportes, Mahalanobis, cobertura). |
| `Anexo_F_alineamiento_MAFFT_126_secs.fasta` | Alineamiento MAFFT E-INS-i de las 121 funcionales + 5 polémicas (126 secuencias). |
| `Anexo_F_dataset_126_sin_alinear.fasta` · `Anexo_F_solo_5_polemicas.fasta` | FASTAs de partida. |
| `Anexo_F_MAFFT_analysis_raw.txt` | Salida en bruto del alineamiento. |
| `Anexo_F_fig01..07_*.png` | Figuras del alineamiento: pérdida de NPA-B (fig01) y NPA-E (fig02), hueco central (fig03), N-terminal truncado (fig04) y extendido (fig05), gaps internos (fig06) y visión global (fig07). |

## Las 5 secuencias y su causa principal de descarte

| Secuencia | Causa principal |
|---|---|
| Fxa6Bg00715 | Deleción interna de 91 aa entre los motivos NPA-B y NPA-E. |
| Fxa6Cg01391 | Clado incongruente en el árbol + rama terminal ~10× la mediana + gaps internos. |
| Fxa5Bg03706 | Pérdida total del motivo NPA-E + truncamiento N-terminal + outlier del PCA. |
| Fxa3Ag00841 | Topología no canónica (8 TMHs) por extensión N-terminal de 104 aa. |
| Fxa5Ag03930 | Pérdida total del motivo NPA-B + truncamiento severo + outlier del PCA. |

Cada secuencia falla al menos un criterio cuantitativo objetivo trazable a un archivo de datos; dos de ellas pierden literalmente uno de los dos motivos NPA canónicos, lo que las descalifica como canales funcionales.
