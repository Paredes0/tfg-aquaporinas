# Anexo B — Detalle de las 26 secuencias excluidas

Citado en §6.1 del TFG.

## Contenido

- `Anexo_B_tabla_26_descartes.csv` — 26 filas x 10 columnas con ID, veredicto Fase 2, TMHs (GFF3/Exonerate), motivos NPA, longitud, filtro principal, métrica complementaria y expresión RNA-seq.

## Composición

- **17 descartes evidentes** por filtro topológico estricto (3 entradas no-membrana + 14 con TMHs <= 5).
- **6 polémicas** analizadas en detalle en Anexo C (Mahalanobis) y Anexo F (evidencia visual).
- **3 candidatas a reanotación futura** con expresión RNA-seq detectable.

## Trazabilidad

Datos derivados de cruzar `tabla_Aquaporinas_traduccion.tabular` (veredictos), `reannotation_candidates.tsv` (expresión) y `Genes_eliminados_NO_aqp.csv` (BLAST de no-AQP).
