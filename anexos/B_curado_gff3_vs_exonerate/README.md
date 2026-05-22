# Anexo B — Curado: GFF3 vs Exonerate (material para reuso y actualización oficial)

Citado en §6.1 del TFG (curado de secuencias) y referenciable desde §5.3 (Curación).

## Por qué este anexo existe

El curado del acuaporinoma decide, gen a gen, qué anotación adoptar entre la oficial del GDR (GFF3 'Benihoppe' v1) y la alternativa generada por Exonerate sobre proteínas de referencia. La decisión usa criterios estructurales (metionina inicial, 6 hélices transmembrana, motivos NPA conservados). Este anexo contiene **todos los datos intermedios** que sustentan esa decisión, de modo que:

- Cualquier persona puede inspeccionar y verificar el material en bruto que respalda cada decisión.
- Otros grupos que quieran **actualizar la anotación oficial de** *F.* x *ananassa* **'Benihoppe' en GDR** con los modelos curados disponen del GFF3 listo para sustituir (`Anexo_B_consenso_aqp_fixed.gff3`) y del FASTA del consenso (`Anexo_B_consenso_aqp_FINAL.fasta`).

## Contenido

### Inputs de la comparación

| Archivo | Contenido |
|---|---|
| `Anexo_B_aquaporin_peptides_GFF3.fasta` | FASTA peptídico de las 129 acuaporinas extraídas del GFF3 oficial 'Benihoppe' v1. |
| `Anexo_B_exonerate_genes_aqp.fasta` | FASTA peptídico de las anotaciones alternativas generadas por Exonerate protein2genome. |
| `Anexo_B_dataset_fusionado_GFF3_Exonerate.fasta` | FASTA combinado (GFF3 + Exonerate) usado como input de los análisis posteriores. |

### Análisis estructurales por fuente

| Archivo | Herramienta |
|---|---|
| `Anexo_B_DeepTMHMM_GFF3.3line` / `Anexo_B_DeepTMHMM_Exonerate.3line` | Predicción de topología transmembrana (DeepTMHMM 1.0). |
| `Anexo_B_Pepstats_GFF3.txt` / `Anexo_B_Pepstats_Exonerate.txt` | Estadísticos peptídicos (Pepstats EMBOSS): pI, carga, %hidrofobicidad, etc. |
| `Anexo_B_DeepLoc_GFF3.csv` / `Anexo_B_DeepLoc_Exonerate.csv` | Predicción de localización subcelular (DeepLoc 2.0). |

### Motivos peptídicos

| Archivo | Contenido |
|---|---|
| `Anexo_B_MEME_combined_15_motivos.txt` | Output de MEME sobre el dataset fusionado (15 motivos M1-M15 incluyendo M2 = NPA loop B). |

### Decisión y artefactos finales

| Archivo | Contenido |
|---|---|
| `Anexo_B_tabla_veredictos_GFF3_vs_Exonerate.tabular` | Tabla maestra: por gen, veredicto Fase 2 (PEP=GFF3 preferido / EXO=Exonerate preferido / MAL=ambas fallan), métricas comparativas, fuente seleccionada. |
| `Anexo_B_identidad_GFF3_vs_Exonerate.png` / `.pdf` | Figura de validación (Figura 2 del TFG): % de identidad peptídica GFF3 vs Exonerate por gen (alineamiento global BLOSUM62). 91 idénticos (100 %) / 38 discrepantes (30–99 %). Reproducible con `scripts/5.2_5.3_homologia_curacion/comparacion_identidad_gff3_exonerate.py` desde `datos/curado/`. |
| `Anexo_B_consenso_aqp_FINAL.fasta` | FASTA peptídico del consenso curado (129 secuencias = mejor opción por gen). |
| `Anexo_B_consenso_aqp_fixed.gff3` | **GFF3 corregido** con los 20 modelos Exonerate-curados sustituyendo a sus contrapartes GFF3 defectuosas. Listo para mergear con el GFF3 oficial completo del genoma. |

## Cómo actualizar el GFF3 oficial de 'Benihoppe' con estos modelos curados

1. Descargar el GFF3 oficial actual desde [GDR](https://www.rosaceae.org/species/fragaria/fragaria_x_ananassa).
2. Identificar las entradas correspondientes a las 20 acuaporinas Exonerate-curadas en `Anexo_B_tabla_veredictos_GFF3_vs_Exonerate.tabular` (filas con `fuente_seq=EXONERATE`).
3. Sustituir esas entradas del GFF3 oficial por las correspondientes de `Anexo_B_consenso_aqp_fixed.gff3`.
4. Ejecutar `scripts/5.2_5.3_homologia_curacion/12_substitute_gff3.py` (requiere ajustar paths a las versiones nuevas) para automatizar la sustitución.

## Trazabilidad

El veredicto de cada gen se basa en:
- **Metionina inicial** presente (criterio Met).
- **6 hélices transmembrana** exactas (criterio TMH=6).
- **Motivos NPA** conservados en los bucles B y E.

Si las tres condiciones se cumplen en GFF3, se prefiere GFF3 (`PEP`). Si fallan en GFF3 pero se cumplen en Exonerate, se prefiere Exonerate (`EXO`). Si fallan en ambas, el gen queda como descarte (`MAL`, ver `anexos/E_descartes/`).
