# Workflows Galaxy del TFG

Dos workflows usados respectivamente en §5.2 (predicción por homología) y §5.5.1 (preprocesamiento RNA-seq) del TFG. Los archivos `.ga` son exportaciones nativas; cualquier instancia de Galaxy puede importarlos.

## 1. `Aquaporin_tblastn-exonerate-filter.ga`

Workflow del apartado **Predicción basada en homología** del TFG.

- **Inputs** (todos depositados junto al `.ga`, salvo el genoma y el GFF3 por tamaño):
  - **Genoma de referencia FASTA** de *Fragaria* x *ananassa* 'Benihoppe' (versión GDR; no incluido por tamaño, accesible en [GDR](https://www.rosaceae.org/Analysis/18085091)). Si se parte del ensamblaje NCBI (`GCA_034370585.1`) en lugar de la versión GDR, hace falta además `tabla_traduccion_NCBI_GDR.tabular` (28 cromosomas; incluida en este directorio) que reescribe los nombres de cromosoma NCBI → GDR. Si se usa GDR directamente —como en el TFG— esta tabla no es necesaria.
  - **Anotación GFF3** de 'Benihoppe' (versión GDR; no incluida por tamaño, misma fuente que el genoma), usada por GFFCompare para evaluar coincidencias con la anotación oficial.
  - **`dataset_blastn_exonerate.fasta`** — depositado en este directorio. Contiene **576 secuencias** de proteínas de acuaporina usadas como **plantillas para el tblastn**, compuestas por:
    - **419 acuaporinas de Rosaceae** descargadas de **NCBI RefSeq** mediante un script en Python con Biopython (módulo `Bio.Entrez`). Este es el mismo conjunto que se utilizó para fijar empíricamente el rango de longitudes 140–380 aa (rango que abarca el 98,81 % de la familia). El conjunto crudo está depositado también en `datos/curado/TODAS_ROSACEAE_aqp_prot_FINAL_V13.fasta` e incluye, como controles internos, secuencias de *Fragaria vesca* (Fv) y *Fragaria* x *ananassa* (Fa).
    - **157 acuaporinomas curados** procedentes de literatura: **35** *Arabidopsis thaliana* (Johanson et al., 2001), **30** *Oryza sativa* (Sakurai et al., 2005), **41** *Malus domestica* (Liu et al., 2019) y **51** *Hevea brasiliensis* (Zou et al., 2015).
    - **Conteo total: 419 + 35 + 30 + 41 + 51 = 576 secuencias.**
    - *Sobre el conteo de cabeceras del archivo:* el FASTA exportado por Galaxy contiene 627 cabeceras. Algunas entradas de *Fragaria vesca* y *Fragaria* x *ananassa* aparecen con dos cabeceras —una con la accesión NCBI y la etiqueta de curación original (p. ej., `>[MIP][GB]_ADJ67992.1 ... aquaporin [Fragaria x ananassa]`) y otra con la misma accesión junto al nombre funcional del gen (p. ej., `>ADJ67992.1 FaPIP2;1`)—; ambas se conservan en el exportado. Estas 51 entradas redundantes no afectan al tblastn, ya que los hits coincidentes se colapsan downstream por la lógica de mejor coincidencia única.
- **Pasos**: tblastn → BEDtools (genera regiones candidatas) → Exonerate protein2genome → filtrado por longitud (140–380 aa) y solapamiento.
- **Outputs**: candidatas no redundantes con coordenadas, ID exonerate y secuencia peptídica.

**URL pública**: `https://usegalaxy.org/published/workflow?id=f74ee88aff766e20`

## 2. `Featurecounts_StringTie_RNAseq.ga`

Workflow del apartado **§5.5.1 — Obtención y procesamiento de datos transcriptómicos** (alternativa Galaxy a los scripts shell locales).

- **Inputs**: BAMs alineados con HISAT2 + GFF3 de 'Benihoppe'.
- **Pasos**: featureCounts (matriz de cuentas) + StringTie (TPM por gen).
- **Outputs**: matriz de cuentas por gen y muestra + matriz TPM.

**URL pública**: `https://usegalaxy.org/published/workflow?id=6d4ebda2fad43946`

## Cómo importar

En Galaxy: `Workflow → Import → Upload File` → seleccionar el `.ga`. También se pueden ejecutar directamente desde las URLs públicas sin descargar.

## Equivalencia con los scripts locales

| Workflow Galaxy | Scripts locales en este repo |
|---|---|
| `Aquaporin_tblastn-exonerate-filter.ga` | `scripts/5.2_5.3_homologia_curacion/auditoria_gff_vs_secuencia.py` (auditoría posterior) |
| `Featurecounts_StringTie_RNAseq.ga` | `scripts/5.5_rna_seq/5.5.1_obtencion_procesamiento/{04_align.sh, 06_count.sh}` |
