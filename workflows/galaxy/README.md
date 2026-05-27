# Workflows Galaxy del TFG

Dos workflows usados respectivamente en §5.2 (predicción por homología) y §5.5.1 (preprocesamiento RNA-seq) del TFG. Los archivos `.ga` son exportaciones nativas; cualquier instancia de Galaxy puede importarlos.

## 1. `Aquaporin_tblastn-exonerate-filter.ga`

Workflow del apartado **Predicción basada en homología** del TFG.

- **Inputs**:
  - Genoma de referencia FASTA de *Fragaria* x *ananassa* 'Benihoppe' (versión GDR). Si se parte del ensamblaje NCBI (`GCA_034370585.1`) en lugar del de GDR, hace falta además el archivo `Galaxy4-[Tabla_Traduccion].tabular` (no incluido en el repo) que reescribe los nombres de cromosoma NCBI → GDR. Si se usa GDR directamente —como en el TFG— esta tabla no es necesaria.
  - Anotación GFF3 de 'Benihoppe' (versión GDR), usada por GFFCompare para evaluar coincidencias con la anotación oficial.
  - **419 secuencias de acuaporinas de Rosaceae** descargadas de UniProt como plantillas BLAST. Sirven también para fijar empíricamente el rango de longitudes (140–380 aa) que se aplica luego como filtro de viabilidad biológica.
  - **Acuaporinomas curados** de las cuatro especies de referencia, también como plantillas BLAST: *Arabidopsis thaliana* (Johanson et al., 2001), *Hevea brasiliensis* (Zou et al., 2015), *Malus domestica* (Liu et al., 2019) y *Oryza sativa* (Sakurai et al., 2005).
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
