# Workflows Galaxy del TFG

Dos workflows usados respectivamente en §5.2 (predicción por homología) y §5.5.1 (preprocesamiento RNA-seq) del TFG. Los archivos `.ga` son exportaciones nativas; cualquier instancia de Galaxy puede importarlos.

## 1. `Aquaporin_tblastn-exonerate-filter.ga`

Workflow del apartado **§5.2 — Predicción basada en homología**.

- **Inputs**: genoma 'Benihoppe' FASTA + las 419 acuaporinas Rosaceae filtradas a 140-380 aa.
- **Pasos**: tblastn → BEDtools (genera regiones candidatas) → Exonerate protein2genome → filtrado por longitud y solapamiento.
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
