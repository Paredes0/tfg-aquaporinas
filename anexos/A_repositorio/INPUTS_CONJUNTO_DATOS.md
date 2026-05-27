# Inputs y conjunto de datos del TFG (§5.1)

Este documento lista los datos primarios y derivados utilizados en el pipeline.

## Datos primarios externos (no incluidos por tamaño — accesibles vía URL pública)

| Dato | Tamaño | Fuente |
|---|---|---|
| Genoma 'Benihoppe' (ensamblaje, `.fa.gz`) | ~600 MB | NCBI [`GCA_034370585.1`](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_034370585.1/) (ASM3437058v1) |
| Anotación GFF3 'Benihoppe' (`.gff3.gz`, v1.0.a1) | ~150 MB | [GDR, Analysis 18085091](https://www.rosaceae.org/Analysis/18085091) |
| 22 muestras RNA-seq paired-end | ~40 GB | NCBI SRA — 4 BioProjects: `PRJNA838938`, `PRJNA715088`, `PRJNA1144869`, `PRJNA632583` |

Lista completa de accesiones SRR (muestra → run → BioProject) en `scripts/5.5_rna_seq/5.5.1_obtencion_procesamiento/SRA.txt`.

## Datos de referencia ligeros (incluidos en el repositorio → Zenodo)

| Dato | Ubicación en el repo | Tamaño | Fuente original |
|---|---|---|---|
| 419 acuaporinas de Rosaceae (NCBI RefSeq); set usado para fijar el rango 140-380 aa de la Figura 4 del TFG | `datos/curado/TODAS_ROSACEAE_aqp_prot_FINAL_V13.fasta` | ~80 KB | Descarga programática vía Biopython (módulo `Bio.Entrez`) a NCBI RefSeq |
| Plantillas BLAST completas: 576 secuencias = 419 Rosaceae NCBI + 35 *At* + 30 *Os* + 41 *Md* + 51 *Hb* (ver nota al pie sobre las 627 entradas brutas del archivo) | `workflows/galaxy/dataset_blastn_exonerate.fasta` | ~230 KB | Composición local en Galaxy a partir de las 419 Rosaceae + los acuaporinomas curados de literatura |
| Tabla de traducción cromosómica NCBI → GDR (28 cromosomas; necesaria solo si se parte del ensamblaje NCBI) | `workflows/galaxy/tabla_traduccion_NCBI_GDR.tabular` | <1 KB | Construida manualmente a partir de los identificadores `CP139748.1`-`CP139775.1` del ensamblaje NCBI |
| *Arabidopsis thaliana* AQPome (n=35) | integrada en el FASTA de plantillas | ~10 KB | Johanson et al. 2001 |
| *Oryza sativa* AQP (n=30) | integrada en el FASTA de plantillas | ~9 KB | Sakurai et al. 2005 |
| *Malus domestica* MIP (n=41) | integrada en el FASTA de plantillas | ~12 KB | Liu et al. 2019 |
| *Hevea brasiliensis* AQP (n=51) | integrada en el FASTA de plantillas | ~15 KB | Zou et al. 2015 |

> *Sobre el FASTA de plantillas BLAST.* El conjunto operativo son las **576** secuencias (419 Rosaceae + 35 *At* + 30 *Os* + 41 *Md* + 51 *Hb*). El archivo exportado por Galaxy presenta **627 cabeceras** porque algunas entradas de *Fragaria vesca* y *Fragaria* x *ananassa* —incluidas dentro de las 419 Rosaceae— aparecen anotadas dos veces: con la accesión NCBI bajo su etiqueta de curación original y con la misma accesión junto al nombre funcional del gen. Las 51 entradas redundantes no afectan al tblastn (los hits coincidentes se colapsan downstream por la lógica de mejor coincidencia única).

## Datos derivados intermedios (depósito Zenodo)

| Dato | Tamaño | Genera |
|---|---|---|
| `design/samples.tsv` (mapeo SRR → tejido → condición) | <1 KB | Subida manual |
| `tabla_Aquaporinas_traduccion.tabular` | ~50 KB | Workflow Galaxy `Aquaporin_tblastn-exonerate-filter.ga` |
| `consenso_aqp.fasta` (consenso curado 129 secs) | ~20 KB | Idem |
| Alineamiento MAFFT 121 funcionales (post-ClipKIT, 430 sites) | ~50 KB | MAFFT E-INS-i externo + ClipKIT |
| `arbol_acuaporinas_2_bueno_sin_parciales.{treefile,iqtree}` | ~50 KB cada | IQ-TREE v3.0.1 |

## Workflows Galaxy (en este repo, `workflows/galaxy/`)

- `Aquaporin_tblastn-exonerate-filter.ga` — usado en §5.2.
- `Featurecounts_StringTie_RNAseq.ga` — usado en §5.5.1 (alternativa a los scripts shell locales).
