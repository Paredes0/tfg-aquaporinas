# Inputs y conjunto de datos del TFG (§5.1)

Este documento lista los datos primarios y derivados utilizados en el pipeline.

## Datos primarios externos (no incluidos por tamaño — accesibles vía URL pública)

| Dato | Tamaño | Fuente |
|---|---|---|
| Genoma 'Benihoppe' FASTA (`.fa.gz`) | ~600 MB | [Genome Database for Rosaceae](https://www.rosaceae.org/species/fragaria/fragaria_x_ananassa) |
| Anotación GFF3 'Benihoppe' (`.gff3.gz`) | ~150 MB | GDR (misma URL) |
| 22 muestras RNA-seq paired-end | ~40 GB | NCBI SRA, BioProject [`PRJNA1010234`](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1010234) |

Lista completa de accesiones SRR en `design/samples.tsv` (depósito Zenodo).

## Datos de referencia ligeros (depósito Zenodo)

| Dato | Tamaño | Fuente original |
|---|---|---|
| 419 acuaporinas Rosaceae UniProt (filtradas 140-380 aa) | ~80 KB | UniProtKB query `family:aquaporin AND organism:Rosaceae` |
| *Arabidopsis thaliana* AQPome (n=35) FASTA | ~10 KB | Johanson et al. 2001 |
| *Hevea brasiliensis* AQP (n=51) FASTA | ~15 KB | Zou et al. 2015 |
| *Malus domestica* MIP (n=41) FASTA | ~12 KB | Liu et al. 2019 |
| *Oryza sativa* AQP (n=30) FASTA | ~9 KB | Sakurai et al. 2005 |

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
