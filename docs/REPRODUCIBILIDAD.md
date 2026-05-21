# Reproducibilidad del TFG — guía paso a paso

> Esta guía explica cómo reproducir los resultados del TFG desde cero, partiendo de los datos primarios públicos.

## Datos primarios necesarios

| Dato | Tamaño | Fuente |
|---|---|---|
| Genoma *F.* × *ananassa* 'Benihoppe' (FASTA + GFF3) | ~600 MB comprimido | [Genome Database for Rosaceae](https://www.rosaceae.org/species/fragaria/fragaria_x_ananassa) |
| 22 muestras RNA-seq paired-end | ~40 GB | NCBI SRA, BioProject `PRJNA1010234` (ver `5.5_rna_seq/5.5.1_obtencion_procesamiento/SRA.txt`) |
| Aquaporinas de referencia Rosaceae | ~50 KB | UniProtKB, query `family:aquaporin AND organism:Rosaceae` |
| 4 datasets curados: At AQPome, Hb (*Hevea*), Md (*Malus*), Os (*Oryza*) | ~30 KB | Johanson 2001, Zou 2015, Liu 2019, Sakurai 2005 |

## Variables de entorno

```bash
# Donde viven los datos del TFG (no-RNA-seq)
export TFG_DATA_ROOT="/mnt/data/resultados_finales"

# Donde vive el pipeline RNA-seq (puede ser otra unidad)
export TFG_RNA_SEQ_ROOT="/mnt/data/RNA-seq_test"

# (Opcional, Linux RNA-seq) Donde van los datos pesados
export TFG_RNASEQ_HDD_DIR="/mnt/hdd/rnaseq_data"
```

## Pasos de reproducción

### 5.2 + 5.3 Predicción por homología + Curación de aquaporinas (~1-2 horas)

```bash
# Predicción por homología (Galaxy o BLAST+ local)
# Resultado: 419 aquaporinas Rosaceae filtradas a 140-380 aa, 3168 loci, 4984 hits, 129 sec.

# Auditoría GFF3 vs Exonerate
cd scripts/5.2_5.3_homologia_curacion
python auditoria_gff_vs_secuencia.py   # → 91 idénticas / 38 discrepancias

# Clasificación integrada (decide GFF3 vs Exonerate por gen)
python clasificacion_integrada_aqp.py  # → tabla_aquaporinas_traduccion.tabular

# Generar consenso FASTA para filogenia
python generar_consenso_fasta.py       # → consenso_aqp.fasta

# Análisis de motivos
python analisis_motivos_unificado.py   # → MEME huella por subfamilia

# PCA integrado + Random Forest
python profiling_final_integrated.py   # → 121 funcionales, PCA Fig. 3
```

### 5.4 Filogenia (~30 min en cluster, ~6 h en local)

```bash
cd scripts/5.4_filogenia

# Construir alineamiento (externo: MAFFT v7 E-INS-i + ClipKIT)
# Construir árbol (externo: IQ-TREE v3.0.1 -m MFP -bb 1000 -bnni -alrt 1000 -abayes)
# → arbol_acuaporinas_2_bueno_sin_parciales.{iqtree,treefile}

# Comparar árboles GFF3 vs Exonerate
python comparar_arboles.py

# Renombrar nodos para iTOL
python rename_tree_nodes.py
python update_prune_ids.py
```

### 5.5.1 Obtención y procesamiento RNA-seq (~6-12 horas, Linux)

```bash
cd scripts/5.5_rna_seq/5.5.1_obtencion_procesamiento

# Activar entorno conda con SRA Toolkit, fastp, HISAT2, samtools, featureCounts
micromamba activate rnaseq_aqp

# Ejecución completa
bash run_pipeline.sh                   # ejecuta 00_setup → 06_count en orden
```

### 5.5.2 Expresión diferencial y abundancia basal (~20 min)

```bash
cd scripts/5.5_rna_seq/5.5.2_de_abundancia

# Análisis R
Rscript 07_de_analysis.R               # DESeq2 DE control vs estrés
Rscript 08_basal_expression.R          # TPM matrix
```

### 5.5.3 Homeólogos y dominancia (~20 min)

```bash
cd scripts/5.5_rna_seq/5.5.3_homeologos

python 11_homeolog_grouping.py         # → 32 grupos homeólogos
Rscript 13_homeolog_expression.R       # → TPM colapsado por grupo
Rscript 14_homeolog_de_analysis.R      # → DE DESeq2 colapsado
python 15_homeolog_efp_viewer.py       # → eFP viewer homeólogos
```

## Verificación

Tras completar la reproducción, ejecuta los tests:

```bash
pytest tests/ -v
```

Si las 112 pruebas pasan, has reproducido las cifras del TFG.

## Variantes

- **Sin reejecución completa**: si solo quieres verificar las cifras, descarga los outputs CSV/TSV finales (~5 MB) y ejecuta `pytest tests/reproducibility/`.
- **Sin RNA-seq**: si solo te interesan los apartados 5.2–5.4, salta los pasos 5.5.x. Los tests del 5.5/5.5.3 se marcarán como `SKIP`.

## Soporte

Si encuentras problemas reproduciendo el pipeline, abre una issue en este repositorio o contacta con `noeparedesalf@gmail.com`.
