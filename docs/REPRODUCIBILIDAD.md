# Reproducibilidad del TFG — guía paso a paso

> Esta guía explica cómo reproducir los resultados del TFG. Hay **dos niveles** de reproducción según de dónde se parta.

## Dos niveles de reproducción

**Nivel 1 — Reproducción del análisis (autocontenida en el repo).**
El curado (§5.2-5.3), la filogenia (§5.4, soportes del Anexo I) y las figuras a partir de las matrices derivadas se reproducen **directamente con el repo**, sin descargas: los datos derivados están en `datos/` y los scripts los leen vía `scripts/common/config.py`. Ejemplo verificado:

```bash
python scripts/5.2_5.3_homologia_curacion/profiling_final_integrated.py
# → reproduce las 121 acuaporinas funcionales (PCA + Random Forest) leyendo de datos/curado/
```

**Nivel 2 — Reproducción desde los datos primarios (la "primera parte").**
Regenerar `datos/` desde cero (genoma + lecturas crudas) requiere los datos primarios pesados, que **no** se incluyen por tamaño. Esta es la parte que un tercero debe montar por su cuenta:

| Punto de entrada | Qué es | Cómo obtenerlo |
|---|---|---|
| **Genoma 'Benihoppe'** | Ensamblaje de referencia alo-octoploide (~600 MB) | NCBI `GCA_034370585.1` (ASM3437058v1); anotación GFF3 oficial en [GDR (Benihoppe v1.0.a1)](https://www.rosaceae.org/Analysis/18085091) |
| **22 muestras RNA-seq** | Lecturas paired-end (~40 GB) | NCBI SRA, 4 BioProjects: `PRJNA838938`, `PRJNA715088`, `PRJNA1144869`, `PRJNA632583` |

### Sobre los workflows Galaxy (`.ga`)

Los dos `.ga` de `workflows/galaxy/` definen la **predicción por homología** (§5.2) y el **preprocesamiento RNA-seq** (§5.5.1). **No se ejecutan por línea de comandos**: son definiciones de workflow de Galaxy y se usan de una de estas formas:

1. **Importarlos en una instancia Galaxy** (web): `Workflow → Import → Upload File` → seleccionar el `.ga`, subir el genoma como input y ejecutar. También se pueden ejecutar desde las URLs públicas (ver `workflows/galaxy/README.md`).
2. **Vía la API de Galaxy** con [BioBlend](https://bioblend.readthedocs.io) (Python), para automatizar la subida del genoma y el lanzamiento.

Para una réplica 100 % en línea de comandos sin Galaxy, habría que reimplementar sus pasos con las herramientas equivalentes (tblastn de BLAST+, BEDtools, Exonerate); el `.ga` documenta los parámetros exactos de cada paso.

## Datos primarios necesarios

| Dato | Tamaño | Fuente |
|---|---|---|
| Genoma *F.* x *ananassa* 'Benihoppe' (ensamblaje) | ~600 MB comprimido | NCBI [`GCA_034370585.1`](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_034370585.1/) (ASM3437058v1) |
| Anotación estructural GFF3 'Benihoppe' (v1.0.a1) | ~150 MB | [GDR, Analysis 18085091](https://www.rosaceae.org/Analysis/18085091) |
| 22 muestras RNA-seq paired-end | ~40 GB | NCBI SRA, 4 BioProjects: `PRJNA838938` / `PRJNA715088` / `PRJNA1144869` / `PRJNA632583` (SRR por muestra en `5.5_rna_seq/5.5.1_obtencion_procesamiento/SRA.txt`) |
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
python analisis_motivos_final.py       # → MEME de las 121 curadas (Figura 8)

# PCA integrado + Random Forest
python profiling_final_integrated.py   # → 121 funcionales, PCA Fig. 3
```

### 5.4 Filogenia (~30 min en cluster, ~6 h en local)

```bash
# Alineamiento (externo: MAFFT v7 E-INS-i + ClipKIT)
# Árbol (externo: IQ-TREE v3.0.1 -m MFP -bb 1000 -bnni -alrt 1000 -abayes)
# → datos/filogenia/arbol_acuaporinas.{iqtree,treefile,contree}
# Los soportes nodales se resumen con:
python anexos/I_soportes_filo/Anexo_I_script_reproducible.py
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

Si las 130 pruebas pasan, has reproducido las cifras del TFG.

## Variantes

- **Sin reejecución completa**: si solo quieres verificar las cifras, descarga los outputs CSV/TSV finales (~5 MB) y ejecuta `pytest tests/reproducibility/`.
- **Sin RNA-seq**: si solo te interesan los apartados 5.2–5.4, salta los pasos 5.5.x. Los tests del 5.5/5.5.3 se marcarán como `SKIP`.

## Soporte

Si encuentras problemas reproduciendo el pipeline, abre una issue en este repositorio o contacta con `noeparedesalf@gmail.com`.
