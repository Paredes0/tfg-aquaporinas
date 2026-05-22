# Acuaporinoma de *Fragaria* x *ananassa* 'Benihoppe' — código y tests del TFG

[![Tests](https://img.shields.io/badge/tests-133%20passed-brightgreen)](tests/)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![DOI](https://img.shields.io/badge/DOI-pendiente-lightgrey)](docs/ZENODO.md)

Repositorio de los scripts utilizados en el Trabajo de Fin de Grado **"Identificación, curado y caracterización transcripcional del acuaporinoma de *Fragaria* x *ananassa* 'Benihoppe'"** (Noé Paredes Alfonso, UCAM, Grado en Biotecnología).

Directores: Dr. José Ramón Acosta Motos y Dr. Juan Nicolás Espinosa.

> Este repositorio contiene únicamente el **código** y los **tests** que demuestran que el pipeline reproduce las cifras citadas en el TFG. Los datos primarios (genoma, FASTQ, BAM, alineamientos masivos) se distribuyen a parte por su tamaño — ver `docs/REPRODUCIBILIDAD.md`.

## Estructura

```
tfg-aquaporinas/
├── scripts/
│   ├── 5.2_5.3_homologia_curacion/   Predicción por homología + Curación (§5.2 + §5.3 TFG v9)
│   ├── 5.4_filogenia/                Reconstrucción filogenética (§5.4 TFG v9)
│   ├── 5.5_rna_seq/                  Análisis RNA-seq (§5.5 TFG v9)
│   │   ├── 5.5.1_obtencion_procesamiento/   Pipeline HISAT2 + featureCounts
│   │   ├── 5.5.2_de_abundancia/             DESeq2 + figuras compuestas
│   │   └── 5.5.3_homeologos/                Agrupamiento + dominancia
│   └── common/                       config.py — rutas centralizadas (apunta a data/)
├── data/                             Datos derivados para reproducir el análisis
│   ├── curado/                       Inputs §5.2-5.3 (FASTAs, DeepTMHMM, MEME, Pepstats, DeepLoc, GDR)
│   ├── filogenia/                    Árbol final (.treefile, .iqtree) + alineamiento ClipKIT
│   └── rna_seq/                      Matrices basal, DE y homeólogos
├── workflows/galaxy/                 Workflows Galaxy publicados (.ga + URLs)
├── annexes/                          Anexos del TFG (A–I)
├── tests/                            133 tests automáticos
└── docs/                             Documentación
```

Los scripts de Python leen sus inputs de `data/` a través de `scripts/common/config.py`, por lo que el análisis de curado y filogenia es reproducible tras un `git clone` (ver `docs/REPRODUCIBILIDAD.md`). Los datos primarios pesados (genoma, FASTQ) no se incluyen: se referencian con sus identificadores públicos.

## Visores interactivos (GitHub Pages)

Dos visores HTML autocontenidos acompañan al TFG y se publican como web estática
desde la carpeta `docs/` (una vez activado GitHub Pages):

- **Atlas eFP de grupos homeólogos** — expresión por tejido (TPM) sobre la planta de fresa.
- **Explorador del PCA** — espacio multidimensional de las 121 acuaporinas funcionales.

Comparten estilo (marfil académico, tipografía serif, paleta de subfamilias) y se
regeneran desde los scripts del repo. Instrucciones de publicación y regeneración
en `docs/GITHUB_PAGES.md`. URL prevista: `https://paredes0.github.io/tfg-aquaporinas/`.

## Cifras clave reproducidas por estos scripts

| Cifra | Apartado | Script | Test |
|---|---|---|---|
| 419 aquaporinas Rosaceae filtradas a 140–380 aa | 5.2 | (descarga inteligente externa) | — |
| 3.168 loci → 4.984 hits Exonerate → 129 secuencias no redundantes | 5.2 | `auditoria_gff_vs_secuencia.py` | `test_cifras_curaduria.py` |
| 91 / 38 idénticas vs discrepantes | 5.3 | `auditoria_gff_vs_secuencia.py` | `test_cifras_curaduria.py` |
| **121 acuaporinas funcionales** (37 NIP, 34 TIP, 32 PIP, 12 SIP, 6 XIP) | 5.3 | `profiling_final_integrated.py` + `clasificacion_integrada_aqp.py` | `test_cifras_curaduria.py` |
| 23 descartes (21 sin 6 TMH + 2 NIPs deleción + 3 aminoacilasas) | 5.3 | `profiling_final_integrated.py` | `test_cifras_curaduria.py` |
| Filogenia final: 281 sec / 430 sites / Q.PLANT+R6 / log L = −45.149,26 | 5.4 | `comparar_arboles.py` | `test_cifras_filogenia.py` |
| Matriz basal TPM sobre 129 acuaporinas, 22 muestras, 6 tejidos | 5.5.1 | `08_basal_expression.R` | `test_cifras_rnaseq.py` |
| 32 grupos homeólogos, 18 cuartetos completos | 5.5.3 | `11_homeolog_grouping.py` | `test_cifras_homeologos.py` |
| Dominancia subgenomas A=9 / B=5 / C=5 / **D=12** | 5.5.3 | `14_homeolog_de_analysis.R` | `test_cifras_homeologos.py` |

## Cómo ejecutar los tests

```bash
# Desde la raíz del repo
pip install pytest pandas biopython numpy scipy scikit-learn
pytest tests/ -v
```

Si todos los tests pasan, las cifras del TFG son reproducibles desde el código de este repo.

## Cómo reproducir el pipeline completo

Ver `docs/REPRODUCIBILIDAD.md` para instrucciones detalladas. Resumen:

1. Descargar los datos primarios (lista en `docs/REPRODUCIBILIDAD.md#datos-primarios`).
2. Configurar `scripts/common/config.py` con las rutas locales.
3. Ejecutar los scripts en orden numérico dentro de cada apartado.

## Cómo citar

Si reutilizas este código o los datos derivados, por favor cita:

> Paredes Alfonso, N. (2026). *Identificación, curado y caracterización transcripcional del acuaporinoma de* Fragaria *x* ananassa *'Benihoppe'*. Trabajo de Fin de Grado, Universidad Católica de Murcia (UCAM). DOI: [pendiente — pendiente del release Zenodo].

Ver `CITATION.cff` para el formato citation file estándar o `docs/ZENODO.md` para el procedimiento de obtención del DOI.

## Auditoría

Los scripts han sido auditados estáticamente para detectar bugs e inconsistencias que pudieran invalidar los resultados. Ver `docs/AUDITORIA_SCRIPTS.md` para el informe completo. Resumen:

- ✅ Lógica de filtros, parsers y clasificación: correcta
- ✅ Cifras citadas en el TFG: reproducibles (verificado por tests)
- ✅ Rutas: los scripts de Python leen de `data/` vía `scripts/common/config.py` (variable `TFG_DATA_ROOT` para apuntar a otra copia); los scripts de figuras de RNA-seq usan `TFG_RNA_SEQ_ROOT`
- ⚠️ Ningún bug detectado invalida los resultados del TFG
