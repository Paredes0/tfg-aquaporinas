<div align="center">

# 🌿 Acuaporinoma de *Fragaria* x *ananassa* 'Benihoppe'

### Código, datos y tests del Trabajo de Fin de Grado

[![Tests](https://github.com/Paredes0/tfg-aquaporinas/actions/workflows/tests.yml/badge.svg)](https://github.com/Paredes0/tfg-aquaporinas/actions/workflows/tests.yml)
[![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)](pyproject.toml)
[![Tests](https://img.shields.io/badge/tests-130%20passing-brightgreen)](tests/)
[![Licencia](https://img.shields.io/badge/licencia-MIT-blue.svg)](LICENSE)
[![DOI](https://img.shields.io/badge/DOI-pendiente%20(Zenodo)-lightgrey)](docs/ZENODO.md)
[![Visores](https://img.shields.io/badge/visores-GitHub%20Pages-2ea44f?logo=github)](https://paredes0.github.io/tfg-aquaporinas/)

*Identificación, curado y caracterización transcripcional del acuaporinoma*
*del genoma alo-octoploide de la fresa cultivada*

</div>

---

**Autor:** Noé Paredes Alfonso · Grado en Biotecnología, Universidad Católica de Murcia (UCAM)
**Directores:** Dr. Juan Nicolás Espinosa y Dr. José Ramón Acosta Motos

Este repositorio contiene el **código**, los **datos derivados** y los **tests de reproducibilidad** del TFG. Tras un `git clone` se pueden **regenerar todas las figuras** y **verificar todas las cifras** citadas en la memoria, sin depender de la máquina original.

> Los datos primarios pesados (genoma ~600 MB, 22 muestras de RNA-seq ~40 GB) no se incluyen por su tamaño: se referencian con sus identificadores públicos (ver [`docs/REPRODUCIBILIDAD.md`](docs/REPRODUCIBILIDAD.md)).

## Dos niveles de reproducción

| Nivel | Punto de partida | Qué reproduce | Dónde corre |
|---|---|---|---|
| **1 — desde el repo** | Datos derivados en `data/` | Figuras, tablas y **todas las cifras** del TFG | Cualquier PC, tras `git clone` |
| **2 — desde datos primarios** | Genoma + 22 FASTQ (SRA) | Matrices de conteo/TPM desde lecturas crudas | Servidor Linux (Galaxy / HISAT2) |

El Nivel 1 es totalmente portable. El Nivel 2 (scripts `.sh` del pipeline de RNA-seq) parte de ~40 GB de datos crudos y se ejecuta en servidor; su salida ya viene incluida en `data/rna_seq/`, así que el Nivel 1 no necesita re-ejecutarlo.

## Estructura

```
tfg-aquaporinas/
├── scripts/
│   ├── 5.2_5.3_homologia_curacion/   Predicción por homología + curado (§6.1)
│   ├── 5.5_rna_seq/                  Análisis RNA-seq (§6.3)
│   │   ├── 5.5.1_obtencion_procesamiento/   Pipeline de servidor (HISAT2 + featureCounts)
│   │   ├── 5.5.2_de_abundancia/             DESeq2 + figuras compuestas
│   │   └── 5.5.3_homeologos/                Agrupamiento + dominancia + visor eFP
│   ├── common/config.py             Rutas centralizadas (apuntan a data/)
│   └── regenerar_figuras.py         Runner: regenera todas las figuras del TFG
├── data/                            Datos derivados, autocontenidos
│   ├── curado/                      Inputs §6.1 (FASTAs, DeepTMHMM, MEME, Pepstats, DeepLoc…)
│   ├── filogenia/                   Árbol final (.treefile, .iqtree, .contree) + alineamiento
│   └── rna_seq/                     Matrices basal, DE y homeólogos
├── workflows/galaxy/                Workflows de Galaxy publicados (.ga + URLs)
├── annexes/                         Anexos A–J del TFG (datos, figuras y scripts reproducibles)
├── tests/                           130 tests (tras git clone → todo en verde)
└── docs/                            Documentación + sitio web (GitHub Pages)
```

> La reconstrucción filogenética (§6.2) se realiza con **IQ-TREE** (herramienta externa, comandos en [`docs/REPRODUCIBILIDAD.md`](docs/REPRODUCIBILIDAD.md)); el árbol resultante viaja en `data/filogenia/` y sus soportes nodales se resumen en el **Anexo I**.

Los scripts de Python leen sus entradas de `data/` a través de `scripts/common/config.py`. La raíz de datos se puede redirigir con la variable de entorno `TFG_DATA_ROOT`.

## Visores interactivos (GitHub Pages)

Dos visores HTML autocontenidos acompañan al TFG, publicados como web estática desde `docs/`:

- **Atlas eFP de grupos homeólogos** — expresión por tejido (TPM) sobre la planta de fresa.
- **Explorador del PCA** — espacio multidimensional de las 121 acuaporinas funcionales.

URL prevista: **https://paredes0.github.io/tfg-aquaporinas/** · Publicación y regeneración en [`docs/GITHUB_PAGES.md`](docs/GITHUB_PAGES.md).

## Cifras clave reproducidas por estos scripts

| Cifra | Apartado | Script | Test |
|---|---|---|---|
| 419 acuaporinas Rosaceae filtradas a 140–380 aa | 6.1 | (descarga externa de NCBI RefSeq) | — |
| 3.168 loci → 4.984 hits Exonerate → 129 secuencias no redundantes | 6.1 | `auditoria_gff_vs_secuencia.py` | `test_cifras_curaduria.py` |
| 91 / 38 idénticas vs discrepantes (GFF3 vs Exonerate) | 6.1 | `auditoria_gff_vs_secuencia.py` | `test_cifras_curaduria.py` |
| **121 acuaporinas funcionales** (37 NIP, 34 TIP, 32 PIP, 12 SIP, 6 XIP) | 6.1 | `clasificacion_integrada_aqp.py` + `profiling_final_integrated.py` | `test_cifras_curaduria.py` |
| 23 descartes (21 sin 6 TMH + 2 NIP con deleción + 3 aminoacilasas) | 6.1 | `profiling_final_integrated.py` | `test_cifras_curaduria.py` |
| Filogenia final: 281 sec / 430 sitios / Q.PLANT+R6 / log L = −45.149,26 | 6.2 | IQ-TREE v3 (externo) → `data/filogenia/` | `test_cifras_filogenia.py` |
| Matriz basal de TPM (129 acuaporinas, 22 muestras, 6 tejidos) | 6.3 | `08_basal_expression.R` | `test_cifras_rnaseq.py` |
| 32 grupos homeólogos, 18 cuartetos completos | 6.3 | `11_homeolog_grouping.py` | `test_cifras_homeologos.py` |
| Dominancia de subgenomas A=9 / B=5 / C=5 / **D=12** | 6.3 | `14_homeolog_de_analysis.R` | `test_cifras_homeologos.py` |

## Cómo ejecutar los tests

```bash
pip install pytest pandas numpy scipy biopython scikit-learn
pytest tests/ -q
```

Los **130 tests** pasan en verde directamente tras un `git clone` (verifican las cifras contra los datos de `data/`). No requieren los datos primarios.

## Cómo regenerar las figuras

`config.py` solo define rutas; **no ejecuta nada**. Para regenerar **todas las figuras** del TFG desde `data/`:

```bash
python scripts/regenerar_figuras.py
```

Genera las **Figuras 4–12** del cuerpo, el visor eFP (Figura 13) y las figuras/tablas de los **Anexos G, I y H**, y las reúne numeradas en `results/figuras_TFG/` (`Figura_04_*.png` … `Figura_13_*.html`), que queda como índice limpio. (La Figura 7 —árbol— se hace en iTOL; las Figuras 1–3 proceden de otras publicaciones.)

## Cómo reproducir el pipeline completo

Ver [`docs/REPRODUCIBILIDAD.md`](docs/REPRODUCIBILIDAD.md). Resumen:

1. Descargar los datos primarios (identificadores en `docs/REPRODUCIBILIDAD.md`).
2. (Opcional) Redirigir `TFG_DATA_ROOT` a tu copia de los datos.
3. Ejecutar los scripts en orden numérico dentro de cada apartado.

## Cómo citar

> Paredes Alfonso, N. (2026). *tfg-aquaporinas: código, datos y anexos del acuaporinoma de* Fragaria *x* ananassa *'Benihoppe'* [Conjunto de datos y software]. Zenodo. https://doi.org/10.5281/zenodo.XXXXXXX

Ver [`CITATION.cff`](CITATION.cff) o [`docs/ZENODO.md`](docs/ZENODO.md) para el procedimiento de obtención del DOI.

## Licencia

Código bajo licencia [MIT](LICENSE). Los datos derivados se distribuyen para uso académico con atribución.
