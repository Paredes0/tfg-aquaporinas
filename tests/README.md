# Tests

> Suite que verifica que los scripts del TFG funcionan y reproducen las cifras citadas en el documento. Todos los datos que necesitan viajan en `datos/`, así que **pasan en verde tras un `git clone`**, sin descargar nada.

## Cómo ejecutar

```bash
pip install pytest pandas numpy scipy biopython scikit-learn
pytest tests/ -q
```

Salida esperada: **130 passed** (0 fallos, 0 saltados).

## Estructura

```
tests/
├── conftest.py              Fixtures globales (datos de datos/ vía config.py) + resumen final
├── fixtures/                Datos sintéticos pequeños (TMHMM, Pepstats, MEME, IQ-TREE, Newick)
├── unit/                    Tests unitarios sobre funciones aisladas
│   ├── test_parsers.py            26 — parseo de TMHMM/Pepstats/MEME/IQ-TREE/Newick
│   ├── test_filters.py            18 — longitud, Z-score, criterios de descarte, cifras
│   ├── test_classification.py     23 — subfamilia y consenso filogenético
│   └── test_smoke_imports.py      18 — cada script compila y sus imports se resuelven
└── reproducibility/         Verifican las cifras del TFG sobre los datos de datos/
    ├── test_cifras_curaduria.py   17 — §6.1: 121 funcionales, reparto, descartes, PCA, RF
    ├── test_cifras_filogenia.py   10 — §6.2: 281 sec, Q.PLANT+R6, log L, soportes
    ├── test_cifras_homeologos.py   9 — §6.3: 32 grupos, dominancia D=12
    └── test_cifras_rnaseq.py       9 — §6.3: 129 AQPs, 6 tejidos, exclusión del outlier
```

## De dónde leen los datos

Los tests de reproducibilidad leen los CSV/TSV/árboles incluidos en `datos/` del repo, a través de `scripts/common/config.py` (`CURADO_DIR`, `FILO_DIR`, `RNASEQ_DIR`). No hacen falta los datos primarios.

Override opcional (si tienes otra copia de `datos/` con la misma estructura):

```bash
# Bash
export TFG_DATA_ROOT="/ruta/a/datos"
# PowerShell
$env:TFG_DATA_ROOT = "D:\ruta\a\datos"
```

Si algún archivo no existiera, el test correspondiente se marca como `SKIP` (no falla).

## Filtrar por categoría

```bash
pytest tests/unit -q              # solo unitarios
pytest tests/reproducibility -q   # solo reproducibilidad
pytest -m reproducibility         # por marcador
```

## Tiempo de ejecución

~2 segundos en total.

## Cómo añadir un nuevo test

Si descubres una cifra del TFG que no está cubierta, añade un test a `tests/reproducibility/`:

```python
@pytest.mark.reproducibility
def test_mi_cifra_nueva(tabla_traduccion):
    df = pd.read_csv(tabla_traduccion, sep='\t')
    n = (df['alguna_columna'] == 'algun_valor').sum()
    assert n == VALOR_DEL_TFG, f"Esperaba {VALOR_DEL_TFG} según el TFG, encontró {n}"
```

## Diagnóstico de fallos

1. **¿Están los datos?** `python scripts/common/config.py` imprime las rutas y comprueba que `datos/{curado,filogenia,rna_seq}` existen.
2. **¿Reejecutaste un script de curado?** `clasificacion_integrada_aqp.py` o `profiling_final_integrated.py` pueden cambiar ligeramente la tabla; los tests usan tolerancia ±2-3 para absorber ajustes menores. Un cambio grande indica un problema real.
3. **Reporta** con el output completo: `pytest tests/test_que_falla.py -v --tb=long`.
