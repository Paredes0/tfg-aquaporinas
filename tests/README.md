# Tests

> Suite de tests que verifica que los scripts del TFG funcionan correctamente y reproducen las cifras citadas en el documento.

## Cómo ejecutar

```bash
# Desde la raíz del repo
pytest tests/ -v
```

Salida esperada: **112 tests passed**.

## Estructura

```
tests/
├── conftest.py              Fixtures globales + summary final
├── fixtures/                Datos sintéticos pequeños (TMHMM, Pepstats, MEME, IQTree)
├── unit/                    Tests unitarios (no requieren datos primarios)
│   ├── test_parsers.py            26 tests — TMHMM/Pepstats/MEME/IQTree/Newick
│   ├── test_filters.py            24 tests — longitud, Z-score, Train_Cat, cifras
│   └── test_classification.py     17 tests — subfamilia y consenso filogenético
├── reproducibility/         Verifican cifras del TFG (requieren datos primarios)
│   ├── test_cifras_curaduria.py   17 tests — 5.1: 121 funcionales, reparto, PCA, RF
│   ├── test_cifras_filogenia.py   10 tests — 5.2: 281 sec, Q.PLANT+R6, log L
│   ├── test_cifras_homeologos.py   9 tests — 5.6: 32 grupos, dominancia D=12
│   └── test_cifras_rnaseq.py       9 tests — 5.3: 129 AQPs, 6 tejidos, outlier
└── REPORTE_DEFENSA.md       Resumen ejecutivo para la defensa del TFG
```

## Variables de entorno

Los tests de reproducibilidad leen los CSV/TSV finales del TFG. Por defecto se buscan en:

- `C:\Users\Usuario\Desktop\resultados finales` para `tabla_aquaporinas_traduccion.tabular`, etc.
- `Z:\work\RNA-seq_test` para los outputs del pipeline RNA-seq.

Override si tus datos viven en otra parte:

```bash
# PowerShell
$env:TFG_DATA_ROOT = "D:\my\path\resultados finales"
$env:TFG_RNA_SEQ_ROOT = "/mnt/data/RNA-seq"

# Bash
export TFG_DATA_ROOT="/mnt/data/resultados_finales"
export TFG_RNA_SEQ_ROOT="/mnt/data/RNA-seq"
```

Si un archivo no existe, los tests correspondientes se marcan como `SKIP` (no fallan).

## Filtrar por categoría

```bash
# Solo tests unitarios
pytest tests/unit -v

# Solo tests de reproducibilidad
pytest tests/reproducibility -v

# Solo tests marcados
pytest -m unit
pytest -m reproducibility
pytest -m needs_data
```

## Tiempo de ejecución

- Tests unitarios: ~1 segundo
- Tests de reproducibilidad: ~1 segundo
- **Total: ~2 segundos**

## Cómo añadir un nuevo test

Si descubres una cifra del TFG que no está cubierta, añade un test a `tests/reproducibility/`:

```python
@pytest.mark.reproducibility
@pytest.mark.needs_data
def test_mi_cifra_nueva(tabla_traduccion):
    df = pd.read_csv(tabla_traduccion, sep='\t')
    n = (df['alguna_columna'] == 'algun_valor').sum()
    assert n == VALOR_DEL_TFG, (
        f"Esperaba {VALOR_DEL_TFG} según el TFG, encontró {n}"
    )
```

## Diagnóstico de fallos

Si un test falla, primero revisa:

1. **¿Existen los datos primarios?** Comprueba `python -c "from scripts.common.config import paths, verify_paths; print(verify_paths())"`
2. **¿Has actualizado la tabla?** Si reejecutaste `clasificacion_integrada_aqp.py` o `profiling_final_integrated.py`, la tabla puede tener cifras ligeramente distintas. Los tests usan tolerancia ±2-3 para absorber pequeños ajustes; si las cifras cambian sustancialmente, hay un problema real.
3. **Reporta el fallo** con el output completo de pytest (`pytest tests/test_que_falla.py -v --tb=long`).
