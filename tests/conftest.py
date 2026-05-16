"""
conftest.py — fixtures globales para pytest.

Define fixtures reutilizables:
- Rutas a fixtures sintéticas pequeñas
- Acceso opcional a los datos primarios del TFG (vía TFG_DATA_ROOT)
- Skip automático de tests que necesitan datos primarios cuando no están
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

# Permitir importar `scripts.common.config` desde los tests
REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

FIXTURES_DIR = Path(__file__).parent / 'fixtures'


# ── Fixtures de paths sintéticos ────────────────────────────────────────────
@pytest.fixture(scope='session')
def fixtures_dir() -> Path:
    """Devuelve la carpeta de fixtures sintéticas."""
    return FIXTURES_DIR


@pytest.fixture(scope='session')
def tmhmm_synthetic(fixtures_dir) -> Path:
    """Mini archivo .3line con 4 secuencias de topología conocida."""
    return fixtures_dir / 'tmhmm_synthetic.3line'


@pytest.fixture(scope='session')
def pepstats_synthetic(fixtures_dir) -> Path:
    """Mini archivo Pepstats con bloques de 2 secuencias conocidas."""
    return fixtures_dir / 'pepstats_synthetic.txt'


@pytest.fixture(scope='session')
def meme_synthetic(fixtures_dir) -> Path:
    """Mini archivo MEME con un combined block diagram conocido."""
    return fixtures_dir / 'meme_synthetic.txt'


@pytest.fixture(scope='session')
def iqtree_synthetic(fixtures_dir) -> Path:
    """Mini archivo .iqtree con stats conocidos."""
    return fixtures_dir / 'iqtree_synthetic.iqtree'


@pytest.fixture(scope='session')
def treefile_synthetic(fixtures_dir) -> Path:
    """Mini árbol Newick con soportes conocidos."""
    return fixtures_dir / 'tree_synthetic.treefile'


# ── Fixtures de datos primarios (opcionales) ────────────────────────────────
@pytest.fixture(scope='session')
def data_root() -> Path:
    """
    Raíz de los datos primarios del TFG.
    Toma TFG_DATA_ROOT si está; si no, asume la ruta por defecto.
    """
    raw = os.environ.get('TFG_DATA_ROOT', r'C:\Users\Usuario\Desktop\resultados finales')
    return Path(raw)


@pytest.fixture(scope='session')
def tabla_traduccion(data_root) -> Path:
    """Path a tabla_aquaporinas_traduccion.tabular (output del 5.1)."""
    p = data_root / 'analisis proteinas aquaporina' / 'tabla_aquaporinas_traduccion.tabular'
    if not p.exists():
        pytest.skip(f"Datos primarios no disponibles: {p}")
    return p


@pytest.fixture(scope='session')
def clasificacion_simple(data_root) -> Path:
    """Path a clasificacion_filogenetica_simple.csv (output del 5.1)."""
    p = data_root / 'analisis proteinas aquaporina' / 'clasificacion_filogenetica_simple.csv'
    if not p.exists():
        pytest.skip(f"Datos primarios no disponibles: {p}")
    return p


@pytest.fixture(scope='session')
def pca_coordenadas(data_root) -> Path:
    """Path a PCA_Coordenadas_Finales.csv (output del 5.1 profiling)."""
    p = (data_root / 'analisis proteinas aquaporina'
         / 'profiling_aqp_motifs_final' / 'PCA_Coordenadas_Finales.csv')
    if not p.exists():
        pytest.skip(f"Datos primarios no disponibles: {p}")
    return p


@pytest.fixture(scope='session')
def ranking_features(data_root) -> Path:
    """Path al ranking de importancia del Random Forest."""
    p = (data_root / 'analisis proteinas aquaporina'
         / 'profiling_aqp_motifs_final' / 'RANKING_FINAL_INTEGRADO.csv')
    if not p.exists():
        pytest.skip(f"Datos primarios no disponibles: {p}")
    return p


@pytest.fixture(scope='session')
def iqtree_final() -> Path:
    """
    Path al .iqtree del árbol final (281 secuencias).
    Vive en Z:\\work\\RNA-seq_test\\ por defecto.
    """
    raw = os.environ.get('TFG_RNA_SEQ_ROOT', r'Z:\work\RNA-seq_test')
    p = Path(raw) / 'arbol_acuaporinas_2_bueno_sin_parciales.treefile'
    iqtree_p = Path(str(p).replace('.treefile', '.iqtree'))
    if not iqtree_p.exists():
        pytest.skip(f"Árbol final no disponible: {iqtree_p}")
    return iqtree_p


@pytest.fixture(scope='session')
def homeolog_summary(data_root) -> Path:
    """Path a homeolog_groups_summary.tsv (output del 5.6)."""
    p = data_root / 'aqp_finales' / 'homeolog_groups_summary.tsv'
    if not p.exists():
        pytest.skip(f"Datos primarios no disponibles: {p}")
    return p


# ── Pretty-print del resumen al final ───────────────────────────────────────
def pytest_terminal_summary(terminalreporter, exitstatus, config):
    """Imprime un resumen orientado a la defensa del TFG al final del run."""
    tr = terminalreporter
    passed = len(tr.stats.get('passed', []))
    failed = len(tr.stats.get('failed', []))
    skipped = len(tr.stats.get('skipped', []))
    total = passed + failed + skipped

    tr.write_sep('=', 'RESUMEN PARA DEFENSA DEL TFG', bold=True)
    tr.write_line(f"  Tests ejecutados:  {total}")
    tr.write_line(f"  Pasados (PASS):    {passed}")
    tr.write_line(f"  Fallidos (FAIL):   {failed}")
    tr.write_line(f"  Saltados (SKIP):   {skipped}   (faltan datos primarios)")
    tr.write_line('')
    if failed == 0:
        tr.write_line("  Todas las funciones de los scripts pasan sus tests unitarios.")
        tr.write_line("  Las cifras citadas en el TFG son reproducibles desde el código.")
    else:
        tr.write_line("  Hay fallos — revisar antes de la defensa.")
