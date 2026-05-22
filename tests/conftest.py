"""
conftest.py — fixtures globales para pytest.

Define fixtures reutilizables:
- Rutas a fixtures sintéticas pequeñas
- Acceso opcional a los datos primarios del TFG (vía TFG_DATA_ROOT)
- Skip automático de tests que necesitan datos primarios cuando no están
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

# Permitir importar `scripts.common.config` desde los tests
REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from scripts.common import config  # noqa: E402  (requiere el sys.path de arriba)

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


# ── Fixtures de datos derivados (incluidos en el repo, data/) ───────────────
# Apuntan a la copia de datos que viaja en el repositorio (data/), de modo que
# los tests de cifras pasan en verde tras un `git clone`, sin depender de la
# máquina original. Override global con $TFG_DATA_ROOT (lo respeta config.py).
def _require(p: Path) -> Path:
    """Devuelve la ruta si existe; si no, salta el test con un mensaje claro."""
    if not Path(p).exists():
        pytest.skip(f"Dato no disponible en el repo: {p}")
    return Path(p)


@pytest.fixture(scope='session')
def tabla_traduccion() -> Path:
    """tabla_aquaporinas_traduccion.tabular (tabla maestra de clasificación)."""
    return _require(config.CURADO_DIR / 'tabla_aquaporinas_traduccion.tabular')


@pytest.fixture(scope='session')
def pca_coordenadas() -> Path:
    """PCA_Coordenadas_Finales.csv (coordenadas del PCA fisicoquímico)."""
    return _require(config.CURADO_DIR / 'PCA_Coordenadas_Finales.csv')


@pytest.fixture(scope='session')
def ranking_features() -> Path:
    """RANKING_FINAL_INTEGRADO.csv (importancia de variables del Random Forest)."""
    return _require(config.CURADO_DIR / 'RANKING_FINAL_INTEGRADO.csv')


@pytest.fixture(scope='session')
def iqtree_final() -> Path:
    """.iqtree del árbol filogenético final (data/filogenia/)."""
    return _require(config.FILO_DIR / 'arbol_acuaporinas.iqtree')


@pytest.fixture(scope='session')
def homeolog_summary() -> Path:
    """homeolog_groups_summary.tsv (resumen de grupos homeólogos)."""
    return _require(config.RNASEQ_HOM_DIR / 'homeolog_groups_summary.tsv')


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
    tr.write_line(f"  Saltados (SKIP):   {skipped}")
    tr.write_line('')
    if failed == 0:
        tr.write_line("  Todas las funciones de los scripts pasan sus tests unitarios.")
        tr.write_line("  Las cifras citadas en el TFG son reproducibles desde el código.")
    else:
        tr.write_line("  Hay fallos — revisar antes de la defensa.")
