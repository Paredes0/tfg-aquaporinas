"""
config.py — Configuración centralizada de rutas del repositorio.

Por defecto, los datos viven en la carpeta `datos/` del PROPIO repositorio, de modo
que el pipeline es autorreproducible tras un `git clone` (no depende de rutas de la
máquina original). Estructura esperada:

    datos/
    ├── curado/      Inputs de §5.2-5.3 (FASTAs, DeepTMHMM, MEME, Pepstats, DeepLoc,
    │                tablas de veredictos, datos GDR). Nombres originales.
    ├── filogenia/   Árbol final (.treefile, .iqtree, .contree) + alineamiento ClipKIT.
    └── rna_seq/     Matrices derivadas (basal/, de/, homeologos/).

Uso desde un script del repo:

    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from scripts.common import config

    df = pd.read_csv(config.CURADO_DIR / 'tabla_aquaporinas_traduccion.tabular', sep='\\t')

Overrides por variable de entorno (opcionales):
    TFG_DATA_ROOT     → apunta a otra copia con la misma estructura data/ (curado, filogenia, rna_seq).
    TFG_RESULTS_DIR   → carpeta donde los scripts escriben sus outputs (por defecto resultados/).
"""
from __future__ import annotations

import os
from pathlib import Path

# ── Raíz del repositorio (este archivo está en scripts/common/) ──────────────
REPO_ROOT: Path = Path(__file__).resolve().parents[2]

# ── Raíz de datos ─────────────────────────────────────────────────────────────
# Por defecto: data/ del repo. Override global con TFG_DATA_ROOT.
_data_env = os.environ.get("TFG_DATA_ROOT")
DATA_DIR: Path = Path(_data_env) if _data_env else REPO_ROOT / "datos"

# ── Subcarpetas de datos ──────────────────────────────────────────────────────
CURADO_DIR: Path = DATA_DIR / "curado"          # §5.2-5.3 predicción + curación
FILO_DIR: Path = DATA_DIR / "filogenia"         # §5.4 reconstrucción filogenética
RNASEQ_DIR: Path = DATA_DIR / "rna_seq"         # §5.5 RNA-seq
RNASEQ_BASAL_DIR: Path = RNASEQ_DIR / "basal"
RNASEQ_DE_DIR: Path = RNASEQ_DIR / "de"
RNASEQ_HOM_DIR: Path = RNASEQ_DIR / "homeologos"

# Los datos GDR (15 secuencias MAKER) se almacenan junto a los de curado.
GDR_DIR: Path = CURADO_DIR

# ── Carpeta de salida (outputs intermedios / resultados regenerados) ──────────
RESULTS_DIR: Path = Path(os.environ.get("TFG_RESULTS_DIR", str(REPO_ROOT / "resultados")))


def ensure_results() -> Path:
    """Crea la carpeta de resultados si no existe y la devuelve."""
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    return RESULTS_DIR


def verify_inputs(strict: bool = False) -> dict[str, bool]:
    """Comprueba que las carpetas de datos existen. Útil para depuración."""
    checks = {
        "data/curado": CURADO_DIR.is_dir(),
        "data/filogenia": FILO_DIR.is_dir(),
        "data/rna_seq": RNASEQ_DIR.is_dir(),
    }
    if strict:
        for name, ok in checks.items():
            if not ok:
                raise FileNotFoundError(f"No existe la carpeta de datos: {name} (bajo {DATA_DIR})")
    return checks


if __name__ == "__main__":
    print(f"REPO_ROOT  = {REPO_ROOT}")
    print(f"DATA_DIR   = {DATA_DIR}")
    print(f"CURADO_DIR = {CURADO_DIR}")
    print(f"FILO_DIR   = {FILO_DIR}")
    print(f"RNASEQ_DIR = {RNASEQ_DIR}")
    print(f"RESULTS_DIR= {RESULTS_DIR}")
    print()
    for name, ok in verify_inputs().items():
        print(f"  [{'OK' if ok else '  '}] {name}")
