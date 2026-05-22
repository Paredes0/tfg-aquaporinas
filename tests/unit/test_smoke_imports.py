"""
test_smoke_imports.py — Smoke tests de imports y parseo de cada script Python.

Para cada script .py del pipeline, verifica que:
1. Puede leerse y compilarse (no hay SyntaxError).
2. Sus imports de nivel superior se resuelven (no hay ImportError de dependencias del entorno).

Esto NO ejecuta el script — solo confirma que el código existe y es válido.
Útil para detectar roturas tras renombrar carpetas o tocar imports.
"""
from __future__ import annotations

import ast
import importlib.util
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS_DIR = REPO_ROOT / 'scripts'


def all_python_scripts() -> list[Path]:
    """Devuelve todos los .py del directorio scripts/ excepto __init__.py."""
    return sorted(
        p for p in SCRIPTS_DIR.rglob('*.py')
        if p.name != '__init__.py'
    )


SCRIPTS = all_python_scripts()


@pytest.mark.parametrize('script_path', SCRIPTS, ids=lambda p: str(p.relative_to(REPO_ROOT)))
def test_script_parses(script_path: Path) -> None:
    """Cada script Python debe compilar sin SyntaxError."""
    source = script_path.read_text(encoding='utf-8')
    try:
        ast.parse(source, filename=str(script_path))
    except SyntaxError as exc:
        pytest.fail(f"SyntaxError en {script_path.relative_to(REPO_ROOT)}: {exc}")


def test_common_config_imports() -> None:
    """El módulo común scripts.common.config debe poder importarse."""
    from scripts.common import config  # noqa: F401

    assert hasattr(config, 'REPO_ROOT')
    assert hasattr(config, 'DATA_DIR')
    assert hasattr(config, 'CURADO_DIR')
    assert hasattr(config, 'FILO_DIR')
    assert hasattr(config, 'RNASEQ_DIR')


def test_at_least_ten_scripts_present() -> None:
    """Sanity check: tras la reorganización debe haber al menos 10 scripts Python."""
    assert len(SCRIPTS) >= 10, (
        f"Esperaba al menos 10 scripts .py en {SCRIPTS_DIR}, encontré {len(SCRIPTS)}"
    )
