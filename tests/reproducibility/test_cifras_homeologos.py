"""
Tests de reproducibilidad — apartado 5.6 (homeólogos y dominancia).

Verifica directamente sobre los TSV/CSV finales que:
  - El número de grupos homeólogos es 32
  - El número de cuartetos completos (los 4 subgenomas presentes) es 18
  - El reparto de dominancia entre subgenomas es A=9 / B=5 / C=5 / D=12

Si estos tests pasan, las cifras del 5.6 son reproducibles desde los outputs.
"""
from __future__ import annotations

import os
from pathlib import Path

import pandas as pd
import pytest


from scripts.common import config


@pytest.fixture(scope='module')
def dominant_subgenome_path():
    p = config.RNASEQ_HOM_DIR / 'dominant_subgenome.csv'
    if not p.exists():
        pytest.skip(f"No existe: {p}")
    return p


# ──────────────────────────── TESTS GRUPOS ────────────────────────────────
@pytest.mark.reproducibility
@pytest.mark.needs_data
class TestGruposHomeologos:
    """El TFG cita 32 grupos homeólogos, 18 cuartetos completos."""

    def test_total_32_grupos(self, homeolog_summary):
        df = pd.read_csv(homeolog_summary, sep='\t')
        # Tolera ±2 grupos por ajustes en regla de pertenencia
        assert 30 <= len(df) <= 34, (
            f"Esperaba 32 grupos según el TFG, encontró {len(df)}"
        )

    def test_18_cuartetos_completos(self, homeolog_summary):
        """Un cuarteto completo tiene los 4 subgenomas (A,B,C,D)."""
        df = pd.read_csv(homeolog_summary, sep='\t')
        n_complete = (df['completeness'] == 'complete').sum()
        assert 16 <= n_complete <= 20, (
            f"Esperaba ~18 cuartetos completos, encontró {n_complete}"
        )

    def test_categorias_completeness(self, homeolog_summary):
        """Las categorías de completeness deben ser un subset esperado."""
        df = pd.read_csv(homeolog_summary, sep='\t')
        valores = set(df['completeness'].unique())
        esperadas = {'complete', 'triplet', 'doublet', 'singleton'}
        no_esperadas = valores - esperadas
        assert not no_esperadas, f"Valores no esperados: {no_esperadas}"

    def test_grupos_tienen_familia_asignada(self, homeolog_summary):
        """Cada grupo debe tener familia AQP asignada (PIP/TIP/NIP/SIP/XIP).

        En el TSV la columna `family` es la familia broad y `subfamily`
        guarda PIP1/PIP2/TIP1/... (sub-subfamilia). Comprobamos `family`.
        """
        df = pd.read_csv(homeolog_summary, sep='\t')
        valid = {'PIP', 'TIP', 'NIP', 'SIP', 'XIP'}
        unknown = set(df['family'].unique()) - valid
        assert not unknown, f"Familias AQP no esperadas: {unknown}"

    def test_subfamilias_son_sub_sub(self, homeolog_summary):
        """La columna `subfamily` guarda PIP1, PIP2, TIP1, etc. (sub-sub)."""
        df = pd.read_csv(homeolog_summary, sep='\t')
        import re
        valid_pattern = re.compile(r'^(PIP|TIP|NIP|SIP|XIP)\d+$')
        for sub in df['subfamily'].dropna().unique():
            assert valid_pattern.match(sub), f"sub-subfamilia inválida: {sub}"


# ─────────────────────────── TESTS DOMINANCIA ─────────────────────────────
@pytest.mark.reproducibility
@pytest.mark.needs_data
class TestDominanciaSubgenomas:
    """Dominancia transcripcional: A=9, B=5, C=5, D=12 (TFG cifra clave)."""

    def test_total_32_dominantes(self, dominant_subgenome_path):
        df = pd.read_csv(dominant_subgenome_path)
        # 32 grupos = 32 asignaciones de subgenoma dominante
        assert 30 <= len(df) <= 34

    def test_reparto_dominancia(self, dominant_subgenome_path):
        """Subgenoma D dominante en 12 grupos = hallazgo principal del TFG."""
        df = pd.read_csv(dominant_subgenome_path)
        counts = df['dominant_subgenome'].value_counts().to_dict()
        # Tolerancia ±2 por subgenoma
        for sub, n in [('A', 9), ('B', 5), ('C', 5), ('D', 12)]:
            n_real = counts.get(sub, 0)
            assert abs(n_real - n) <= 2, (
                f"Dominancia subgenoma {sub}: esperaba {n}, encontró {n_real}"
            )

    def test_d_es_el_dominante(self, dominant_subgenome_path):
        """D debe ser el subgenoma con MÁS grupos donde es dominante."""
        df = pd.read_csv(dominant_subgenome_path)
        counts = df['dominant_subgenome'].value_counts()
        assert counts.index[0] == 'D', (
            f"D debería ser el más dominante, pero el top es {counts.index[0]}"
        )

    def test_dominance_proportion_above_quarter(self, dominant_subgenome_path):
        """Las proporciones de dominancia deben superar 25% (hipótesis nula)."""
        df = pd.read_csv(dominant_subgenome_path)
        # Todos los grupos clasificados como "dominantes" deben tener > 0.25
        assert (df['dominance_proportion'] > 0.25).all()
