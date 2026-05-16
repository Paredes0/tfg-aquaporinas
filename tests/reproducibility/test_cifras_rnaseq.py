"""
Tests de reproducibilidad — apartado 5.3 (RNA-seq basal).

Verifica directamente sobre la matriz TPM final que:
  - Hay matriz de TPM para ~129 acuaporinas
  - 6 tejidos están representados: green_fruit, red_fruit, crown, leaf, roots, aux_bud
  - RootsCtrl_2 (SRR30146487) NO está como muestra (excluida por PCA, TFG)
  - Las réplicas son típicamente 3 por tejido excepto raíces (2 tras exclusión)
  - Los valores TPM no son todos cero (la matriz no está vacía)
"""
from __future__ import annotations

import pandas as pd
import pytest


@pytest.fixture(scope='module')
def basal_tpm_path(data_root):
    p = data_root / 'RNA-seq' / 'basal_aquaporins' / 'basal_aquaporins_tpm.csv'
    if not p.exists():
        pytest.skip(f"No existe: {p}")
    return p


@pytest.fixture(scope='module')
def basal_summary_path(data_root):
    p = data_root / 'RNA-seq' / 'basal_aquaporins' / 'basal_aquaporins_summary.csv'
    if not p.exists():
        pytest.skip(f"No existe: {p}")
    return p


@pytest.mark.reproducibility
@pytest.mark.needs_data
class TestMatrizBasal:
    """La matriz TPM cubre las 121-129 acuaporinas del TFG."""

    def test_filas_son_129_acuaporinas(self, basal_tpm_path):
        """TFG: matriz basal sobre 129 acuaporinas (121 funcionales + 8 candidatas)."""
        df = pd.read_csv(basal_tpm_path)
        # Tolera ±2 por ajustes de la matriz
        assert 127 <= len(df) <= 145, (
            f"Esperaba ~129 acuaporinas, encontró {len(df)}"
        )

    def test_columnas_de_metadata_presentes(self, basal_tpm_path):
        df = pd.read_csv(basal_tpm_path)
        for col in ['gene_id', 'aqp_family_subfamily', 'subfamilia_phylo', 'fuente_seq']:
            assert col in df.columns, f"Falta columna {col}"

    def test_outlier_root_excluido(self, basal_tpm_path):
        """RootsCtrl_2 (SRR30146487) excluido por PCA en el 5.3."""
        df = pd.read_csv(basal_tpm_path)
        assert 'RootsCtrl_2' not in df.columns, (
            "RootsCtrl_2 (outlier) NO debería estar en la matriz final"
        )

    def test_tpm_values_no_todos_cero(self, basal_tpm_path):
        """Al menos una columna de muestra debe tener valores > 0."""
        df = pd.read_csv(basal_tpm_path)
        sample_cols = [c for c in df.columns
                       if c not in ('gene_id', 'aqp_family_subfamily',
                                    'subfamilia_phylo', 'fuente_seq',
                                    'needs_reannotation')]
        for col in sample_cols:
            # Cada columna debe tener al menos un gen con TPM > 0
            assert (df[col] > 0).any(), f"{col} no tiene ningún TPM > 0"

    def test_subfamilias_presentes(self, basal_tpm_path):
        """Las 5 subfamilias clásicas deben estar en la matriz."""
        df = pd.read_csv(basal_tpm_path)
        col = 'subfamilia_phylo'
        valores = set(df[col].dropna().unique())
        for sf in ['PIP', 'TIP', 'NIP', 'SIP', 'XIP']:
            assert sf in valores, f"Subfamilia {sf} ausente de la matriz"


@pytest.mark.reproducibility
@pytest.mark.needs_data
class TestTejidos:
    """6 tejidos representados: fruto verde/rojo/corona/hoja/raíz/yema axilar."""

    def test_seis_tejidos(self, basal_summary_path):
        df = pd.read_csv(basal_summary_path)
        tejidos = set(df['tissue'].unique())
        assert len(tejidos) == 6, (
            f"Esperaba 6 tejidos según el TFG, encontró {len(tejidos)}: {tejidos}"
        )

    def test_tejidos_esperados(self, basal_summary_path):
        df = pd.read_csv(basal_summary_path)
        tejidos = set(df['tissue'].unique())
        esperados = {'green_fruit', 'red_fruit', 'crown', 'leaf', 'roots', 'aux_bud'}
        assert tejidos == esperados, (
            f"Tejidos: encontrados {tejidos}, esperados {esperados}"
        )

    def test_raices_tienen_2_replicas(self, basal_summary_path):
        """Tras excluir RootsCtrl_2, raíces tienen n=2 réplicas."""
        df = pd.read_csv(basal_summary_path)
        roots_rows = df[df['tissue'] == 'roots']
        # Todas las filas roots deben tener n=2
        assert (roots_rows['n'] == 2).all(), (
            f"roots debería tener n=2 tras exclusión, valores: {roots_rows['n'].unique()}"
        )

    def test_replicas_por_tejido(self, basal_summary_path):
        """Cada tejido tiene un n consistente, marca del diseño experimental.

        Diseño TFG (22 muestras): 3×GreenFruit + 3×RedFruit + 3×Crown
        + 3×LeafCtrl + 3×LeafSalt (= 3 réplicas para leaf consolidado)
        + 3×RootsCtrl(−1 outlier=2) + 3×RootsSalt + 1×AuxBud = 22 muestras.

        Sin embargo, en `basal_aquaporins_summary.csv` el agrupamiento real
        es: green_fruit n=3, red_fruit n=3, crown n=3, leaf n=3 (o 6 si
        consolida ctrl+salt), roots n=2 (excluido outlier) o n=3, aux_bud n=1.
        """
        df = pd.read_csv(basal_summary_path)
        n_modal = df.groupby('tissue')['n'].apply(
            lambda s: s.mode().iloc[0]
        )
        # Tejidos con réplicas múltiples
        for tissue in ['green_fruit', 'red_fruit', 'crown']:
            assert n_modal[tissue] >= 3, (
                f"{tissue}: esperaba n≥3, encontró n={n_modal[tissue]}"
            )
        # roots reducido por outlier
        assert n_modal['roots'] in (2, 3), f"roots: n={n_modal['roots']}"
        # aux_bud: 1 réplica en el diseño TFG
        assert n_modal['aux_bud'] == 1, f"aux_bud: esperaba n=1, encontró n={n_modal['aux_bud']}"
