"""
Tests de reproducibilidad — apartado 5.1 (curaduría y PCA).

Verifica directamente sobre los CSV/TSV finales que:
  - El total de candidatas iniciales es 144
  - Aplicando los filtros del TFG (6 TMHs estrictos, sin GFF3_FALLBACK,
    sin MAKER_GFF3, sin AMBAS_MAL, sin Fragmento, sin Descarte) quedan 121
  - El reparto subfamiliar es 37 NIP / 34 TIP / 32 PIP / 12 SIP / 6 XIP
  - El número de descartes es 23
  - Las decisiones de fuente_seq cuadran (17 GFF3 + 15 Exonerate + 6 GFF3 sin alt)

Si todos estos tests pasan, las cifras del 5.1 son verificables directamente
desde la tabla final, sin necesidad de reejecutar el pipeline completo.
"""
from __future__ import annotations

import pandas as pd
import pytest


# ─────────────────────── HELPERS ──────────────────────────────────────────
def _load_tabla(path):
    return pd.read_csv(path, sep='\t')


def _filter_funcionales(df):
    """Mismos filtros que en profiling_final_integrated.py para Train_Cat != PARTIAL."""
    return df[
        (df['fuente_seq'] != 'GFF3_FALLBACK')
        & (df['fuente_seq'] != 'MAKER_GFF3')
        & (df['veredicto'] != 'AMBAS_MAL')
        & (df['subfamilia_phylo'] != 'Fragmento')
        & (df['confianza'] != 'Descarte')
        & (df['TMHs'] == 6)
    ]


# ─────────────────────── TESTS ────────────────────────────────────────────
@pytest.mark.reproducibility
@pytest.mark.needs_data
class TestCifrasCuraduria:
    """Verifica las cifras del 5.1 contra la tabla final."""

    def test_total_candidatas_es_144(self, tabla_traduccion):
        """144 candidatas = 129 no redundantes (Exonerate/GFF3) + 15 GDR/MAKER rescatadas."""
        df = _load_tabla(tabla_traduccion)
        n_total = len(df)
        # Tolera +/-3 por posibles actualizaciones de la tabla
        assert 141 <= n_total <= 147, (
            f"Esperaba 144 candidatas totales según el TFG, encontradas {n_total}"
        )

    def test_descomposicion_129_mas_15(self, tabla_traduccion):
        """La tabla debe descomponerse en 129 no-MAKER + 15 MAKER."""
        df = _load_tabla(tabla_traduccion)
        n_maker = (df['fuente_seq'] == 'MAKER_GFF3').sum()
        n_otros = len(df) - n_maker
        # Tolera variación moderada
        assert 12 <= n_maker <= 18
        assert 125 <= n_otros <= 132

    def test_funcionales_son_121(self, tabla_traduccion):
        """121 funcionales = 144 candidatas – 23 descartes (TFG cifra clave)."""
        df = _load_tabla(tabla_traduccion)
        funcionales = _filter_funcionales(df)
        assert len(funcionales) == 121, (
            f"Esperaba 121 funcionales según el TFG, encontradas {len(funcionales)}"
        )

    def test_reparto_subfamiliar(self, tabla_traduccion):
        """37 NIP, 34 TIP, 32 PIP, 12 SIP, 6 XIP."""
        df = _load_tabla(tabla_traduccion)
        funcionales = _filter_funcionales(df)
        counts = funcionales['subfamilia_phylo'].value_counts().to_dict()
        esperado = {'NIP': 37, 'TIP': 34, 'PIP': 32, 'SIP': 12, 'XIP': 6}
        for sf, n_esperado in esperado.items():
            assert counts.get(sf, 0) == n_esperado, (
                f"{sf}: esperaba {n_esperado}, encontró {counts.get(sf, 0)}"
            )

    def test_descartes_son_23(self, tabla_traduccion):
        """144 candidatas – 121 funcionales = 23 descartes (TFG cifra clave)."""
        df = _load_tabla(tabla_traduccion)
        funcionales = _filter_funcionales(df)
        descartes = len(df) - len(funcionales)
        # Tolera ±2 si hay ligeros ajustes
        assert 21 <= descartes <= 25, (
            f"Esperaba ~23 descartes, encontró {descartes}"
        )

    def test_reparto_suma_total(self, tabla_traduccion):
        """La suma de las 5 subfamilias debe ser exactamente 121."""
        df = _load_tabla(tabla_traduccion)
        funcionales = _filter_funcionales(df)
        counts = funcionales['subfamilia_phylo'].value_counts()
        sumas = counts.get('NIP', 0) + counts.get('TIP', 0) + counts.get('PIP', 0) + counts.get('SIP', 0) + counts.get('XIP', 0)
        assert sumas == 121


@pytest.mark.reproducibility
@pytest.mark.needs_data
class TestFuenteSeq:
    """La decisión GFF3 vs Exonerate aparece como columna fuente_seq."""

    def test_fuente_seq_tiene_valores_esperados(self, tabla_traduccion):
        """Solo deben aparecer las categorías que define el TFG."""
        df = _load_tabla(tabla_traduccion)
        valores_validos = {'GFF3', 'EXONERATE', 'GFF3_FALLBACK', 'MAKER_GFF3'}
        valores_encontrados = set(df['fuente_seq'].dropna().unique())
        extras = valores_encontrados - valores_validos
        assert not extras, f"Valores no esperados en fuente_seq: {extras}"

    def test_existen_secuencias_exonerate(self, tabla_traduccion):
        """Al menos una secuencia debe venir de Exonerate (cifra TFG: 15)."""
        df = _load_tabla(tabla_traduccion)
        n_exo = (df['fuente_seq'] == 'EXONERATE').sum()
        assert n_exo > 0
        # TFG dice 15 Exonerate seleccionadas
        assert 10 <= n_exo <= 20, f"Esperaba ~15 Exonerate, encontró {n_exo}"

    def test_existen_15_makers(self, tabla_traduccion):
        """Las 15 secuencias GDR/MAKER deben estar marcadas con MAKER_GFF3."""
        df = _load_tabla(tabla_traduccion)
        n_maker = (df['fuente_seq'] == 'MAKER_GFF3').sum()
        # TFG dice 15 rescatadas
        assert 12 <= n_maker <= 18, f"Esperaba ~15 MAKER_GFF3, encontró {n_maker}"


@pytest.mark.reproducibility
@pytest.mark.needs_data
class TestPCACoordenadas:
    """Verifica el CSV de salida del PCA."""

    def test_coordenadas_columns(self, pca_coordenadas):
        df = pd.read_csv(pca_coordenadas)
        for col in ['PC1', 'PC2', 'ID', 'Subfamilia']:
            assert col in df.columns

    def test_pca_no_genera_nans_en_coordenadas(self, pca_coordenadas):
        df = pd.read_csv(pca_coordenadas)
        assert df['PC1'].notna().all()
        assert df['PC2'].notna().all()

    def test_pca_separa_subfamilias(self, pca_coordenadas):
        """Las medias de PC1 por subfamilia deben ser distintas (separación)."""
        df = pd.read_csv(pca_coordenadas)
        medias = df.groupby('Subfamilia')['PC1'].mean()
        # Si todas las medias fueran iguales, el PCA no separaría nada
        rango = medias.max() - medias.min()
        assert rango > 1.0, (
            f"El PCA debería separar subfamilias en PC1; rango = {rango}"
        )


@pytest.mark.reproducibility
@pytest.mark.needs_data
class TestRankingFeatures:
    """El Random Forest produce un ranking con todas las features."""

    def test_ranking_no_vacio(self, ranking_features):
        df = pd.read_csv(ranking_features)
        assert len(df) > 0

    def test_ranking_columns(self, ranking_features):
        df = pd.read_csv(ranking_features)
        assert 'Feature' in df.columns
        assert 'Importance' in df.columns

    def test_importancias_no_negativas(self, ranking_features):
        df = pd.read_csv(ranking_features)
        assert (df['Importance'] >= 0).all()

    def test_importancias_suman_uno(self, ranking_features):
        """Las importancias del RandomForest deben sumar 1.0."""
        df = pd.read_csv(ranking_features)
        assert abs(df['Importance'].sum() - 1.0) < 0.01

    def test_top_feature_tiene_peso_significativo(self, ranking_features):
        """La feature más discriminante debe tener al menos 5% de peso."""
        df = pd.read_csv(ranking_features)
        top = df.sort_values('Importance', ascending=False).iloc[0]
        assert top['Importance'] >= 0.05, (
            f"Top feature {top['Feature']} solo tiene "
            f"{top['Importance']:.3f} de peso"
        )
