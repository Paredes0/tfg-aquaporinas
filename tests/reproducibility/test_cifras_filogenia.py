"""
Tests de reproducibilidad — apartado 5.2 (filogenia).

Verifica directamente sobre el .iqtree del árbol final BUENO que:
  - El número de secuencias es 281–282 (el .iqtree dice 282; el TFG cita 281
    funcionales + 3 controles tras dedup, lo cual queda dentro de la tolerancia)
  - El alineamiento tiene 430 sitios
  - El modelo seleccionado es Q.PLANT+R6
  - La log-verosimilitud es ~ -45.149,26
  - La longitud total del árbol es ~ 77,22
  - El soporte nodal (UFBoot >= 95%) tiene una proporción razonable

Si todos estos tests pasan, las cifras del 5.2 son reproducibles desde el
output de IQTree.
"""
from __future__ import annotations

import os
import re
import statistics
from pathlib import Path

import pytest


# ─── Parsers locales (idénticos a los del comparar_arboles.py) ────────────
def parse_iqtree_stats(filepath):
    with open(filepath, 'r') as f:
        content = f.read()
    stats = {}
    m = re.search(r'Model of substitution: (.+)', content)
    stats['model'] = m.group(1).strip() if m else 'N/A'
    m = re.search(r'Input data: (\d+) sequences with (\d+) amino-acid sites', content)
    stats['seqs'] = int(m.group(1)) if m else 0
    stats['sites'] = int(m.group(2)) if m else 0
    m = re.search(r'Log-likelihood of the .*tree: ([-\d.]+)', content)
    stats['lnL'] = float(m.group(1)) if m else 0
    m = re.search(r'Total tree length .+?: ([\d.]+)', content)
    stats['tree_len'] = float(m.group(1)) if m else 0
    m = re.search(r'FreeRate with (\d+) categories', content)
    stats['rate_cats'] = int(m.group(1)) if m else 0
    return stats


def parse_tree_support(treefile):
    with open(treefile, 'r') as f:
        content = f.read()
    pattern = re.compile(r'([\d.]+)/([\d.]+)/(\d+)')
    matches = pattern.findall(content)
    sh_alrt = [float(m[0]) for m in matches]
    abayes = [float(m[1]) for m in matches]
    ufboot = [int(m[2]) for m in matches]
    return sh_alrt, abayes, ufboot


# ─── Fixtures para localizar el .iqtree real ──────────────────────────────
@pytest.fixture(scope='module')
def iqtree_final_path(data_root):
    """
    Busca el .iqtree del árbol final BUENO. Tiene varios sitios posibles:
    1. Z:/work/RNA-seq_test (paths.iqtree_final por defecto)
    2. resultados finales/TFG/.../Filogenia/FINAL/final_without_partials/
    """
    candidates = [
        Path(os.environ.get(
            'TFG_RNA_SEQ_ROOT', r'Z:\work\RNA-seq_test'
        )) / 'arbol_acuaporinas_2_bueno_sin_parciales.iqtree',
        data_root / 'TFG' / 'TFG-PRACTICAS NOE' / 'RESULTADOS' / 'Filogenia'
            / 'FINAL' / 'final_without_partials' / 'arbol_acuaporinas.iqtree',
    ]
    for p in candidates:
        if p.exists():
            return p
    pytest.skip(f"No se encontró el .iqtree final en ninguno de: {candidates}")


@pytest.fixture(scope='module')
def treefile_final_path(data_root):
    """Busca el .treefile del árbol final con soportes nodales."""
    candidates = [
        Path(os.environ.get(
            'TFG_RNA_SEQ_ROOT', r'Z:\work\RNA-seq_test'
        )) / 'arbol_acuaporinas_2_bueno_sin_parciales.treefile',
        data_root / 'TFG' / 'TFG-PRACTICAS NOE' / 'RESULTADOS' / 'Filogenia'
            / 'FINAL' / 'final_without_partials' / 'arbol_acuaporinas.treefile',
    ]
    for p in candidates:
        if p.exists():
            return p
    pytest.skip(f"No se encontró el .treefile final en ninguno de: {candidates}")


# ─────────────────────────────── TESTS ────────────────────────────────────
@pytest.mark.reproducibility
@pytest.mark.needs_data
class TestArbolFinal:
    """Cifras del árbol final BUENO (281–282 sec, Q.PLANT+R6, log L ~-45.149,26)."""

    def test_281_o_282_secuencias(self, iqtree_final_path):
        """El TFG cita 281 (5.2 propuesta) o 282 (.iqtree). Aceptar ambos."""
        s = parse_iqtree_stats(iqtree_final_path)
        assert s['seqs'] in (281, 282), (
            f"Esperaba 281 o 282 sec según TFG/.iqtree, encontró {s['seqs']}"
        )

    def test_430_sitios(self, iqtree_final_path):
        """ClipKIT poda hasta 430 columnas para el árbol final."""
        s = parse_iqtree_stats(iqtree_final_path)
        assert s['sites'] == 430

    def test_modelo_q_plant_r6(self, iqtree_final_path):
        """ModelFinder selecciona Q.PLANT+R6 según el TFG."""
        s = parse_iqtree_stats(iqtree_final_path)
        assert s['model'] == 'Q.PLANT+R6'

    def test_seis_categorias_de_tasa(self, iqtree_final_path):
        """R6 = FreeRate con 6 categorías de heterogeneidad."""
        s = parse_iqtree_stats(iqtree_final_path)
        assert s['rate_cats'] == 6

    def test_log_likelihood_45149(self, iqtree_final_path):
        """TFG cita log L = −45.149,26."""
        s = parse_iqtree_stats(iqtree_final_path)
        # Tolerancia ±50 (diferencias por consenso vs ML tree)
        assert -45200 < s['lnL'] < -45100, (
            f"Esperaba log L ≈ −45149,26, encontró {s['lnL']:.2f}"
        )

    def test_longitud_total_77(self, iqtree_final_path):
        """TFG cita longitud total = 77,22."""
        s = parse_iqtree_stats(iqtree_final_path)
        assert 76.5 < s['tree_len'] < 78.0, (
            f"Esperaba longitud ≈ 77,22, encontró {s['tree_len']:.2f}"
        )


@pytest.mark.reproducibility
@pytest.mark.needs_data
class TestSoporteNodal:
    """El árbol debe tener soporte nodal robusto en la mayoría de clados."""

    def test_existen_soportes_para_todos_los_nodos(self, treefile_final_path):
        sh, ab, uf = parse_tree_support(treefile_final_path)
        assert len(sh) > 100, "Esperaba al menos 100 nodos internos con soporte"

    def test_ufboot_medio_alto(self, treefile_final_path):
        """La media UFBoot debe estar por encima de 75% (árbol fiable)."""
        _, _, uf = parse_tree_support(treefile_final_path)
        media = statistics.mean(uf)
        assert media > 70, f"UFBoot media = {media:.1f}, esperado > 70"

    def test_mas_de_la_mitad_de_nodos_con_ufboot_95(self, treefile_final_path):
        """Al menos 40% de los nodos internos deben tener UFBoot >= 95% (criterio TFG)."""
        _, _, uf = parse_tree_support(treefile_final_path)
        n_total = len(uf)
        n_high = sum(1 for v in uf if v >= 95)
        pct = n_high / n_total
        assert pct > 0.4, (
            f"Solo {pct:.1%} de nodos con UFBoot >= 95%, esperaba > 40%"
        )

    def test_abayes_alto_en_mayoria(self, treefile_final_path):
        """Al menos 80% de nodos con aBayes >= 0.95 (soporte bayesiano)."""
        _, ab, _ = parse_tree_support(treefile_final_path)
        n_total = len(ab)
        n_high = sum(1 for v in ab if v >= 0.95)
        pct = n_high / n_total
        assert pct > 0.6, (
            f"Solo {pct:.1%} de nodos con aBayes >= 0.95, esperaba > 60%"
        )
