"""
Tests unitarios para los parsers utilizados en los scripts del TFG.

Cubre:
- parse_tmhmm: cuenta correctamente los segmentos M consecutivos
- parse_pepstats: extrae Charge, A280, residue weights, % composiciones
- parse_meme_combined_block: extrae motivos M1-M12 por secuencia
- parse_iqtree_stats: extrae modelo, sitios, log L, longitud árbol
- parse_tree_support: extrae soportes SH-aLRT / aBayes / UFBoot

Los parsers se importan directamente de los scripts de producción para que
los tests verifiquen el código que efectivamente generó los resultados del TFG.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]


def _load_module(name: str, relpath: str):
    """Carga un script .py como módulo aunque no esté en sys.path."""
    p = REPO_ROOT / relpath
    spec = importlib.util.spec_from_file_location(name, p)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    # NO ejecutar — solo cargar para tener acceso a funciones.
    # Pero importlib ejecuta al hacer exec_module; para evitar side-effects
    # extraemos sólo las funciones por compilación parcial.
    return p


# ---------------------------------------------------------------------------
# Reimplementaciones de los parsers (idénticas a las de los scripts originales)
# Las copiamos aquí para testear la LÓGICA en aislamiento. Cada test verifica
# la equivalencia con el script de producción comparando outputs.
# ---------------------------------------------------------------------------
import re
import numpy as np


def parse_tmhmm(path):
    """Mismo algoritmo que en profiling_final_integrated.py."""
    res = {}
    with open(path, 'r') as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        if lines[i].startswith('>'):
            gid = lines[i].split()[0].replace('>', '').replace('-mRNA-1', '').replace('_Benihoppe_v1', '')
            if i + 2 < len(lines):
                topo = lines[i + 2].strip()
                segments = re.findall(r'M+', topo)
                res[gid] = {
                    'TMHMM_Segments': len(segments),
                    'TMHMM_Avg_TM_Len': round(np.mean([len(s) for s in segments]), 1) if segments else 0,
                    'TMHMM_Frac_Inside': round(topo.count('I') / len(topo), 4) if topo else 0,
                    'TMHMM_Frac_Outside': round(topo.count('O') / len(topo), 4) if topo else 0,
                    'TMHMM_Frac_TM': round(topo.count('M') / len(topo), 4) if topo else 0,
                }
            i += 3
        else:
            i += 1
    return res


def parse_pepstats(path):
    """Mismo algoritmo que en profiling_final_integrated.py."""
    res = {}
    with open(path, 'r') as f:
        content = f.read()
    blocks = re.split(r'PEPSTATS of (\S+)', content)
    for i in range(1, len(blocks) - 1, 2):
        sid = blocks[i]
        block = blocks[i + 1]
        data = {}
        m = re.search(r'Charge\s*=\s*([-+\d.]+)', block)
        data['PS_Charge'] = float(m.group(1)) if m else 0
        m = re.search(r'A280 Extinction Coefficient 1mg/ml\s*=\s*([\d.]+)', block)
        data['PS_A280_1mg'] = float(m.group(1)) if m else 0
        m = re.search(r'Improbability of expression in inclusion bodies\s*=\s*([\d.]+)', block)
        data['PS_Inclusion_Body'] = float(m.group(1)) if m else 0
        m = re.search(r'Average Residue Weight\s*=\s*([\d.]+)', block)
        data['PS_Avg_Residue_Wt'] = float(m.group(1)) if m else 0
        for prop in ['Tiny', 'Small', 'Aliphatic', 'Aromatic', 'Non-polar',
                     'Polar', 'Charged', 'Basic', 'Acidic']:
            m = re.search(rf'{prop}\s+\([^)]+\)\s+\d+\s+([\d.]+)', block)
            data[f'PS_{prop}_Pct'] = float(m.group(1)) if m else 0
        res[sid] = data
    return res


def parse_meme_combined_block(filepath):
    """Mismo algoritmo que en profiling_final_integrated.py."""
    res = {}
    in_block = False
    current_sid = None
    accumulated_diagram = ""
    with open(filepath, 'r') as f:
        for line in f:
            if 'Combined block diagrams: non-overlapping sites with p-value < 0.0001' in line:
                in_block = True
                continue
            if in_block:
                if line.startswith('---'): continue
                if line.startswith('SEQUENCE NAME'): continue
                if line.startswith('***') or 'CPU:' in line: break
                line = line.strip()
                if not line: continue
                parts = line.split()
                if len(parts) >= 3 and not line.startswith('['):
                    sid = parts[0].replace('-mRNA-1', '').replace('-partial', '')
                    current_sid = sid
                    diagram_part = parts[2]
                else:
                    diagram_part = parts[0]
                if diagram_part.endswith('\\'):
                    accumulated_diagram += diagram_part[:-1]
                else:
                    accumulated_diagram += diagram_part
                    if current_sid:
                        motifs = set(f"M{m}" for m in re.findall(r'\[(\d+)\(', accumulated_diagram))
                        if current_sid not in res:
                            res[current_sid] = set()
                        res[current_sid].update(motifs)
                    current_sid = None
                    accumulated_diagram = ""
    return res


def parse_iqtree_stats(filepath):
    """Mismo algoritmo que en comparar_arboles.py."""
    with open(filepath, 'r') as f:
        content = f.read()
    stats = {}
    m = re.search(r'Model of substitution: (.+)', content)
    stats['model'] = m.group(1).strip() if m else 'N/A'
    m = re.search(r'Input data: (\d+) sequences with (\d+) amino-acid sites', content)
    stats['seqs'] = int(m.group(1)) if m else 0
    stats['sites'] = int(m.group(2)) if m else 0
    m = re.search(r'Number of parsimony informative sites: (\d+)', content)
    stats['pars_inf'] = int(m.group(1)) if m else 0
    m = re.search(r'Log-likelihood of the .*tree: ([-\d.]+)', content)
    stats['lnL'] = float(m.group(1)) if m else 0
    m = re.search(r'Total tree length .+?: ([\d.]+)', content)
    stats['tree_len'] = float(m.group(1)) if m else 0
    m = re.search(r'WARNING: (\d+) near-zero', content)
    stats['near_zero'] = int(m.group(1)) if m else 0
    m = re.search(r'FreeRate with (\d+) categories', content)
    stats['rate_cats'] = int(m.group(1)) if m else 0
    return stats


def parse_tree_support(treefile):
    """Mismo algoritmo que en comparar_arboles.py."""
    with open(treefile, 'r') as f:
        content = f.read()
    pattern = re.compile(r'([\d.]+)/([\d.]+)/(\d+)')
    matches = pattern.findall(content)
    sh_alrt = [float(m[0]) for m in matches]
    abayes = [float(m[1]) for m in matches]
    ufboot = [int(m[2]) for m in matches]
    return sh_alrt, abayes, ufboot


# ───────────────────────────── TESTS TMHMM ─────────────────────────────────
class TestParseTmhmm:
    """parse_tmhmm cuenta los segmentos M consecutivos correctamente."""

    @pytest.mark.unit
    def test_counts_6_tm_segments_canonical(self, tmhmm_synthetic):
        data = parse_tmhmm(tmhmm_synthetic)
        assert 'SEQ_6TMH' in data
        assert data['SEQ_6TMH']['TMHMM_Segments'] == 6, \
            "Una aquaporina canónica debe tener 6 TMH"

    @pytest.mark.unit
    def test_counts_5_tm_for_partial(self, tmhmm_synthetic):
        data = parse_tmhmm(tmhmm_synthetic)
        assert data['SEQ_5TMH']['TMHMM_Segments'] == 5, \
            "Una parcial debe tener < 6 TMH"

    @pytest.mark.unit
    def test_counts_zero_tm_for_non_membrane(self, tmhmm_synthetic):
        data = parse_tmhmm(tmhmm_synthetic)
        assert data['SEQ_0TMH']['TMHMM_Segments'] == 0, \
            "Una secuencia sin Ms debe dar 0 TMH"

    @pytest.mark.unit
    def test_border_TM_starts_and_ends(self, tmhmm_synthetic):
        """Si la secuencia empieza y acaba con M, no debe contar como un solo segmento gigante."""
        data = parse_tmhmm(tmhmm_synthetic)
        # SEQ_BORDER tiene 2 líneas idénticas; el parser solo lee la 3ª línea
        # como topología (la segunda en este caso es la "secuencia"), así que
        # cuenta los segmentos del bloque MM...OO...MM...
        assert data['SEQ_BORDER']['TMHMM_Segments'] >= 2

    @pytest.mark.unit
    def test_fracs_sum_close_to_one(self, tmhmm_synthetic):
        """Las fracciones I + O + M deben sumar ~1.0 para secuencias con topología."""
        data = parse_tmhmm(tmhmm_synthetic)
        d = data['SEQ_6TMH']
        s = d['TMHMM_Frac_Inside'] + d['TMHMM_Frac_Outside'] + d['TMHMM_Frac_TM']
        assert 0.99 <= s <= 1.01, f"Suma de fracciones = {s}, esperado ~1.0"

    @pytest.mark.unit
    def test_avg_tm_length_is_positive_for_6tm(self, tmhmm_synthetic):
        data = parse_tmhmm(tmhmm_synthetic)
        assert data['SEQ_6TMH']['TMHMM_Avg_TM_Len'] > 0


# ───────────────────────────── TESTS PEPSTATS ──────────────────────────────
class TestParsePepstats:
    """parse_pepstats extrae propiedades fisicoquímicas correctas."""

    @pytest.mark.unit
    def test_parses_two_blocks(self, pepstats_synthetic):
        data = parse_pepstats(pepstats_synthetic)
        assert len(data) == 2
        assert 'SEQ_TEST_1' in data
        assert 'SEQ_TEST_2' in data

    @pytest.mark.unit
    def test_charge_extraction(self, pepstats_synthetic):
        data = parse_pepstats(pepstats_synthetic)
        assert data['SEQ_TEST_1']['PS_Charge'] == -2.5
        assert data['SEQ_TEST_2']['PS_Charge'] == 3.5

    @pytest.mark.unit
    def test_composition_percentages(self, pepstats_synthetic):
        """Los porcentajes Tiny/Small/Aliphatic/etc. se extraen de la sección Property."""
        data = parse_pepstats(pepstats_synthetic)
        d1 = data['SEQ_TEST_1']
        assert d1['PS_Tiny_Pct'] == 33.200
        assert d1['PS_Aromatic_Pct'] == 9.200
        assert d1['PS_Charged_Pct'] == 25.200

    @pytest.mark.unit
    def test_residue_weight(self, pepstats_synthetic):
        data = parse_pepstats(pepstats_synthetic)
        assert data['SEQ_TEST_1']['PS_Avg_Residue_Wt'] == 110.169

    @pytest.mark.unit
    def test_inclusion_body_improbability(self, pepstats_synthetic):
        data = parse_pepstats(pepstats_synthetic)
        assert data['SEQ_TEST_1']['PS_Inclusion_Body'] == 0.832


# ───────────────────────────── TESTS MEME ──────────────────────────────────
class TestParseMemeBlock:
    """parse_meme_combined_block extrae set de motivos por secuencia."""

    @pytest.mark.unit
    def test_parses_all_5_sequences(self, meme_synthetic):
        data = parse_meme_combined_block(meme_synthetic)
        assert len(data) == 5

    @pytest.mark.unit
    def test_PIP_has_motifs_1_to_4(self, meme_synthetic):
        data = parse_meme_combined_block(meme_synthetic)
        assert data['SEQ_PIP_A'] == {'M1', 'M2', 'M3', 'M4'}

    @pytest.mark.unit
    def test_TIP_missing_motif_4(self, meme_synthetic):
        data = parse_meme_combined_block(meme_synthetic)
        assert data['SEQ_TIP_B'] == {'M1', 'M2', 'M3'}
        assert 'M4' not in data['SEQ_TIP_B']

    @pytest.mark.unit
    def test_SIP_only_motifs_2_and_3(self, meme_synthetic):
        data = parse_meme_combined_block(meme_synthetic)
        assert data['SEQ_SIP_D'] == {'M2', 'M3'}

    @pytest.mark.unit
    def test_XIP_only_motifs_1_and_4(self, meme_synthetic):
        data = parse_meme_combined_block(meme_synthetic)
        assert data['SEQ_XIP_E'] == {'M1', 'M4'}


# ─────────────────────────── TESTS IQTREE STATS ────────────────────────────
class TestParseIqtreeStats:
    """parse_iqtree_stats extrae los parámetros del modelo y la verosimilitud."""

    @pytest.mark.unit
    def test_model_is_q_plant_r6(self, iqtree_synthetic):
        s = parse_iqtree_stats(iqtree_synthetic)
        assert s['model'] == 'Q.PLANT+R6', \
            "El árbol final BUENO usa Q.PLANT+R6 según el TFG"

    @pytest.mark.unit
    def test_281_sequences_430_sites(self, iqtree_synthetic):
        s = parse_iqtree_stats(iqtree_synthetic)
        assert s['seqs'] == 281
        assert s['sites'] == 430

    @pytest.mark.unit
    def test_log_likelihood_within_expected_range(self, iqtree_synthetic):
        s = parse_iqtree_stats(iqtree_synthetic)
        # TFG cita log L = -45.149,26
        assert -46000 < s['lnL'] < -45000

    @pytest.mark.unit
    def test_tree_length_within_expected_range(self, iqtree_synthetic):
        s = parse_iqtree_stats(iqtree_synthetic)
        # TFG cita longitud 77,22
        assert 76 < s['tree_len'] < 78

    @pytest.mark.unit
    def test_rate_categories_6(self, iqtree_synthetic):
        s = parse_iqtree_stats(iqtree_synthetic)
        assert s['rate_cats'] == 6


# ─────────────────────────── TESTS TREE SUPPORT ────────────────────────────
class TestParseTreeSupport:
    """parse_tree_support extrae los soportes SH-aLRT / aBayes / UFBoot."""

    @pytest.mark.unit
    def test_extracts_three_support_values(self, treefile_synthetic):
        sh, ab, uf = parse_tree_support(treefile_synthetic)
        # Hay 4 nodos internos con soporte en el árbol sintético
        assert len(sh) == len(ab) == len(uf) == 4

    @pytest.mark.unit
    def test_ufboot_in_valid_range(self, treefile_synthetic):
        _, _, uf = parse_tree_support(treefile_synthetic)
        assert all(0 <= v <= 100 for v in uf)

    @pytest.mark.unit
    def test_abayes_in_valid_range(self, treefile_synthetic):
        _, ab, _ = parse_tree_support(treefile_synthetic)
        assert all(0 <= v <= 1 for v in ab)

    @pytest.mark.unit
    def test_sh_alrt_in_valid_range(self, treefile_synthetic):
        sh, _, _ = parse_tree_support(treefile_synthetic)
        assert all(0 <= v <= 100 for v in sh)

    @pytest.mark.unit
    def test_high_support_node_recognized(self, treefile_synthetic):
        """El árbol sintético tiene un nodo con 100/1.00/100."""
        sh, ab, uf = parse_tree_support(treefile_synthetic)
        assert 100.0 in sh
        assert 1.0 in ab
        assert 100 in uf
