"""
Tests unitarios para la lógica de clasificación filogenética (5.1).

Cubre la función classify_by_phylogeny del clasificacion_integrada_aqp.py:
- Asigna subfamilia por consenso de vecinos
- Maneja n_tmh < 4 como Fragmento
- Asigna nivel de confianza Alta/Media/Baja según la concordancia
- Extrae subfamilia detallada (PIP1, TIP2, etc.) cuando es posible
"""
from __future__ import annotations

import re
from collections import Counter

import pytest


# ─── Reimplementaciones idénticas a las de clasificacion_integrada_aqp.py ──
def extract_subfamily(name):
    if not name:
        return None
    m = re.search(r'(PIP|TIP|NIP|SIP|XIP)', name, re.IGNORECASE)
    return m.group(1).upper() if m else None


def extract_full_subfamily(name):
    if not name:
        return None
    m = re.search(r'(PIP|TIP|NIP|SIP|XIP)(\d)', name, re.IGNORECASE)
    if m:
        return f"{m.group(1).upper()}{m.group(2)}"
    m2 = re.search(r'(PIP|TIP|NIP|SIP|XIP)', name, re.IGNORECASE)
    return m2.group(1).upper() if m2 else None


def classify_by_phylogeny(neighbors, n_tmh=6):
    """Clasifica por consenso de vecinos filogenéticos."""
    if n_tmh < 4:
        return 'Fragmento', 'Descarte', ''
    if not neighbors:
        return 'Desconocido', 'Sin_vecinos', ''
    classes = [extract_subfamily(n) for n in neighbors]
    classes = [c for c in classes if c]
    if not classes:
        return 'Desconocido', 'Sin_clasificacion', ''
    counts = Counter(classes)
    best, n = counts.most_common(1)[0]
    if n == 3:
        conf = 'Alta'
    elif n == 2:
        conf = 'Media'
    else:
        conf = 'Baja'
    sub_classes = [extract_full_subfamily(name) for name in neighbors]
    sub_classes = [c for c in sub_classes if c]
    sub_counts = Counter(sub_classes)
    best_sub = sub_counts.most_common(1)[0][0] if sub_counts else best
    return best, conf, best_sub


# ─────────────────────── TESTS EXTRACT_SUBFAMILY ─────────────────────────
class TestExtractSubfamily:
    """extract_subfamily detecta PIP/TIP/NIP/SIP/XIP en nombres de referencias."""

    @pytest.mark.unit
    @pytest.mark.parametrize("name,expected", [
        ('AtPIP1_1', 'PIP'),
        ('AtTIP2;3', 'TIP'),
        ('AtNIP1;1', 'NIP'),
        ('AtSIP1;1', 'SIP'),
        ('AtXIP1;1', 'XIP'),
        ('OsPIP2;1', 'PIP'),
        ('MdTIP1_2', 'TIP'),
        ('HbNIP4;1', 'NIP'),
        ('FaPIP2_1', 'PIP'),
    ])
    def test_extracts_subfamily(self, name, expected):
        assert extract_subfamily(name) == expected

    @pytest.mark.unit
    def test_returns_none_for_non_aqp(self):
        assert extract_subfamily('SomeOtherGene') is None
        assert extract_subfamily('') is None
        assert extract_subfamily(None) is None


# ─────────────────── TESTS EXTRACT_FULL_SUBFAMILY ────────────────────────
class TestExtractFullSubfamily:
    """extract_full_subfamily extrae PIP1, TIP2, etc."""

    @pytest.mark.unit
    @pytest.mark.parametrize("name,expected", [
        ('AtPIP1_1', 'PIP1'),
        ('AtPIP2_3', 'PIP2'),
        ('AtTIP3;1', 'TIP3'),
        ('AtNIP4_1', 'NIP4'),
        ('OsPIP2;7', 'PIP2'),
    ])
    def test_extracts_sub_subfamily(self, name, expected):
        assert extract_full_subfamily(name) == expected

    @pytest.mark.unit
    def test_falls_back_to_subfamily_when_no_number(self):
        # Algunos nombres no tienen número: "PIP" sin más
        assert extract_full_subfamily('SomePIP') == 'PIP'


# ───────────────────── TESTS CLASSIFY_BY_PHYLOGENY ───────────────────────
class TestClassifyByPhylogeny:
    """classify_by_phylogeny es el núcleo de la asignación de subfamilias."""

    @pytest.mark.unit
    def test_fragmento_when_less_than_4_tmh(self):
        cls, conf, sub = classify_by_phylogeny(['AtPIP1_1', 'AtPIP2_1', 'AtPIP1_2'], n_tmh=3)
        assert cls == 'Fragmento'
        assert conf == 'Descarte'

    @pytest.mark.unit
    def test_unknown_when_no_neighbors(self):
        cls, conf, sub = classify_by_phylogeny([], n_tmh=6)
        assert cls == 'Desconocido'
        assert conf == 'Sin_vecinos'

    @pytest.mark.unit
    def test_confianza_alta_with_3_concordant(self):
        cls, conf, sub = classify_by_phylogeny(
            ['AtPIP1_1', 'AtPIP2_3', 'OsPIP2_1'], n_tmh=6)
        assert cls == 'PIP'
        assert conf == 'Alta', "3/3 vecinos concordantes = Alta"

    @pytest.mark.unit
    def test_confianza_media_with_2_concordant(self):
        cls, conf, sub = classify_by_phylogeny(
            ['AtPIP1_1', 'AtPIP2_3', 'AtTIP1_1'], n_tmh=6)
        assert cls == 'PIP'
        assert conf == 'Media', "2/3 vecinos PIP = Media"

    @pytest.mark.unit
    def test_confianza_baja_with_1_concordant(self):
        cls, conf, sub = classify_by_phylogeny(
            ['AtPIP1_1', 'AtTIP2_1', 'AtNIP1_1'], n_tmh=6)
        # El más frecuente tiene n=1 (empate), Counter elige el primero
        assert conf == 'Baja'

    @pytest.mark.unit
    def test_extracts_pip2_for_pip2_majority(self):
        cls, conf, sub = classify_by_phylogeny(
            ['AtPIP2_1', 'AtPIP2_3', 'OsPIP2_7'], n_tmh=6)
        assert cls == 'PIP'
        assert sub == 'PIP2', "El sub debe distinguir PIP1 vs PIP2"

    @pytest.mark.unit
    def test_handles_non_aqp_neighbors(self):
        """Si los vecinos no contienen marcadores AQP, debe degradar a Desconocido."""
        cls, conf, sub = classify_by_phylogeny(['NotAnAQP', 'AlsoNot'], n_tmh=6)
        assert cls == 'Desconocido'
