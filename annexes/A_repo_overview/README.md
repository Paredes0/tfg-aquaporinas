# Anexo A — Visión global del repositorio

Este repositorio contiene el código, los workflows y los anexos del TFG *"Identificación, curado y caracterización transcripcional del acuaporinoma de* Fragaria *x* ananassa *'Benihoppe'"* (Noé Paredes Alfonso, UCAM, 2026).

## Estructura

- `scripts/`: pipeline reproducible (5.2 a 5.5 del TFG).
- `workflows/galaxy/`: workflows Galaxy publicados.
- `annexes/`: 9 anexos auto-contenidos (B-I).
- `tests/`: 122 tests automáticos que verifican las cifras del TFG.
- `docs/`: reproducibilidad, auditoría de scripts, instrucciones Zenodo.

## Cifras verificadas por la suite de tests

- 121 acuaporinas funcionales (37 NIP / 34 TIP / 32 PIP / 12 SIP / 6 XIP).
- Árbol filogenético final: 281 secuencias, 430 sitios, Q.PLANT+R6, log L = −45.149,26.
- 32 grupos homeólogos, 18 cuartetos completos, dominancia subgenómica A=9 / B=5 / C=5 / D=12.

## Cómo ejecutar los tests

```bash
pip install pytest pandas biopython numpy scipy scikit-learn
pytest tests/ -v
```

Tiempo total: ~3 segundos. Resultado esperado: 122/122 PASS.

## Inputs y datos primarios

Ver `INPUTS_CONJUNTO_DATOS.md`.
