# Pendiente para completar el repositorio `tfg-aquaporinas`

> Estado a 2026-05-22. El sub-proyecto de estructuración de anexos (#37) está cerrado y mergeado a `main` (133/133 tests verdes). Este documento lista lo que **aún falta** para dejar el repositorio 100 % listo para la entrega y para su publicación con DOI.

---

## 1. Acciones humanas de publicación (Noé) — bloqueantes para el DOI

Estas acciones no las puede hacer Claude; requieren tu cuenta y tus credenciales.

- [ ] **Crear el repositorio en GitHub** `tfg-aquaporinas` como **público** (vacío, sin README inicial).
- [ ] **Conectar GitHub ↔ Zenodo** una sola vez en [zenodo.org/account/settings/github/](https://zenodo.org/account/settings/github/) y activar el switch del repo. Detalle paso a paso en `docs/ZENODO.md`.
- [ ] **Primer push**: `git push -u origin main`.
- [ ] **Tag de release**: `git tag -a v1.0.0-tfg-entregado -m "TFG entregado para defensa"` + `git push origin --tags`. Zenodo emitirá el DOI automáticamente.
- [ ] **Subir los datos derivados a Zenodo** (~250 KB). Lista completa en `annexes/A_repo_overview/INPUTS_CONJUNTO_DATOS.md` (sección "Datos de referencia ligeros" + "Datos derivados intermedios").
- [ ] **Activar GitHub Pages** una vez (Settings → Pages → Deploy from a branch → `main` / `/docs`). Publica los visores interactivos (eFP de homeólogos + PCA). Pasos en `docs/GITHUB_PAGES.md`. URL prevista: `https://paredes0.github.io/tfg-aquaporinas/`.

## 2. Placeholders a sustituir tras obtener el DOI

- [ ] **`CITATION.cff`** (línea 32): sustituir `<usuario>` por tu usuario real de GitHub en `repository-code`.
- [ ] **`README.md`**: sustituir el badge `DOI-pendiente` por el DOI real y actualizar la línea "Cómo citar" (quitar "[pendiente — pendiente del release Zenodo]").
- [ ] **`docs/ZENODO.md`**: el `10.5281/zenodo.XXXXXXX` es ilustrativo; no requiere edición, pero conviene anotar ahí el DOI real cuando llegue.
- [ ] **§10 del TFG.docx** (fuera de este repo): insertar el párrafo "Data availability" + DOI + URL. Texto propuesto en `resultados finales/docs/superpowers/specs/2026-05-21-anexos-tfg-design.md`, sección 8.

## 3. Contenido del repo aún incompleto

### Anexo E — Soportes filogenéticos ✅ COMPLETADO (2026-05-22)

Generado a partir del árbol final (282 hojas, 277 nodos internos, Q.PLANT+R6, log L = −45.149,26). `annexes/E_soportes_filo/` ya contiene:

- [x] `Anexo_E_script_reproducible.py` — parsea soportes SH-aLRT/aBayes/UFboot y genera tabla, histograma y subtabla.
- [x] `Anexo_E_tabla_277_nodos_soportes.csv` — 277 nodos con sus 3 soportes + indicadores de umbral.
- [x] `Anexo_E_figura_histograma_soportes.{png,pdf}` — 116 nodos (41,9 %) con triple soporte alto.
- [x] `Anexo_E_tabla_nodos_subfamilias.csv` — PIP/NIP/SIP triple alto, TIP monofilético, XIP basal.

La cifra de 116 nodos (41,9 %) coincide exactamente con la documentada en `README_MAPA_ANEXOS.md`. **No quedan huecos de contenido en el repo.**

## 4. Verificaciones opcionales recomendadas antes del push

- [ ] **Revisar `.gitignore`**: los anexos usaron `git add -f` para incluir `.fasta/.png/.pdf` (que el `.gitignore` excluye por defecto). Confirmar que el `.gitignore` no vuelve a excluirlos en futuros `git add` y que ningún dato primario pesado (>10 MB) entró por error. Comprobar tamaño total del repo antes del push (`git count-objects -vH`).
- [ ] **Revisar tamaño de `annexes/`**: confirmar que el repo sigue siendo razonable para GitHub (idealmente < 100 MB). Las figuras PNG y el alineamiento MAFFT son lo más pesado.
- [ ] **Comprobar que `consenso_aqp_fixed.gff3`** (en `annexes/H_curado_gff3_vs_exonerate/`) abre correctamente y tiene las 20 entradas Exonerate-curadas esperadas.

## 5. Estado de los scripts (completo — sin pendientes)

El pipeline está cubierto y sincronizado con la numeración del TFG v9:

| Carpeta | Cubre §TFG v9 | Scripts |
|---|---|---|
| `scripts/5.2_5.3_homologia_curacion/` | §5.2 Predicción + §5.3 Curación | 7 (incluye `12_substitute_gff3.py` rescatado) |
| `scripts/5.4_filogenia/` | §5.4 Reconstrucción filogenética | 3 |
| `scripts/5.5_rna_seq/5.5.1_obtencion_procesamiento/` | §5.5.1 | 10 (.sh) |
| `scripts/5.5_rna_seq/5.5.2_de_abundancia/` | §5.5.2 | R + eFP + 5 scripts de figuras |
| `scripts/5.5_rna_seq/5.5.3_homeologos/` | §5.5.3 | 4 |
| `scripts/common/` | infraestructura | 2 |

**No falta ningún script** del pipeline actual. Los 3 scripts de reanotación se descartaron (decisión del TFG v9); `12_substitute_gff3.py` se rescató porque genera el GFF3 corregido (material del Anexo H).

Cobertura de tests: **133/133** (45 reproducibilidad + 67 unit + 21 smoke). Limitaciones de cobertura documentadas honestamente en `tests/REPORTE_DEFENSA.md` (sección "Qué garantizan estos tests y qué no").

### 5.1. Reproducibilidad de rutas — refactor a config.py (2026-05-22) ✅

- `scripts/common/config.py` reescrito: ahora apunta a `data/` del repo por defecto; los scripts ya no tienen rutas a la máquina original.
- `data/` añadido (curado/, filogenia/, rna_seq/) con los datos derivados (verificados idénticos a los originales por md5).
- `profiling_final_integrated.py` y `analisis_motivos_unificado.py` leen de `config.CURADO_DIR`. **Verificado**: el PCA reproduce las 121 funcionales leyendo de `data/curado/` (exit 0).
- Los 5 `compose_fig*.py`, los dos visores eFP y el script del Anexo C se refactorizaron a `config.py` (ver §5.2): leen de `data/` y escriben en `results/`.

### 5.2. Verificación end-to-end pendiente (fase "luego") ⏳

Falta **ejecutar cada script y confirmar que reproduce sus cifras**, no solo el PCA:
- [ ] `analisis_motivos_unificado.py` — necesita `clasificacion_filogenetica_simple.csv`, que genera `clasificacion_integrada_aqp.py`; confirmar el orden de ejecución (clasificación → motivos).
- [ ] Scripts de `5.4_filogenia/` (comparar_arboles, rename_tree_nodes, update_prune_ids) — verificar que leen los treefiles de `data/filogenia/`.
- [ ] Scripts R (`07_de_analysis.R`, `08_basal_expression.R`, homeólogos) — requieren entorno R + DESeq2.
- [x] **Refactor de los 5 `compose_fig*.py` a `config.RNASEQ_DIR`** ✅ (2026-05-22): leen de `data/rna_seq/` y escriben en `results/figuras_rnaseq/`; las 5 figuras (10 PDF/PNG) se regeneran desde el repo.
- [x] **Refactor de los dos visores eFP** (`10_generate_efp_viewer.py`, `15_homeolog_efp_viewer.py`) a config.py ✅ (2026-05-22): SVG `fxa_vectorizado.svg` incorporado a `data/rna_seq/`. (El eFP individual queda obsoleto; no se publica.)
- [x] **Refactor del script del Anexo C** a config.py ✅ (2026-05-22): lee `PCA_Coordenadas_Finales.csv` de `results/profiling_aqp_motifs_final/` y la tabla de `data/curado/`; regenera las dos figuras (paneles + elipses superpuestas) idénticas a las commiteadas. El antiguo Anexo F se fusionó en el Anexo C.
- [ ] Scripts R (`07_de_analysis.R`, `08_basal_expression.R`, homeólogos) — requieren entorno R + DESeq2.
- [ ] Auditoría completa reejecutando todo el repo (acordada con Noé).

## 6. Sub-proyecto relacionado (fuera de este repo)

- [ ] **Sub-proyecto #38 — Secuencias parciales con expresión**: verificar si las 8 candidatas parciales identificadas en §6.1 tienen expresión en RNA-seq (ya se conocen 2: `Fxa7Dg01388` PIP TPM máx 43,14 y `Fxa6Dg03790` SIP TPM máx 3,01). Si la tienen, redactar un párrafo en §6.3 del TFG mencionándolas como candidatas a reanotación futura. Es trabajo de redacción del TFG, no del repo, pero su material (`reannotation_candidates.tsv`) está referenciado en `annexes/B_descartes/`.

---

## Resumen de prioridades

| Prioridad | Tarea | Bloqueante |
|---|---|---|
| 🔴 Alta | Crear GitHub + Zenodo + push + tag (sección 1) | Para obtener el DOI |
| 🔴 Alta | Sustituir DOI/usuario placeholders (sección 2) | Tras el DOI |
| 🟡 Media | Generar Anexo E (sección 3) | `.treefile` accesible |
| 🟡 Media | Verificar tamaño repo + `.gitignore` (sección 4) | Antes del push |
| 🟢 Baja | Sub-proyecto #38 (sección 6) | Independiente, redacción TFG |
