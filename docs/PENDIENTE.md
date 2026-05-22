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

## 2. Placeholders a sustituir tras obtener el DOI

- [ ] **`CITATION.cff`** (línea 32): sustituir `<usuario>` por tu usuario real de GitHub en `repository-code`.
- [ ] **`README.md`**: sustituir el badge `DOI-pendiente` por el DOI real y actualizar la línea "Cómo citar" (quitar "[pendiente — pendiente del release Zenodo]").
- [ ] **`docs/ZENODO.md`**: el `10.5281/zenodo.XXXXXXX` es ilustrativo; no requiere edición, pero conviene anotar ahí el DOI real cuando llegue.
- [ ] **§10 del TFG.docx** (fuera de este repo): insertar el párrafo "Data availability" + DOI + URL. Texto propuesto en `resultados finales/docs/superpowers/specs/2026-05-21-anexos-tfg-design.md`, sección 8.

## 3. Contenido del repo aún incompleto

### Anexo E — Soportes filogenéticos (PENDIENTE de generar)

`annexes/E_soportes_filo/` contiene solo un `README.md` placeholder. Faltan 4 archivos, **bloqueados** por la disponibilidad del `.treefile` final (vive en `Z:\work\RNA-seq_test\arbol_acuaporinas_2_bueno_sin_parciales.treefile`):

- [ ] `Anexo_E_tabla_277_nodos_soportes.csv` — un nodo por fila: SH-aLRT (1000 iter), UFBoot (1000 iter), aBayes + indicadores de triple soporte.
- [ ] `Anexo_E_figura_histograma_soportes.png` — distribución de los 3 estadísticos con líneas de umbrales (UFBoot ≥ 95 %, SH-aLRT ≥ 80 %, aBayes ≥ 0,95).
- [ ] `Anexo_E_tabla_nodos_subfamilias.csv` — subtabla con los nodos que delimitan PIP, TIP, NIP, SIP, XIP.
- [ ] `Anexo_E_script_reproducible.py` — script (ete3/dendropy) que regenera tabla y figura desde el `.treefile`.

**Cómo desbloquearlo**: copiar el `.treefile` final a local (o montar `Z:`), escribir el script de extracción de soportes y ejecutarlo. Es el único hueco de contenido real del repo.

## 4. Verificaciones opcionales recomendadas antes del push

- [ ] **Revisar `.gitignore`**: los anexos usaron `git add -f` para incluir `.fasta/.png/.pdf` (que el `.gitignore` excluye por defecto). Confirmar que el `.gitignore` no vuelve a excluirlos en futuros `git add` y que ningún dato primario pesado (>10 MB) entró por error. Comprobar tamaño total del repo antes del push (`git count-objects -vH`).
- [ ] **Revisar tamaño de `annexes/`**: confirmar que el repo sigue siendo razonable para GitHub (idealmente < 100 MB). Las figuras PNG y el alineamiento MAFFT son lo más pesado.
- [ ] **Comprobar que `consenso_aqp_fixed.gff3`** (en `annexes/I_curado_gff3_vs_exonerate/`) abre correctamente y tiene las 20 entradas Exonerate-curadas esperadas.

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

**No falta ningún script** del pipeline actual. Los 3 scripts de reanotación se descartaron (decisión del TFG v9); `12_substitute_gff3.py` se rescató porque genera el GFF3 corregido (material del Anexo I).

Cobertura de tests: **133/133** (45 reproducibilidad + 67 unit + 21 smoke). Limitaciones de cobertura documentadas honestamente en `tests/REPORTE_DEFENSA.md` (sección "Qué garantizan estos tests y qué no").

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
