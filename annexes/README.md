# Anexos del TFG — mapa y citas

Esta carpeta agrupa los anexos materializables del TFG con numeración correlativa A–J. Cada anexo es auto-contenido: contiene sus datos derivados, figuras, scripts reproducibles y memoria técnica donde aplique.

| Anexo | Carpeta | Contenido | Citado en TFG v10 |
|---|---|---|---|
| **A** | `A_repo_overview/` | Visión global del repo + inventario de inputs (§5.1) | §10 (disponibilidad de datos y código) |
| **B** | `B_descartes/` | Tabla de las 26 secuencias excluidas del análisis | §6.1 |
| **C** | `C_pca_robustez/` | Prueba de robustez del PCA: elipses A vs B **superpuestas** y en **paneles separados** + distancias de Mahalanobis | §6.1 |
| **D** | `D_meme_motivos/` | Análisis MEME de las 121 (datos + frecuencia de motivos = Fig 5) + motivo NPA del bucle B en SIP | §6.2.3 (Fig 5) |
| **E** | `E_soportes_filo/` | Distribución SH-aLRT/UFBoot/aBayes sobre 277 nodos | §6.2.2 |
| **F** | `F_evidencia_5polemicas/` | Evidencia visual del descarte de las 5 polémicas | §6.1 |
| **G** | `G_figuras_TFG/` | Versiones alta resolución de las figuras del cuerpo | §6.1 a §6.3 |
| **H** | `H_curado_gff3_vs_exonerate/` | Inputs y outputs del curado GFF3 vs Exonerate (FASTAs, DeepTMHMM, Pepstats, DeepLoc, MEME, tabla veredictos, `Anexo_H_consenso_aqp_fixed.gff3` para actualización oficial) + figura de identidad (Fig 2) | §6.1 (Fig 2) |
| **I** | `I_rnaseq_datos/` | Matrices RNA-seq (TPM basal, DE, homeólogos) + informe MultiQC | §6.3 (Fig 6–9) |
| **J** | `J_parametros_fisicoquimicos/` | Tabla por gen: pI, Mw, nº TMH, variante NPA-B/E, localización | §6.1 y §6.2.3 |

## Política de archivos

- Convención de nombres: `Anexo_<LETRA>_<descripción>.<ext>`.
- Datos primarios pesados (>100 MB) no se incluyen: ver `A_repo_overview/INPUTS_CONJUNTO_DATOS.md`.

## Historial de cambios

- **2026-05-22**: el antiguo **Anexo F** (paneles del PCA) se fusionó en el **Anexo C**, que ya contenía las dos vistas de la misma prueba (elipses superpuestas + paneles separados). Los anexos posteriores se renumeraron: G→F, H→G, I→H.
- **2026-05-22**: el **Anexo D** se amplió de "MEME SIP" a "MEME de las 121" (carpeta `D_meme_sip/` → `D_meme_motivos/`), incorporando los datos MEME completos y la figura de frecuencia de motivos (Fig 5). Se añadieron dos anexos para cubrir contenido citado en el TFG que faltaba: **Anexo I** (datos RNA-seq + MultiQC) y **Anexo J** (parámetros fisicoquímicos por gen).
