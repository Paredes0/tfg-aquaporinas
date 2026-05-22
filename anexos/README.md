# Anexos del TFG — mapa y citas

Esta carpeta agrupa los anexos materializables del TFG con numeración correlativa A–J, **ordenada por la primera aparición de cada anexo en el TFG v12** (regla APA). Cada anexo es auto-contenido: contiene sus datos derivados, figuras, scripts reproducibles y memoria técnica donde aplique.

| Anexo | Carpeta | Contenido | Citado en TFG v12 |
|---|---|---|---|
| **A** | `A_repositorio/` | Visión global del repo + inventario de inputs + reproducibilidad (incluye workflows de Galaxy) | §5 (disponibilidad de datos), §5.2 (Fig 3) |
| **B** | `B_curado_gff3_vs_exonerate/` | Curado GFF3 vs Exonerate: FASTAs, DeepTMHMM, Pepstats, DeepLoc, MEME, tabla de veredictos, GFF3 corregido + figura de identidad | §5.3 / §6.1 (Fig 5) |
| **C** | `C_datos_rnaseq/` | Matrices RNA-seq (TPM basal, DE, homeólogos) + informe MultiQC | §5.5 / §6.3 (Figs 9–12) |
| **D** | `D_figuras_TFG/` | Versiones alta resolución de las Figuras 4–12 del cuerpo | §6.1 a §6.3 |
| **E** | `E_descartes/` | Tabla de las 26 secuencias excluidas del análisis | §6.1 |
| **F** | `F_evidencia_5polemicas/` | Evidencia visual del descarte de las secuencias polémicas (alineamientos MAFFT) | §6.1 |
| **G** | `G_robustez_pca/` | Robustez del PCA: elipses A vs B (superpuestas + paneles) + distancias de Mahalanobis | §6.1 |
| **H** | `H_parametros_fisicoquimicos/` | Tabla por gen: pI, Mw, nº TMH, isoforma, modelo (GFF3/Exonerate), variantes NPA | §6.1 y §6.2.3 |
| **I** | `I_soportes_filo/` | Distribución SH-aLRT/UFBoot/aBayes sobre los 277 nodos | §6.2.2 |
| **J** | `J_meme_motivos/` | Análisis MEME de las 121 (frecuencia de motivos = Fig 8) + motivo NPA del bucle B en SIP | §6.2.3 (Fig 8) |

## Política de archivos

- Convención de nombres: `Anexo_<LETRA>_<descripción>.<ext>`.
- Datos primarios pesados (>100 MB) no se incluyen: ver `A_repositorio/INPUTS_CONJUNTO_DATOS.md`.

## Historial de cambios

- **2026-05-22**: el antiguo Anexo F (paneles del PCA) se fusionó en el de robustez del PCA; se añadieron los anexos de datos RNA-seq, parámetros fisicoquímicos y se amplió el de MEME.
- **2026-05-22 (renumeración APA)**: las letras se reasignaron por **orden de primera cita en el TFG v12**. Mapa respecto a la numeración temática previa: overview A→A, curado H→B, RNA-seq I→C, figuras G→D, descartes B→E, polémicas F→F, PCA C→G, fisicoquímico J→H, soportes E→I, MEME D→J. La propuesta `resultados finales/ANEXOS_TFG.md` documenta el §10 del TFG con esta numeración.
