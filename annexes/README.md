# Anexos del TFG — mapa y citas

Esta carpeta agrupa los anexos materializables del TFG con numeración correlativa A–H. Cada anexo es auto-contenido: contiene sus datos derivados, figuras, scripts reproducibles y memoria técnica donde aplique.

| Anexo | Carpeta | Contenido | Citado en TFG v9 |
|---|---|---|---|
| **A** | `A_repo_overview/` | Visión global del repo + inventario de inputs (§5.1) | §10 (general) |
| **B** | `B_descartes/` | Tabla de las 26 secuencias excluidas del análisis | §6.1 |
| **C** | `C_pca_robustez/` | Prueba de robustez del PCA (versiones A vs B + Mahalanobis) | §6.1 |
| **D** | `D_meme_sip/` | Conservación del motivo NPA del bucle B en SIPs | §6.2.3 (pie Fig 5) |
| **E** | `E_soportes_filo/` | Distribución SH-aLRT/UFBoot/aBayes sobre 277 nodos | §6.2.2 *(pendiente)* |
| **F** | `F_pca_paneles/` | Comparación gráfica del PCA en paneles separados | §6.1 |
| **G** | `G_evidencia_5polemicas/` | Evidencia visual del descarte de las 5 polémicas | §6.1 |
| **H** | `H_figuras_TFG/` | Versiones alta resolución de las figuras del cuerpo | §6.1 a §6.3 |
| **I** | `I_curado_gff3_vs_exonerate/` | Inputs y outputs del curado GFF3 vs Exonerate (FASTAs, DeepTMHMM, Pepstats, DeepLoc, MEME, tabla veredictos, consenso_aqp_fixed.gff3 para actualización oficial) | §6.1 |

## Política de archivos

- Convención de nombres: `Anexo_<LETRA>_<descripción>.<ext>`.
- Datos primarios pesados (>100 MB) no se incluyen: ver `A_repo_overview/INPUTS_CONJUNTO_DATOS.md`.
- Anexo E está pendiente: requiere acceso al `.treefile` del árbol final.
