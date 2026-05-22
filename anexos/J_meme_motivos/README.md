# Anexo J — Análisis de motivos peptídicos (MEME)

Citado en §6.2.3 del TFG. Sustenta la **Figura 8** (frecuencia de los 15 motivos MEME por sub-subfamilia) y aporta el detalle cuantitativo de la pérdida del motivo NPA del bucle B en la subfamilia SIP.

## Contenido

| Archivo | Contenido |
|---|---|
| `Anexo_J_MEME_ALL_AQP_121.txt` | Salida completa de MEME sobre las 121 acuaporinas curadas (15 motivos M1–M15, matrices de probabilidad y diagramas por secuencia). Fuente de datos de la Figura 8. |
| `Anexo_J_figura_frecuencia_motivos.png` | Mapa de calor de frecuencia de los 15 motivos por sub-subfamilia (versión de la Figura 8). |
| `Anexo_J_tabla_motivos_NPA_SIP.csv` | Variante del motivo NPA (bucles B y E) en las 12 SIP. |
| `Anexo_J_PSSM_motivo_M2.csv` | Matriz de probabilidad (PSSM) del motivo M2 (NPA del bucle B). |
| `Anexo_J_memoria_tecnica.md` | Memoria técnica con el detalle por isoforma SIP. |

La Figura 8 se reproduce ejecutando `scripts/5.2_5.3_homologia_curacion/analisis_motivos_final.py`, que lee `datos/curado/ALL_AQP.txt` y `datos/curado/tabla_aqp_ordenada.csv` vía `scripts/common/config.py`.

## Cifras de control (uniformes con el TFG)

- 121 acuaporinas curadas · 15 motivos (M1–M15).
- Motivo NPA del **bucle B** en SIP: ninguna de las 12 conserva el NPA canónico → variantes **NPL ×3** (SIP2;2), **NPT ×2** (SIP1;1) y **NPS ×7** (SIP1;3). Coincide con el apartado 6.2.3.
