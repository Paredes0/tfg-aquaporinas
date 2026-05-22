# Anexo J — Análisis de motivos peptídicos (MEME) de las 121 acuaporinas de *Fragaria* x *ananassa* 'Benihoppe'

> Material de apoyo del apartado 6.2.3. Sustenta la **Figura 5** (frecuencia de los
> 15 motivos por sub-subfamilia) y aporta el detalle cuantitativo de la pérdida del
> motivo NPA del bucle B en la subfamilia SIP.

## 0. Contenido del anexo

| Archivo | Contenido |
|---|---|
| `Anexo_J_MEME_ALL_AQP_121.txt` | Salida completa de MEME sobre las 121 acuaporinas curadas (15 motivos M1–M15, matrices de probabilidad y diagramas por secuencia). Es la fuente de datos de la Figura 5. |
| `Anexo_J_figura_frecuencia_motivos.png` | Mapa de calor de frecuencia de los 15 motivos por sub-subfamilia (versión de la Figura 5). |
| `Anexo_J_tabla_motivos_NPA_SIP.csv` | Variante del motivo NPA (bucles B y E) en las 12 SIP. |
| `Anexo_J_PSSM_motivo_M2.csv` | Matriz de probabilidad (PSSM) del motivo M2 (NPA del bucle B). |
| `Anexo_J_memoria_tecnica.md` | Este documento. |

La frecuencia de motivos por sub-subfamilia (Figura 5) se reproduce ejecutando
`scripts/5.2_5.3_homologia_curacion/analisis_motivos_final.py`, que lee
`data/curado/ALL_AQP.txt` y `data/curado/tabla_aqp_ordenada.csv` vía
`scripts/common/config.py`.

---

## 1. Detalle: conservación del motivo NPA del bucle B (M2) en la subfamilia SIP

Material de apoyo cuantitativo para la afirmación del apartado 6.2.3: ninguna de las
12 SIP curadas conserva el motivo NPA canónico en el bucle B, sino variantes
degeneradas (NPL, NPT, NPS).

### 1.1. Tabla de motivos NPA por isoforma SIP

Archivo: `Anexo_J_tabla_motivos_NPA_SIP.csv` (12 filas × 7 columnas). Origen del dato:
columnas `motivo_B` y `motivo_E` de `data/curado/tabla_aqp_ordenada.csv`, obtenidas
durante la curaduría estructural del apartado 6.1.

| Subgrupo SIP | Nº isoformas | NPA bucle B | NPA bucle E |
|---|---|---|---|
| FaSIP2;2 | 3 | **NPL** (todas) | NPA (todas) |
| FaSIP1 (NPT) | 2 | **NPT** | NPA (ambas) |
| FaSIP1;3 (NPS) | 7 | **NPS** | NPA (3) / NPI (4) |
| **Total** | **12** | **0 NPA canónico** | 8 NPA + 4 NPI |

**Cifra clave:** ninguna de las 12 SIP presenta la secuencia NPA estricta en el bucle B;
en su lugar, tres variantes degeneradas estructuradas por clado (SIP2;2 con NPL, dos
SIP1 con NPT, siete SIP1;3 con NPS).

### 1.2. Matriz de probabilidades del motivo M2 (PSSM)

Archivo: `Anexo_J_PSSM_motivo_M2.csv` (21 posiciones × 20 aminoácidos). Origen del dato:
sección `letter-probability matrix: alength= 20 w= 21 nsites= 116 E= 1.2e-1345` de
`Anexo_J_MEME_ALL_AQP_121.txt`.

Cada fila corresponde a una posición del motivo M2 (consenso `AGISGGHINPAVTFGLALARH`).
Cada columna corresponde a un aminoácido. Los valores son la probabilidad estimada por
MEME a partir de las 116 secuencias del descubrimiento.

**Posición crítica para la discusión (fila 11, "A" del MNPAV):**

| A | C | D | E | F | G | H | I | K | **L** | M | N | P | Q | R | **S** | T | V | W | Y |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **90,5 %** | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | **0 %** | 0 | 0 | 0 | 0 | 0 | **9,5 %** | 0 | 0 | 0 | 0 |

MEME aprende que esa posición admite A (90,5 %) o S (9,5 %), pero nunca L. Las isoformas
SIP2;2 (con NPL) caen fuera del rango aceptado y por eso no se detectan como instancias
del motivo M2. Las dos isoformas SIP1 con NPT no contribuyen a la matriz, pero pasan el
escaneo posterior gracias a la conservación de las demás 20 posiciones del motivo. Las
siete isoformas SIP1;3 con NPS sí están incluidas en el descubrimiento (contribuyen al
9,5 % de S).

### 1.3. Concordancia con la literatura

La degeneración del motivo NPA del bucle B en SIP es un rasgo distintivo descrito en
otras plantas:

- **Wudick et al. (2009)**, *New Phytologist*: SIP es la única subfamilia *plant*-AQP con el motivo NPA degenerado en su primer bucle (NPT, NPL, NPC u otras variantes en lugar de NPA estricto).
- **Nicolas-Espinosa y Carvajal (2022)**, *The Plant Genome*: BoiSIP1;1 de brócoli presenta NPT/NPA, en línea con el patrón descrito aquí en *F.* x *ananassa*.

El presente trabajo aporta el patrón completo en una especie alo-octoploide: las 12 SIP
muestran las tres variantes (NPL/NPT/NPS) estructuradas por clado homeólogo, sin ningún
caso de NPA canónico en el bucle B.
