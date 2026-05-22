# Anexo W — Conservación del motivo NPA del bucle B (M2) en la subfamilia SIP de *Fragaria* × *ananassa*

> Material de apoyo cuantitativo para la afirmación del apartado 6.2.3: ninguna de las 12 SIPs curadas conserva el motivo NPA canónico en el bucle B, sino variantes degeneradas (NPL, NPT, NPS).

## 1. Tabla de motivos NPA por isoforma SIP

Archivo: `AnexoW_tabla_motivos_NPA_SIP.csv` (12 filas × 7 columnas).

Origen del dato: columnas `motivo_B` y `motivo_E` de `aqp_finales/tabla_aqp_ordenada.csv`, obtenidas durante la curaduría estructural del apartado 6.1.

| Subgrupo SIP | Nº isoformas | NPA bucle B | NPA bucle E |
|---|---|---|---|
| FaSIP2;2 | 3 | **NPL** (todas) | NPA (todas) |
| FaSIP1 (NPT) | 2 | **NPT** | NPA (ambas) |
| FaSIP1;3 (NPS) | 7 | **NPS** | NPA (3) / NPI (4) |
| **Total** | **12** | **0 NPA canónico** | 8 NPA + 4 NPI |

**Cifra clave:** ninguna de las 12 SIPs presenta la secuencia NPA estricta en el bucle B. En su lugar, se observan tres variantes degeneradas estructuradas por clado (SIP2;2 con NPL, dos SIP1 con NPT, siete SIP1;3 con NPS).

## 2. Matriz de probabilidades del motivo M2 (PSSM)

Archivo: `AnexoW_PSSM_motivo_M2_NPA_B.csv` (21 posiciones × 20 aminoácidos).

Origen del dato: sección `letter-probability matrix: alength= 20 w= 21 nsites= 116 E= 1.2e-1345` del archivo `MEME_FINAL/ALL_AQP.txt`, líneas 660-680.

Cada fila corresponde a una posición del motivo M2 (consenso `AGISGGHINPAVTFGLALARH`). Cada columna corresponde a un aminoácido. Los valores son la probabilidad estimada por MEME a partir de las 116 secuencias del descubrimiento.

**Posición crítica para la discusión (fila 11, "A" del MNPAV):**

| A | C | D | E | F | G | H | I | K | **L** | M | N | P | Q | R | **S** | T | V | W | Y |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **90,5 %** | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | **0 %** | 0 | 0 | 0 | 0 | 0 | **9,5 %** | 0 | 0 | 0 | 0 |

MEME aprende que esa posición admite A (90,5 %) o S (9,5 %), pero nunca L. Las isoformas SIP2;2 (con NPL) caen fuera del rango aceptado y por eso no se detectan como instancias del motivo M2 ni en el descubrimiento ni en el escaneo posterior.

Las dos isoformas SIP1 con NPT no contribuyen a la matriz (su T no aparece en la columna T de la fila 11), pero sí pasan el escaneo posterior gracias a la conservación de las demás 20 posiciones del motivo. Las siete isoformas SIP1;3 con NPS sí están incluidas en el descubrimiento original (contribuyen al 9,5 % de S).

## 3. Concordancia con la literatura

La degeneración del motivo NPA del bucle B en SIP es un rasgo distintivo descrito en otras plantas:

- **Wudick et al. (2009)**, *New Phytologist*: SIP es la única subfamilia plant-AQP que presenta el motivo NPA degenerado en su primer loop (NPT, NPL, NPC u otras variantes en lugar de NPA estricto).
- **Nicolas-Espinosa y Carvajal (2022)**, *The Plant Genome*: BoiSIP1;1 de brócoli presenta NPT/NPA, en línea con el patrón aquí descrito en *F.* × *ananassa*.

El presente trabajo aporta el patrón completo en una especie alo-octoploide: las 12 SIPs muestran las tres variantes (NPL/NPT/NPS) estructuradas por clado homeólogo, sin ningún caso de NPA canónico en el bucle B.

## 4. Reproducción

El cálculo de frecuencias de motivo por subfamilia (que produce el heatmap de la Figura 5) se obtiene ejecutando `MEME_FINAL/analisis_motivos_meme.py` sobre `MEME_FINAL/ALL_AQP.txt`. El script `MEME_FINAL/_calcular_frecuencias.py` reproduce las cifras concretas reportadas en el apartado 6.2.3.
