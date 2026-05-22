# Pruebas cuantitativas y visuales del descarte de las 5 acuaporinas polémicas

> Documento consolidado para el Anexo del TFG. Aporta evidencia numérica y visual del descarte de las 5 secuencias que conservan 6 (o casi 6) hélices transmembrana nominales y por tanto requieren más justificación que un simple recuento topológico.
>
> **Fuentes de datos cruzadas**:
> - `tabla_Aquaporinas_traduccion.tabular` — veredicto, TMHs, motivos NPA, longitud, pI, Mw.
> - `MEME_FINAL/ALL_AQP.txt` — motivos MEME globales.
> - `PCA_Coordenadas_Finales.csv` — coordenadas PC1, PC2 y Mahalanobis al cuadrado.
> - `fxa_final.iqtree`, `fxa_final.treefile` — árbol exploratorio (304 secuencias = 144 candidatas + 3 controles + 157 referencias).
> - `audit_5_polemicas.fasta` — alineamiento MAFFT E-INS-i (offset 0,23) de las 121 funcionales + 5 polémicas (126 secuencias × 463 columnas).
>
> **Umbrales estadísticos**:
> - χ²(2 gl, α=0,05) = 5,991 (umbral PCA-Mahalanobis).
> - Longitudes de rama terminal del árbol exploratorio: mediana = 0,0323 ; P90 = 0,3143 ; P95 = 0,5104 ; P99 = 0,8555.
> - Cobertura de columnas core del alineamiento (≥ 80 % conservadas en las 121 funcionales): umbral inferior funcional = 88,3 % ; P5 funcional = 94,4 % ; mediana funcional = 99,6 %.

---

## 1. Tabla resumen ejecutivo

| ID en árbol | Veredicto tabla | TMHs | NPA-B alineamiento | NPA-E alineamiento | Long. real | Rama terminal | Soportes (SH-aLRT / aBayes / UFBoot) | PCA d² | Cobertura core | Causa principal del descarte |
|---|---|---:|---|---|---:|---:|---|---:|---:|---|
| `Fxa6Bg00715-maker` | MAKER_GFF3 | 6 | ✓ NPAV | ✓ NPARS | 225 aa | 0,0134 | 0 / 0,333 / 58 | dentro elipse | **79,7 %** | **91 aa de deleción interna** entre NPA-B y NPA-E (cols 216-306) |
| `Fxa6Cg01391-maker` | MAKER_GFF3 | 6 | ✓ NPAV | ✓ NPARS | 297 aa | **0,3286** | 0 / 0,333 / 95 | dentro elipse | 96,1 % | **Clado equivocado en árbol** (PIP en lugar de NIP) + **rama 10× la mediana** + 8 bloques de gap interno (134 cols totales) |
| `Fxa5Bg03706-partial` | AMBAS_MAL | 6 | ✓ NPAV | ❌ **8 cols gap** | 242 aa | 0,1794 | 96 / 1 / 100 | **33,03** | 90,0 % | **Pérdida total del motivo NPA-E** + truncamiento N-terminal de 85 aa + PCA outlier extremo |
| `Fxa3Ag00841-partial` | AMBAS_MAL | **8** | ✓ NPAI | ✓ NPARS | 254 aa | 0,0287 | 0 / 0,333 / 66 | n/d | 89,2 % | **Topología no canónica (8 TMHs)**: extensión N-terminal de 104 aa antes del cuerpo canónico |
| `Fxa5Ag03930-partial` | MANUAL_CURATED | 6* | ❌ **toda en gap** | ✓ NPARS | 178 aa (árbol) | 0,0297 | 96 / 1 / 100 | **22,45** | **67,5 %** | **Pérdida total del motivo NPA-B** + truncamiento severo (178 aa) + PCA outlier |

\* Las 6 TMHs reportadas por DeepTMHMM corresponden a la versión GFF3_FALLBACK de 220 aa; la versión efectivamente usada en el árbol y en este alineamiento tiene 178 aa y carece del N-terminal completo.

**Hallazgo clave**: cada una de las 5 secuencias **falla al menos un criterio cuantitativo objetivo** trazable a un archivo de datos. Dos de las cinco (Fxa5Bg03706 y Fxa5Ag03930) pierden literalmente uno de los dos motivos NPA canónicos, lo que las descalifica como canales funcionales independientemente de cualquier otro filtro.

---

## 2. Posición de los motivos NPA en el alineamiento MAFFT

Construido el consenso sobre las 121 funcionales:

| Motivo | Columnas en el alineamiento | Consenso | Conservación media |
|---|---|---|---|
| **NPA loop B** (asa B citoplasmática) | 192–194 | `GGHINPAVTFGLA` | 93 % |
| **NPA loop E** (asa E extracelular) | 371–373 | `GASMNPARSFGPAV` | 98 % |

Estos dos motivos forman el **filtro selectivo dual** del canal y son la signatura funcional canónica de las acuaporinas (Wang et al., 2020).

---

## 3. Análisis individual con captura del alineamiento

### 3.1. `Fxa6Bg00715-maker` — Deleción interna masiva entre los dos NPA

**Veredicto tabla**: NIP MAKER_GFF3 ; **TMHs**: 6 ; **Longitud**: 225 aa ; **PCA d²**: dentro de elipse.

**Posición en árbol exploratorio**: dentro del clado FaNIP1 cromosoma 6 (junto a `FaNIP1_1`, `Fxa6Ag00747`, `Fxa6Dg00623`), con rama corta (0,0134) pero soportes locales nulos (SH-aLRT = 0 %, aBayes = 0,333).

**Hallazgo MAFFT**: **91 aa de gap continuo entre las columnas 216 y 306**, en la región intermedia entre los dos motivos NPA donde residen las hélices TM2, TM3 y el bucle C. La proteína conserva los dos NPA pero el cuerpo central del canal está vacío.

- 9 bloques internos ≥ 5 aa.
- Suma total de gap interno: 219 columnas (47,3 %).
- Cobertura sobre columnas core: **79,7 %** — 8,6 puntos por debajo de la peor secuencia funcional retenida (88,3 %).

**Veredicto**: las 6 TMHs reportadas por DeepTMHMM cuentan hélices residuales del extremo C; el alineamiento revela que falta la mitad central del canal. Secuencia objetivamente más fragmentada que cualquiera de las 121 funcionales.

**Visualización**: ver `capture_03_Fxa6Bg00715_huecocentral.png` — el contraste con las 3 NIPs canónicas es evidente: las 3 referencias muestran residuos en toda la región 216-306; Fxa6Bg00715-maker está prácticamente vacía.

---

### 3.2. `Fxa6Cg01391-maker` — Clado equivocado en el árbol filogenético

**Veredicto tabla**: NIP MAKER_GFF3 ; **TMHs**: 6 ; **Longitud**: 297 aa ; **PCA d²**: dentro de elipse NIP.

**Posición en árbol exploratorio**: **dentro del clado PIP2;5 del cromosoma 6**, junto a `Fxa6Ag01549`, `Fxa6Cg01396`, `Fxa6Ag01550`, `Fxa6Cg01397`, `Fxa6Bg01475`, `Fxa6Bg01476`. **La tabla la asigna a NIP**; el árbol filogenético la coloca en PIP. La incongruencia es prueba directa de problema estructural.

**Hallazgo MAFFT**: cobertura core 96,1 % (dentro del rango funcional), pero **34 aa de gap interno en cols 60-93** y otros bloques significativos (27 aa cols 217-243, 16 aa cols 283-298), sumando 134 columnas en gap interno (28,9 % de la secuencia). Conserva ambos NPA (NPA-B en col 192 y NPA-E en col 371) con variantes aceptables.

**Métricas filogenéticas anómalas**:
- **Longitud de rama terminal**: 0,3286 — **frente a la mediana del árbol de 0,0323, es 10× mayor** (percentil 90 del árbol entero).
- **Soporte SH-aLRT = 0 %** y **aBayes = 0,333** — la posición filogenética no es estadísticamente significativa.
- **Clado de destino incongruente** con su anotación tabular.
- **BLAST contra UniProt**: mejor hit `P43286|PIP21_ARATH` (AtPIP2;1) — el BLAST también la asocia a PIP, no a NIP.

**Decisión sobre Fxa6Cg01391 — argumentación del descarte**:

Este es el caso más sutil de los cinco. No falla los filtros estructurales primarios (6 TMHs, dos NPA, cobertura razonable, PCA dentro de elipse), pero falla cuando se aplica el filtro filogenético. Las opciones son tres:

| Opción | Argumento a favor | Argumento en contra |
|---|---|---|
| **A. Mantener como NIP** | La tabla y el archivo de motivos MEME la asignan a NIP. | El árbol la coloca en PIP con rama anómala; el BLAST también la asocia a PIP. Las dos fuentes filogenéticas convergen en que NO es NIP. |
| **B. Reclasificarla a PIP atípica** | El árbol y el BLAST coinciden en PIP. Conserva 6 TMHs y NPA. | Procede de MAKER_GFF3 (anotación que fue descartada en bloque por criterio metodológico). La rama anómala (0,33) y soportes nulos (SH-aLRT = 0) indican que **su posición exacta dentro de PIP tampoco es defendible**. La reclasificación sería arbitraria sin re-anotación manual del modelo. |
| **C. Descartar como quimera génica** ✅ | (1) Procede de MAKER_GFF3, que se descartó en bloque por criterio metodológico uniforme. (2) Incongruencia funcional: tabla NIP / árbol PIP / BLAST PIP → la anotación funcional es ambigua. (3) Soportes locales nulos: ni siquiera con la asignación PIP el árbol resuelve su posición. (4) Rama en P90 del árbol: divergencia anómala. | Más conservador: se pierde una secuencia que podría ser funcional. Pero el principio de uniformidad metodológica prevalece sobre el rescate selectivo. |

**Recomendación**: **mantener el descarte como quimera génica** (opción C). Es el único caso de las 5 que solo se caza por el filtro filogenético — y precisamente esto justifica que el TFG argumente la complementariedad de los tres filtros (topológico, MAFFT y multivariante) más el árbol como cuarto filtro estructural-evolutivo. Mencionarla explícitamente en el TFG **refuerza** el discurso metodológico en lugar de debilitarlo.

**Visualización**: ver `capture_06_Fxa6Cg01391_gaps.png` — muestra el gap de 34 aa en cols 60-93 frente a las tres NIPs canónicas, que sí tienen residuos en esa región.

---

### 3.3. `Fxa5Bg03706-partial` — Pérdida total del motivo NPA-E

**Veredicto tabla**: NIP AMBAS_MAL ; **TMHs**: 6 ; **Longitud**: 242 aa ; **PCA d²**: 33,03.

**Posición en árbol exploratorio**: dentro de un subclado de parciales NIP1 (junto a `Fxa5Ag03930-partial` y `Fxa5Cg03343-partial`), con soporte alto del nodo parental (96/1/100) pero formando un subclado periférico de fragmentos.

**Hallazgo MAFFT — el dato decisivo**: en las columnas 367-374 del alineamiento, **donde las 121 funcionales conservan la firma `GASMNPAR`, Fxa5Bg03706 tiene 8 columnas consecutivas en gap**. El motivo NPA-E está físicamente ausente.

La etiqueta `motivo_E = SPV` que aparece en `tabla_Aquaporinas_traduccion.tabular` se generó por un alineador menos sensible que el MAFFT E-INS-i: ese algoritmo intentó "rescatar" la secuencia local SPV que aparece más adelante y la asignó a la posición del NPA-E. Con el alineamiento profundo del MAFFT con offset 0,23, queda claro que la región del motivo está en gap.

**Métricas adicionales**:
- **PCA d² = 33,03** frente al umbral 5,991 — **5,5× por encima del umbral**.
- 85 aa de gap N-terminal (cols 6-90).
- Cobertura core 90 % — por debajo del P5 funcional (94,4 %).

**Veredicto**: secuencia con topología canónica (6 TMHs) pero **sin filtro selectivo C-terminal**. No puede funcionar como canal de agua sin la firma NPA dual.

**Visualización**: ver `capture_02_NPA_E.png` — Fxa5Bg03706-partial tiene 8 celdas blancas (gap) en la región del NPA-E mientras todas las demás muestran NPARS.

---

### 3.4. `Fxa3Ag00841-partial` — Topología no canónica (8 TMHs)

**Veredicto tabla**: NIP AMBAS_MAL ; **TMHs**: 8 ; **Longitud**: 254 aa ; **PCA d²**: n/d.

**Posición en árbol exploratorio**: clado FaNIP1 cromosoma 6 (junto a las MAKER_GFF3 problemáticas y `Fxa3Dg00714-partial-maker`), rama corta (0,0287) pero soportes locales nulos (SH-aLRT = 0 %).

**Hallazgo MAFFT**: la secuencia empieza solo en la **columna 105 del alineamiento**, lo que significa que **la versión sin alinear contiene 104 aa adicionales** respecto al consenso canónico — exactamente los aminoácidos necesarios para acomodar las dos hélices TM extra que DeepTMHMM detecta.

- 6 bloques internos ≥ 5 aa, todos moderados (≤ 18 aa).
- Cobertura core 89,2 % — justo en el límite del rango funcional.
- Conserva ambos NPA (NPA-B en col 192 con variante NPAI, NPA-E en col 371 canónico).

**Veredicto**: secuencia con el cuerpo central correcto pero con un fragmento extra de 104 aa N-terminales que añade dos hélices transmembrana espurias. Probable error de anotación del modelo génico (exones N-terminales incorrectamente fusionados) o región repetitiva mal interpretada.

**Visualización**: ver `capture_05_Fxa3Ag00841_Nterm_extra.png` — Fxa3Ag00841-partial empieza varias columnas después que las referencias NIP, mostrando que su zona alineable comienza tarde respecto a las funcionales.

---

### 3.5. `Fxa5Ag03930-partial` — Pérdida total del motivo NPA-B + truncamiento severo

**Veredicto tabla**: NIP MANUAL_CURATED ; **TMHs**: 6 (versión GFF3_FALLBACK de 220 aa) ; **Longitud en árbol**: 178 aa ; **PCA d²**: 22,45.

**Discrepancia de longitud**: la tabla registra 220 aa para esta entrada (versión GFF3_FALLBACK), pero la secuencia que entró al árbol exploratorio y al MAFFT tiene **únicamente 178 aa** — los primeros ≈ 42 residuos del N-terminal están truncados respecto a la versión GFF3.

**Posición en árbol exploratorio**: subclado de parciales NIP1, junto a `Fxa5Bg03706-partial` y `Fxa5Cg03343-partial`. Soporte del nodo parental 96/1/100.

**Hallazgo MAFFT — el dato decisivo**: en las columnas 184-204 del alineamiento, **donde las 121 funcionales conservan `HISGAHFNPAVTIAFA...`, Fxa5Ag03930 tiene 21 columnas consecutivas en gap**. El motivo NPA-B está físicamente ausente.

**Métricas adicionales**:
- **PCA d² = 22,45** frente al umbral 5,991 — **3,7× por encima del umbral**.
- 8 bloques internos ≥ 5 aa, varios muy grandes: 85 aa (cols 6-90), 67 aa (cols 106-172), 35 aa (cols 177-211).
- Cobertura sobre columnas core: **67,5 %** — **20,8 puntos por debajo del mínimo funcional** (88,3 %). Es objetivamente la secuencia más fragmentada del análisis.

**Veredicto**: la secuencia analizada filogenéticamente (178 aa) es esencialmente un fragmento de los dos tercios C-terminales del canal: conserva NPA-E pero carece de NPA-B y de los dominios transmembrana 1-3. La etiqueta "6 TMHs" de DeepTMHMM se aplicó a la versión GFF3_FALLBACK de 220 aa; sobre los 178 aa efectivos, la secuencia no satisface el filtro topológico canónico.

**Visualizaciones**:
- `capture_01_NPA_B.png` — Fxa5Ag03930-partial tiene su fila completamente blanca (gap total) en toda la región del NPA-B.
- `capture_04_Nterm_truncado.png` — la magnitud del truncamiento N-terminal es evidente comparada con las tres NIPs canónicas.

---

## 4. Distribución de las 5 polémicas en métricas clave del árbol exploratorio

Para contextualizar las cifras:

| Estadístico del árbol exploratorio (304 hojas) | Valor |
|---|---:|
| Mediana de longitud de rama terminal | 0,0323 |
| Media | 0,1163 |
| Percentil 90 | 0,3143 |
| Percentil 95 | 0,5104 |
| Percentil 99 | 0,8555 |
| Máximo | 2,1003 |

**Comparativa de las 5**:

| ID | Rama terminal | Percentil aproximado | Categoría |
|---|---:|---|---|
| Fxa6Bg00715-maker | 0,0134 | < P40 | corta |
| Fxa3Ag00841-partial | 0,0287 | ≈ P50 | corta |
| Fxa5Ag03930-partial | 0,0297 | ≈ P50 | corta |
| Fxa5Bg03706-partial | 0,1794 | ≈ P75 | moderada |
| **Fxa6Cg01391-maker** | **0,3286** | **≈ P90** | **anómala** |

Solo Fxa6Cg01391 sería identificable únicamente por longitud de rama. Las otras 4 se detectan por **otras métricas independientes** (NPA ausentes en MAFFT, soportes locales nulos, topología no canónica, PCA outlier severo). Esta diversidad de fallos prueba que el panel de filtros del TFG (DeepTMHMM + MAFFT + PCA + árbol filogenético) **detecta tipos de incongruencia distintos y complementarios**.

---

## 5. Cobertura efectiva sobre las 231 columnas core del alineamiento

Las **231 columnas core** son aquellas donde ≥ 80 % de las 121 funcionales tienen residuo (no gap) — esencialmente los segmentos transmembrana, los bucles que portan los motivos NPA y las regiones de empaquetamiento estructural.

| Conjunto | Mediana | P5 | Mínimo |
|---|---:|---:|---:|
| 121 funcionales | 99,6 % | 94,4 % | 88,3 % |

**Comparativa de las 5 polémicas**:

| ID | Cobertura | Lectura |
|---|---:|---|
| Fxa6Cg01391-maker | 96,1 % | Dentro del rango funcional (cazada por árbol) |
| Fxa5Bg03706-partial | 90,0 % | Por debajo del P5 funcional |
| Fxa3Ag00841-partial | 89,2 % | Rozando el mínimo funcional |
| Fxa6Bg00715-maker | 79,7 % | 8,6 puntos por debajo del mínimo funcional |
| Fxa5Ag03930-partial | **67,5 %** | **20,8 puntos por debajo** — la más fragmentada |

---

## 6. Veredicto consolidado

Cada una de las 5 polémicas se sustenta en **al menos dos métricas cuantitativas independientes** trazables a archivos de origen. Ninguna se descarta arbitrariamente.

| ID | Métrica principal | Métrica complementaria | Filtro que la caza |
|---|---|---|---|
| **Fxa6Bg00715-maker** | Cobertura core 79,7 % (mínimo func: 88,3 %) | Gap interno de 91 aa entre los dos NPA | **MAFFT** |
| **Fxa6Cg01391-maker** | Rama terminal 0,3286 (P90 del árbol) | Clado erróneo en árbol (PIP vs NIP tabular) | **Árbol filogenético** |
| **Fxa5Bg03706-partial** | NPA-E ausente (8 cols gap en el motivo) | PCA d² = 33 (5,5× umbral) | **MAFFT** + **PCA** |
| **Fxa3Ag00841-partial** | TMHs = 8 (topología no canónica) | Extensión N-terminal de 104 aa en alineamiento | **DeepTMHMM** + **MAFFT** |
| **Fxa5Ag03930-partial** | NPA-B ausente (21 cols gap en el motivo) | Cobertura core 67,5 % (mínimo func: 88,3 %) | **MAFFT** + **PCA** |

---

## 7. Capturas del alineamiento incluidas

| Archivo | Contenido | Polémica(s) ilustrada(s) |
|---|---|---|
| `capture_01_NPA_B.png` | Región del motivo NPA-B (cols 184-204) con las 3 NIPs canónicas como referencia y las 5 polémicas | Fxa5Ag03930 (gap total) |
| `capture_02_NPA_E.png` | Región del motivo NPA-E (cols 363-385) | Fxa5Bg03706 (8 cols gap sobre el motivo) |
| `capture_03_Fxa6Bg00715_huecocentral.png` | Región central entre los dos NPA (cols 216-310) | Fxa6Bg00715 (91 aa gap continuo) |
| `capture_04_Nterm_truncado.png` | Región N-terminal completa (cols 1-200) | Fxa5Ag03930 + Fxa5Bg03706 (truncamiento severo) |
| `capture_05_Fxa3Ag00841_Nterm_extra.png` | Inicio del cuerpo canónico (cols 91-170) | Fxa3Ag00841 (zona alineable empieza tarde) |
| `capture_06_Fxa6Cg01391_gaps.png` | Región N-proximal (cols 56-145) | Fxa6Cg01391 (gap de 34 aa cols 60-93) |
| `capture_07_global_5polemicas.png` | Vista global de las 5 polémicas (alineamiento completo, 463 cols, NPA-B y NPA-E destacados) | Las 5 a la vez |

---

## 8. Dataset entregable

- **`audit_5_polemicas/dataset_121_funcionales_mas_5_polemicas.fasta`** — 126 secuencias (121 funcionales + 5 polémicas) sin alinear.
- **`audit_5_polemicas/5_polemicas.fasta`** — solo las 5 polémicas, sin alinear.
- **`audit_5_polemicas/audit_5_polemicas.fasta`** — alineamiento MAFFT E-INS-i, offset 0,23, 463 columnas. Es el archivo usado para todo el análisis cuantitativo aquí.

Headers de las 5 en los FASTA:
- `>Fxa6Bg00715-maker` (sufijo `-maker` = MAKER_GFF3)
- `>Fxa6Cg01391-maker` (sufijo `-maker` = MAKER_GFF3)
- `>Fxa5Bg03706-partial` (sufijo `-partial` = AMBAS_MAL)
- `>Fxa3Ag00841-partial` (sufijo `-partial` = AMBAS_MAL)
- `>Fxa5Ag03930-partial` (sufijo `-partial` = MANUAL_CURATED descartada en revisión PCA)
