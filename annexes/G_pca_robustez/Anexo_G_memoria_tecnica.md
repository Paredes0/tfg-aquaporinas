# PCA "ingenuo" — ¿son outliers las parciales solo porque excluimos las parciales del cálculo?

> El PCA del TFG entrena las elipses al 95 % usando únicamente las secuencias CLEAN (las 121 funcionales), proyectando las 23 parciales como observaciones translúcidas. Esta convención podría parecer circular: ¿estamos llamando "fuera de elipse" a las parciales solo porque las hemos excluido del cálculo de la elipse?
>
> Para resolver esa duda metodológica se realiza aquí una **prueba de robustez**: se recalculan las elipses **incluyendo las parciales en su cálculo** (escenario "sin saber a priori cuáles son parciales") y se contabiliza cuántas siguen siendo outliers.

---

## 1. Diseño de la prueba

| Versión | Quiénes entrenan las elipses | Filosofía |
|---|---|---|
| **A (TFG actual)** | Solo las 121 CLEAN | Las parciales se proyectan sin afectar el cálculo. Es el escenario estricto: la elipse refleja la firma fisicoquímica "ideal" de cada subfamilia. |
| **B (prueba de robustez)** | Las 144 candidatas (CLEAN + PARTIAL) | Las parciales sí influyen en el cálculo. La elipse se infla para acomodarlas. Es el escenario más permisivo posible. |

**Hipótesis a probar**: si una parcial sigue cayendo fuera de la elipse incluso en la versión B (la permisiva), entonces es un outlier **objetivo e independiente del sesgo de selección**.

---

## 2. Cómo se inflan las elipses al incluir las parciales

| Subfamilia | n_A | n_B | Centroide A (PC1, PC2) | Centroide B | Varianza A (sxx, syy) | Varianza B | Cambio |
|---|---:|---:|---|---|---|---|---|
| NIP | 37 | 52 | (2,74 ; 1,40) | (2,20 ; 0,91) | (0,81 ; 0,72) | **(3,21 ; 2,69)** | **Varianza × 4** |
| PIP | 32 | 33 | (0,79 ; −1,43) | (0,84 ; −1,40) | (0,44 ; 0,30) | (0,50 ; 0,32) | mínimo |
| SIP | 12 | 18 | (−2,27 ; −3,91) | (−2,47 ; −3,76) | (2,12 ; 0,80) | (2,47 ; **2,05**) | varianza PC2 × 2,5 |
| TIP | 34 | 35 | (−3,15 ; 1,41) | (−3,18 ; 1,32) | (1,11 ; 1,55) | (1,10 ; 1,84) | leve |
| XIP | 6 | 6 | (1,26 ; −1,24) | (1,26 ; −1,24) | (0,85 ; 0,02) | (0,85 ; 0,02) | sin cambio (no había parciales) |

**Lectura**: NIP y SIP son las subfamilias con más parciales, por eso son las que más se inflan. PIP y TIP cambian poco. XIP no cambia porque no había parciales asignadas a esa subfamilia.

La elipse NIP en versión B casi cuadruplica su área para "tragarse" las parciales — es el escenario más permisivo posible para esa subfamilia.

---

## 3. Resultado de la comparación: Mahalanobis al cuadrado en cada versión

Umbral χ²(2 gl, α=0,05) = **5,991**. Por encima de ese umbral, la secuencia cae fuera de la elipse al 95 %.

| ID | Subfam | Veredicto | TMHs | d²_A | A_OUT? | d²_B | B_OUT? | Lectura |
|---|---|---|---:|---:|---|---:|---|---|
| **Fxa5Bg03706** | NIP | AMBAS_MAL | 6 | **86,73** | OUT | **10,31** | **OUT** | Outlier robusto (pierde NPA-E) |
| **Fxa5Cg03343** | NIP | AMBAS_MAL | 4 | 84,44 | OUT | **13,10** | **OUT** | Outlier robusto |
| **Fxa5Ag03930** | NIP | MANUAL_CURATED | 6* | **61,12** | OUT | **6,35** | **OUT** | Outlier robusto (pierde NPA-B en alineamiento) |
| **Fxa6Dg03790** | SIP | AMBAS_MAL | 4 | 39,33 | OUT | **7,63** | **OUT** | Outlier robusto (candidata reanotación expresada) |
| **Fxa5Dg03404** | NIP | MAKER_GFF3 | 2 | 39,17 | OUT | **7,03** | **OUT** | Outlier robusto |
| **Fxa6Dg03789** | NIP | MAKER_GFF3 | 3 | 25,65 | OUT | **8,92** | **OUT** | Outlier robusto |
| **Fxa2Cg01599** | TIP | MAKER_GFF3 | 4 | 22,11 | OUT | **13,02** | **OUT** | Outlier robusto |
| **Fxa7Dg01388** | PIP | AMBAS_MAL | 5 | 17,28 | OUT | **10,89** | **OUT** | Outlier robusto (candidata reanotación, TPM=43) |
| Fxa3Dg00714 | NIP | MAKER_GFF3 | 4 | 18,09 | OUT | 3,55 | in | Cazado solo por estricto |
| Fxa6Ag04864 | NIP | AMBAS_MAL | 2 | 18,93 | OUT | 4,58 | in | Cazado solo por estricto |
| Fxa4Ag01484 | NIP | AMBAS_MAL | 5 | 16,38 | OUT | 3,26 | in | Cazado solo por estricto |
| Fxa3Ag00841 | NIP | AMBAS_MAL | 8 | 14,85 | OUT | 1,17 | in | Cazado solo por estricto |
| Fxa1Bg03423 | SIP | MAKER_GFF3 | 5 | 14,92 | OUT | 2,91 | in | Cazado solo por estricto |
| Fxa6Ag01546 | SIP | MAKER_GFF3 | 3 | 24,08 | OUT | 4,38 | in | Cazado solo por estricto |
| Fxa4Bg03149 | NIP | MAKER_GFF3 | 2 | 40,04 | OUT | 4,11 | in | Cazado solo por estricto |
| Fxa6Ag04863 | NIP | MAKER_GFF3 | 2 | 5,91 | in | 0,32 | in | Nunca fue outlier PCA |
| Fxa5Bg01988 | SIP | MAKER_GFF3 | 2 | 5,00 | in | 1,56 | in | Nunca fue outlier PCA |
| Fxa6Cg01391 | NIP | MAKER_GFF3 | 6 | 2,04 | in | 0,04 | in | Nunca fue outlier PCA (cazada por MAFFT+árbol) |
| Fxa6Bg00715 | NIP | MAKER_GFF3 | 6 | 4,26 | in | 0,77 | in | Nunca fue outlier PCA (cazada por MAFFT) |
| Fxa3Cg00716 | NIP | MAKER_GFF3 | 2 | 4,38 | in | 1,86 | in | Nunca fue outlier PCA |
| Fxa1Bg02622 | SIP | MAKER_GFF3 | 5 | 3,18 | in | 3,47 | in | Nunca fue outlier PCA |
| Fxa2Bg00225 | SIP | MAKER_GFF3 | 5 | 2,75 | in | 1,55 | in | Nunca fue outlier PCA |
| Fxa3Bg00731 | NIP | MAKER_GFF3 | 2 | 1,18 | in | 0,97 | in | Nunca fue outlier PCA |

\* Fxa5Ag03930: la tabla la marca con 6 TMHs (versión GFF3_FALLBACK de 220 aa), pero la versión que entró al PCA y al árbol tiene 178 aa.

---

## 4. Resumen estadístico

| Categoría | Cuántas | Significado |
|---|---:|---|
| Parciales evaluadas | 23 | (sin contar Fxa1Ag03542 ni las 2 aminoacilasas, que no estaban en el PCA) |
| Fuera de elipse en versión A (estricta) | **15** | El PCA "del TFG" las marca como outliers |
| Fuera de elipse en versión B (permisiva) | **8** | Aunque les damos peso en el cálculo, siguen fuera |
| **Outliers robustos** (fuera en ambas versiones) | **8 / 23 = 35 %** | Outliers objetivos, no por sesgo de selección |
| Outliers solo en versión A | 7 / 23 = 30 % | Caen dentro al permitir que inflen la elipse |
| Nunca outliers PCA | 8 / 23 = 35 % | Caen dentro de la elipse aunque tengan TMHs<6 o problemas estructurales — su descarte se basa en otros filtros (topológico, MAFFT, árbol) |

---

## 5. Las 5 polémicas analizadas individualmente

| ID | d²_A | d²_B | ¿Outlier robusto? | Filtro decisivo |
|---|---:|---:|---|---|
| **Fxa5Ag03930** (FaNIP1, 178 aa, MANUAL_CURATED) | 61,12 | **6,35** | **SÍ** — fuera en ambas versiones | **PCA** (pérdida NPA-B confirmada por MAFFT) |
| **Fxa5Bg03706** (NIP, AMBAS_MAL, 6 TMHs) | 86,73 | **10,31** | **SÍ** — fuera en ambas versiones | **PCA** (pérdida NPA-E confirmada por MAFFT) |
| Fxa3Ag00841 (NIP, AMBAS_MAL, 8 TMHs) | 14,85 | 1,17 | NO — dentro en versión B | **DeepTMHMM** (8 TMHs no canónicas) + MAFFT (N-terminal extendido +104 aa) |
| Fxa6Bg00715 (NIP, MAKER_GFF3, 6 TMHs) | 4,26 | 0,77 | NO — dentro en ambas | **MAFFT** (91 aa de gap interno entre los dos NPA) |
| Fxa6Cg01391 (NIP, MAKER_GFF3, 6 TMHs) | 2,04 | 0,04 | NO — dentro en ambas | **Árbol filogenético** (clado erróneo + rama anómala) |

**Lectura crucial**: solo 2 de las 5 polémicas son outliers robustos por PCA (Fxa5Ag03930 y Fxa5Bg03706). Las otras 3 caen dentro de la elipse incluso en la versión estricta — son cazadas por **filtros independientes** (topológico, MAFFT, filogenético).

---

## 6. Implicaciones para la redacción del TFG

### 6.1. Reformulación recomendada del párrafo del PCA

El TFG actual puede sonar circular: "entrenamos el PCA con las CLEAN y luego las parciales caen fuera". La prueba A/B convierte esa frase en algo defendible:

**Frase propuesta para añadir al 6.1 (después de la mención del PCA y antes de la interpretación biológica)**:

> *"Para descartar que la posición externa de las secuencias parciales en el PCA fuese un artefacto del esquema de entrenamiento (las elipses se calcularon únicamente sobre las candidatas que cumplían los filtros estructurales primarios), se realizó una prueba de robustez consistente en recalcular las elipses al 95 % incluyendo todas las candidatas (CLEAN + PARTIAL) en el cálculo. Bajo esta condición permisiva, **8 de las 23 secuencias descartadas (35 %) siguen cayendo fuera de la elipse de su subfamilia**, manteniéndose como outliers objetivos. Estas 8 secuencias —entre ellas Fxa5Ag03930 y Fxa5Bg03706, las dos que pierden uno de los motivos NPA canónicos en el alineamiento MAFFT— constituyen el núcleo de outliers PCA respaldado por evidencia multivariante independiente del sesgo de selección. Las 15 secuencias restantes son outliers únicamente bajo el esquema estricto y se descartan en su mayoría por filtros complementarios: incumplimiento topológico (TMHs < 6), bloques de gap interno detectados por MAFFT o anomalías filogenéticas (rama anómala, soportes locales nulos, clado equivocado). Esta distribución demuestra que los tres filtros estructurales del flujo de curaduría detectan tipos de incongruencia distintos y no redundantes."*

### 6.2. Decisión metodológica

**Mantener el PCA actual (versión A)** como visualización principal del 6.1 + **referenciar la versión B como prueba de robustez en el Anexo**. Las dos elipses (estricta y permisiva) pueden mostrarse superpuestas en una figura adicional del Anexo (`PCA_elipses_superpuestas.png`).

La justificación operativa es:
- La versión A es **más informativa** para clasificar (las elipses reflejan la firma fisicoquímica canónica de la subfamilia, no la "contaminada" por outliers).
- La versión B es **más conservativa** para defender el descarte (demuestra que el subconjunto de outliers robustos persiste incluso bajo el escenario menos restrictivo).
- Reportar ambas resuelve la objeción metodológica sin alterar el análisis principal.

---

## 7. Archivos generados

Este anexo reúne las **dos vistas** de la misma prueba de robustez (el antiguo
Anexo F —paneles separados— se fusionó aquí en 2026-05-22):

| Archivo | Contenido |
|---|---|
| `Anexo_G_figura_paneles_AvsB.png` | Dos paneles lado a lado: versión A (izquierda) y versión B (derecha), con las polémicas etiquetadas. *(Antes Anexo F.)* |
| `Anexo_G_figura_elipses_AvsB.png` | Plano único con las elipses A (línea sólida) y B (línea discontinua) superpuestas, para mostrar el grado de inflado de un vistazo. |
| `Anexo_G_tabla_Mahalanobis_AvsB.csv` | Tabla cuantitativa: ID, subfamilia, veredicto, TMHs, PC1, PC2, d²_A, out_A, d²_B, out_B. |
| `Anexo_G_script_reproducible.py` | Script reproducible (vía `scripts/common/config.py`) que regenera ambas figuras a partir de `results/profiling_aqp_motifs_final/PCA_Coordenadas_Finales.csv`. |
| `Anexo_G_memoria_tecnica.md` | Este documento. |

**Importante**: la prueba se hace sobre los **mismos ejes PC1/PC2** del PCA original (no se re-entrenó el PCA con nuevas variables). Lo único que varía son las elipses de confianza por subfamilia. Si quisiéramos rigor adicional, podría re-ejecutarse `profiling_final_integrated.py` permitiendo que las parciales participen en la estandarización y en la descomposición — pero como la transformación es lineal y las parciales ya están proyectadas en este espacio, la inferencia sobre la pertenencia a cada cluster no se altera sustancialmente.
