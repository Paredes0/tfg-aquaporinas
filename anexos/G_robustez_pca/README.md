# Anexo G — Robustez del análisis de componentes principales (PCA)

Citado en §6.1 del TFG. Prueba de robustez del PCA fisicoquímico (Figura 6 del cuerpo): el PCA del TFG entrena las elipses de confianza al 95 % usando solo las 121 secuencias funcionales y proyecta las 23 parciales como observaciones. Este anexo comprueba que las parciales atípicas lo siguen siendo **aunque se incluyan en el cálculo de las elipses** (escenario permisivo), descartando un sesgo circular de selección.

## Contenido

| Archivo | Contenido |
|---|---|
| `Anexo_G_memoria_tecnica.md` | Memoria: diseño de la prueba (versión A vs B), inflado de las elipses por subfamilia y resultado. |
| `Anexo_G_tabla_Mahalanobis_AvsB.csv` | Distancia de Mahalanobis al cuadrado de cada candidata en las dos versiones. |
| `Anexo_G_figura_elipses_AvsB.png` | Elipses al 95 % superpuestas (A: solo funcionales; B: con parciales). |
| `Anexo_G_figura_paneles_AvsB.png` | Las dos versiones en paneles separados. |
| `Anexo_G_script_reproducible.py` | Regenera las figuras y la tabla vía `scripts/common/config.py`. |

## Diseño de la prueba

- **Versión A (la del TFG):** las elipses se entrenan solo con las 121 funcionales (CLEAN).
- **Versión B (prueba de robustez):** las 144 candidatas (CLEAN + parciales) entrenan las elipses, que se inflan para acomodarlas (escenario más permisivo posible).

Si una parcial sigue cayendo fuera de la elipse incluso en la versión B, es un outlier **objetivo e independiente** del criterio de selección. El umbral es la χ² con 2 grados de libertad (α = 0,05) → d² = 5,991.
