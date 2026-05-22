# Anexo I — Distribución de soportes filogenéticos

Citado en §6.2.2 del TFG como respaldo cuantitativo de la robustez del árbol filogenético definitivo (Q.PLANT+R6, 282 hojas, 430 sitios, log L = −45.149,26).

## Contenido

| Archivo | Qué es |
|---|---|
| `Anexo_I_script_reproducible.py` | Script que recorre el `.treefile` final, extrae los tres soportes por nodo, genera la tabla, el histograma y la subtabla de subfamilias. |
| `Anexo_I_tabla_277_nodos_soportes.csv` | Un nodo interno por fila (277 en total): SH-aLRT (%), aBayes, UFBoot (%) e indicadores de superación de umbral + triple soporte alto. |
| `Anexo_I_figura_histograma_soportes.png` / `.pdf` | Histograma de los tres estadísticos sobre los 277 nodos, con líneas verticales en los umbrales convencionales. |
| `Anexo_I_tabla_nodos_subfamilias.csv` | Soporte del nodo (MRCA) que define cada subfamilia canónica. |

## Resultados

Sobre los **277 nodos internos** del árbol definitivo:

- **SH-aLRT ≥ 80 %**: 198 nodos.
- **aBayes ≥ 0,95**: 196 nodos.
- **UFBoot ≥ 95 %**: 122 nodos.
- **Triple soporte alto** (los tres umbrales a la vez): **116 nodos (41,9 %)**.

Los nodos que delimitan las subfamilias canónicas se encuentran entre los de máxima confianza: PIP, NIP y SIP presentan triple soporte máximo (SH-aLRT 100 / aBayes 1 / UFBoot 100) y TIP forma un clado monofilético (SH-aLRT 97,6 / aBayes 1 / UFBoot 69). La subfamilia XIP queda como grupo basal: su nodo de origen (MRCA) coincide con la raíz del árbol no enraizado, un punto que por definición no porta valor de soporte; su monofilia se aprecia en el árbol completo (Figura 7 del TFG).

## Reproducción

```bash
# Opcional: apuntar a tu copia del .treefile final
export TFG_TREEFILE=/ruta/a/arbol_acuaporinas.treefile
python Anexo_I_script_reproducible.py
```

El script localiza el árbol vía la variable `TFG_TREEFILE` o, en su defecto, en las rutas locales conocidas. El árbol de referencia es el de 282 hojas (`final_without_partials`), que produce 277 nodos internos y coincide con los parámetros citados en §6.2.2.

## Notas técnicas

- Formato de soporte en el `.treefile`: `)SH-aLRT/aBayes/UFboot:branch_length` (p. ej. `)98.9/1/100:0.23`), tal como lo anota IQ-TREE con `-alrt 1000 -abayes -bb 1000`.
- El mapeo de cada hoja de *Fragaria* a su subfamilia se toma de `homeolog_groups.tsv` (columna `subfamily`, de la que se extrae la subfamilia canónica); las hojas de referencia (At/Os/Md/Hb) se asignan por el nombre.

### Sobre las dos versiones del árbol (282 vs 281 hojas)

Existen dos versiones de la reconstrucción definitiva que difieren en **una sola hoja**:

- **282 hojas → 277 nodos internos** (`final_without_partials/arbol_acuaporinas.treefile`): incluye la rama de la secuencia parcial `Fxa6Dg03790`. Es la versión del `.iqtree` del que proceden los parámetros citados en §6.2.2 (log L = −45.149,26) y los estadísticos de soporte (277 nodos, 116 con triple soporte alto, 41,9 %; medianas UFBoot 91 / SH-aLRT 92,2 / aBayes 1,000).
- **281 hojas → 276 nodos** (`Z:\work\RNA-seq_test\arbol_acuaporinas_2_bueno_sin_parciales.treefile`): la misma reconstrucción tras podar `Fxa6Dg03790`. Sus estadísticos son prácticamente idénticos (276 nodos, 116 triple soporte, 42,0 %; medianas UFBoot 91,5 / SH-aLRT 92,2 / aBayes 1,000).

**Este anexo usa la versión de 282 hojas (277 nodos)** porque es la que respalda exactamente las cifras citadas en el cuerpo del TFG (§6.2.2). La secuencia parcial `Fxa6Dg03790` se excluyó del conjunto curado de 121 funcionales pero su rama se conservó en el árbol del que se calcularon los soportes, ya que podarla no altera de forma apreciable los resultados (la diferencia entre ambas versiones es de un nodo y menos de una décima en los porcentajes). Esta decisión es coherente con la nota del proyecto: "281 secuencias / parámetros del `.iqtree` de 282 que es esencialmente la misma reconstrucción".
