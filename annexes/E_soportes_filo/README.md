# Anexo E — Distribución de soportes filogenéticos (PENDIENTE)

Este anexo está planificado pero su contenido aún no se ha generado. Citado en §6.2.2 del TFG.

## Contenido planificado

- `Anexo_E_tabla_277_nodos_soportes.csv` — un nodo por fila, columnas SH-aLRT (1000 iter), UFBoot (1000 iter), aBayes + indicadores de triple soporte.
- `Anexo_E_figura_histograma_soportes.png` — distribución de los 3 estadísticos con líneas de umbrales (UFBoot ≥ 95%, SH-aLRT ≥ 80%, aBayes ≥ 0,95).
- `Anexo_E_tabla_nodos_subfamilias.csv` — subtabla con los nodos que delimitan PIP, TIP, NIP, SIP, XIP.
- `Anexo_E_script_reproducible.py` — script ete3/dendropy que regenera tabla y figura desde el `.treefile` final.

## Bloqueo

El `.treefile` final (con anotación tripartita SH-aLRT/aBayes/UFboot) vive en `Z:\work\RNA-seq_test\arbol_acuaporinas_2_bueno_sin_parciales.treefile`. Cuando el árbol esté accesible localmente, se ejecuta el script Python y se materializa el contenido. Tarea diferida a sub-proyecto independiente o resuelta antes del release final.
