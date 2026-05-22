# Visores interactivos en GitHub Pages

Los visores HTML del TFG se publican como una pequeña web estática servida por
**GitHub Pages** desde la carpeta `docs/` de este repositorio. No hace falta
servidor propio ni base de datos: son archivos HTML autocontenidos.

## Qué se publica

| Archivo en `docs/` | Contenido | Lo genera |
|---|---|---|
| `index.html` | Página de portada que enlaza a los visores | (escrito a mano) |
| `efp_homeologos.html` | Atlas eFP de expresión por tejido de los 32 grupos homeólogos | `scripts/5.5_rna_seq/5.5.3_homeologos/15_homeolog_efp_viewer.py` |
| `pca_interactivo.html` | Explorador del PCA de las 121 acuaporinas | `scripts/5.2_5.3_homologia_curacion/profiling_final_integrated.py` |
| `.nojekyll` | Le dice a Pages que sirva los HTML tal cual (sin procesarlos con Jekyll) | (archivo vacío) |

> El visor eFP por gen individual se retiró por estar obsoleto; su script
> (`10_generate_efp_viewer.py`) se conserva en el repositorio pero su salida ya
> no se publica.

Los tres comparten un mismo estilo académico (fondo marfil, tipografía serif
Fraunces y la paleta oficial de subfamilias) para que la portada, el eFP y el
PCA se vean como una sola pieza.

## Cómo activar GitHub Pages (una sola vez)

1. Sube el repositorio a GitHub (ver `docs/PENDIENTE.md`, sección 1).
2. En GitHub, entra en **Settings → Pages**.
3. En **Source**, elige **Deploy from a branch**.
4. Selecciona la rama **`main`** y la carpeta **`/docs`**. Pulsa **Save**.
5. Espera ~1 minuto. La web quedará disponible en:
   `https://paredes0.github.io/tfg-aquaporinas/`

A partir de ahí, cada `git push` a `main` que toque `docs/` actualiza la web
automáticamente.

## Cómo regenerar los visores

Si cambian los datos o el código, vuelve a generar los HTML y cópialos a `docs/`:

```bash
# Atlas eFP de homeólogos  ->  results/efp_viewer_homeologs.html
python scripts/5.5_rna_seq/5.5.3_homeologos/15_homeolog_efp_viewer.py
cp results/efp_viewer_homeologs.html docs/efp_homeologos.html

# Explorador PCA  ->  results/profiling_aqp_motifs_final/PCA_INTERACTIVO_FINAL.html
python scripts/5.2_5.3_homologia_curacion/profiling_final_integrated.py
cp results/profiling_aqp_motifs_final/PCA_INTERACTIVO_FINAL.html docs/pca_interactivo.html
```

(Ambos scripts leen sus entradas de `data/` vía `scripts/common/config.py`, así
que funcionan tras un `git clone` sin configurar rutas.)

## Comprobación local antes de publicar

Pages sirve los archivos por HTTP, no como `file://`. Para verlos igual en local:

```bash
cd docs
python -m http.server 8000
# abre http://127.0.0.1:8000 en el navegador
```
