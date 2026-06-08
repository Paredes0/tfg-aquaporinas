# Sistema de diseño · Deck de defensa TFG (acuaporinoma)

> Reglas objetivas para que **toda** diapositiva use bien el espacio y sea legible a distancia.
> Lienzo de referencia: **1920 × 1080**. Frame útil ≈ 1696 × 876 px (padding 112 lateral, 96 arriba, 108 abajo).
> Estas reglas **mandan sobre cualquier `font-size` inline**. Si una diapo contradice esto, la diapo está mal.

## 1 · Suelos de tamaño (fondo claro)

| Rol | Tamaño | Token / color |
|---|---|---|
| Título de diapo | 64 px serif | `.title` · `--ink` (#14182a) |
| Subtítulo / intro bajo el título | 34–40 px | `.subtitle` · `--ink-soft` |
| Cuerpo / punto clave | **28–32 px · nunca < 24** | `--ink-soft` (#3a4055) |
| Texto secundario / pie dentro de tarjeta | 22–24 px · **nunca < 20** | `--ink-muted` (#6c7080) |
| Etiqueta mono / eyebrow / `card-title` | 18–22 px | acento o `--ink-muted` |
| Chrome `topbar` | 22 px | `--ink-muted` |
| `footer-rule` (citas) | 16 px | `--ink-muted` |
| Cifra hero (`megafig`) | 96–168 px | acento |

Regla rápida: **si dudas, súbelo**. El error por defecto del deck era poner texto a 16–24 px teniendo sitio de sobra.

## 2 · Color y contraste

- **Prohibido `--ink-dim` (#a0a4b4) para texto que se lee.** Solo para hairlines/decoración. (Contraste ~2,3:1 sobre blanco = ilegible.)
- Jerarquía de tinta sobre claro: titulares `--ink` → cuerpo `--ink-soft` → secundario `--ink-muted`.
- Acentos sobre **fondo claro**: usar las variantes oscuras ya definidas — water `#2980b9`, fruit `#c0392b`, leaf `#229954`, amber `#d68910`. (Las brillantes `#4ec5e0/#ff5577/…` son solo para `data-theme="dark"`.)
- **Hex brillantes de subfamilia** (PIP `#E74C3C`, TIP `#3498DB`, NIP `#2ECC71`, SIP `#F39C12`, XIP `#9B59B6`): para **rellenos, puntos, bordes**. Como **texto fino sobre blanco** van en peso 600; el ámbar SIP `#F39C12` casi no se ve en texto → engrosar o usar `--c-amber`.
- Paleta de subfamilia (fills/puntos) es inamovible: concordancia con figuras del TFG.

## 3 · Uso del espacio

- El contenido debe ocupar **~75–85 % de la altura útil** del frame. Nada de medias diapos vacías (caso índice) ni texto apretado contra el borde inferior (caso 5.1 comprimido).
- Diapos "ligeras" (índice, separadoras de resultado): se **escala la tipografía y el aire**, no se deja hueco muerto.
- `.frame` ya centra vertical (`justify-content: center`). Para diapos con mucho contenido, repartir en bandas con `gap` 24–44.
- Una idea grande por diapo de resultado; el detalle, a notas de orador.

## 4 · Estructura repetida (no romper)

- `topbar` arriba (eyebrow + contador `NN / TOTAL`), `footer-rule` abajo (citas + contador). El contador SIEMPRE con denominador = nº total de diapos.
- `data-label="NN · …"` por `<section>`, en orden.
- `#speaker-notes`: array JSON **alineado por índice** con las `<section>`. Si se añade/quita/reordena una diapo → reconstruir el array y renumerar contadores + `data-label`.
- Separadoras de bloque: `data-theme="dark"` + `scene-overlay-bar`.
- Bindings: `data-scene`, `data-embed`, `data-aquaporin-target`, `data-embed-target`. No tocar sin actualizar `deck-orchestrator.js`.
- Color de bloque unificado (water) en índice / roadmap / separadoras — NO rojo en el Bloque IV.

## 5 · Convenciones de contenido (heredadas)

- Registro formal de defensa TFG; sin jerga ni coloquialismos; títulos = títulos, no frases.
- La diapo lleva puntos clave; la explicación llana la da el orador.
- `x` latina (no `×`). Decimales con coma. *F.* x *ananassa* en cursiva.
- Paleta y cifras consolidadas: ver `CLAUDE.md` (no contradecir).

## 6 · Checklist antes de dar una diapo por buena

1. ¿Algún texto < 24 px que sea de leer? → subir.
2. ¿Algún texto en `--ink-dim`? → cambiar a `--ink-muted`/`--ink-soft`.
3. ¿El contenido llena ~80 % de la altura? → si no, escalar.
4. ¿`topbar`/`footer`/contador/`data-label` correctos y coherentes con el total?
5. Overflow Playwright = `[]` (sin desbordes) a 1920×1080.
