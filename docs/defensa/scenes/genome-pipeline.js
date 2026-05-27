// scenes/genome-pipeline.js — Bloque I, slide 05.
//
// V4 (2026-05-26): correccion de timing y continuidad visual.
//   * Las cartas de especies SIEMPRE aterrizan antes de que arranquen las
//     flechas/conexiones, y el grid de cromosomas aparece junto con los
//     conectores de salida (no antes ni despues).
//   * Timeline ~6,5 s, con cascada chr-a-chr de los puntos.
//   * Conector Exonerate → grid: ahora un abanico de 5 lineas finas que parten
//     del borde derecho de la carta Exonerate y abren hacia la fila superior
//     del grid (en lugar del unico arco curvo a una sola cromatica vertical
//     antigua).
//   * Bajo el grid: nuevo bloque de leyenda con (a) titular "129 candidatas en
//     28 cromosomas (7 × 4 subgenomas)", (b) aclaracion del cariotipo 2n=8x=56,
//     (c) detalle Chr 6 / Subg. A, y (d) disclaimer pequeno en cursiva que dice
//     "Conteo por cromosoma × subgenoma exacto · posicion intracromosomica
//     esquematica" — porque sabemos a que cromosoma y subgenoma pertenece cada
//     locus (gene_id Fxa{chr}{subg}gXXXXX), pero NO la coordenada genomica
//     dentro del cromosoma.
//
// 4 acuaporinomas referencia → tblastn → BEDtools/Exonerate → 129 candidatas
// distribuidas en 28 cromosomas del genoma alo-octoploide (subgenomas A/B/C/D
// por cromosomas 1-7). F. x ananassa es 2n=8x=56: 56 cromosomas en el set
// diploide, 28 distintos en el set haploide (7 cromosomas base × 4
// subgenomas, cada uno en pareja). Aqui representamos las 28 haploides.
export function init(slide) {
  if (!window.gsap) {
    const script = document.createElement('script');
    script.src = 'vendor/gsap.min.js';
    document.head.appendChild(script);
  }

  ensurePipelineStyles();

  const container = document.createElement('div');
  container.style.cssText = 'position:absolute;inset:0;width:100%;height:100%;z-index:0;display:flex;align-items:center;justify-content:center;';
  slide.insertBefore(container, slide.firstChild);

  const svg = createPipelineSVG();
  container.appendChild(svg);

  // Carga asincrona de los datos reales de distribucion; mientras tanto, los
  // huecos del grid quedan vacios. Si la carga falla, recurrimos a un fallback
  // con los totales por subgenoma del propio JSON hardcodeado aqui.
  let disposed = false;
  loadLociDistribution()
    .then(data => {
      if (disposed) return;
      buildChromosomeGrid(svg, data);
      buildExonerateFan(svg);   // abanico de salida ahora que el grid existe
      setTimeout(() => animatePipeline(svg), 200);
    })
    .catch(err => {
      console.warn('[genome-pipeline] fetch de distribucion fallido, fallback', err);
      if (disposed) return;
      buildChromosomeGrid(svg, FALLBACK_DATA);
      buildExonerateFan(svg);
      setTimeout(() => animatePipeline(svg), 200);
    });

  return function dispose() {
    disposed = true;
    container.remove();
  };
}

// ---------------------------------------------------------------------------
// Inyeccion de CSS reutilizable
// ---------------------------------------------------------------------------
function ensurePipelineStyles() {
  if (document.getElementById('genome-pipeline-styles')) return;
  const style = document.createElement('style');
  style.id = 'genome-pipeline-styles';
  style.textContent = `
    .gp-backdrop { fill: #1a2138; opacity: 0.55; filter: url(#gp-blur); }
    .gp-card     { fill: #1a2138; stroke: rgba(78,197,224,0.45); stroke-width: 1; }
    .gp-card-glow{ filter: url(#gp-cardglow); }
    .gp-mono     { font-family: 'IBM Plex Mono', monospace; }
    .gp-serif    { font-family: 'Cormorant Garamond', Georgia, serif; }
  `;
  document.head.appendChild(style);
}

// ---------------------------------------------------------------------------
// Datos de la escena
// ---------------------------------------------------------------------------
// Queries del tblastn: 4 acuaporinomas curados + conjunto Rosaceae (NCBI RefSeq)
const SPECIES = [
  { name: 'Arabidopsis thaliana',  short: 'A. thaliana',      glyph: '\u{1F33F}', count: 35,  y: 235 },
  { name: 'Oryza sativa',          short: 'O. sativa',        glyph: '\u{1F33E}', count: 30,  y: 355 },
  { name: 'Malus domestica',       short: 'M. domestica',     glyph: '\u{1F34F}', count: 41,  y: 475 },
  { name: 'Hevea brasiliensis',    short: 'H. brasiliensis',  glyph: '\u{1F332}', count: 51,  y: 595 },
  { name: 'Rosaceae · NCBI RefSeq', short: 'Rosaceae RefSeq', glyph: '\u{1F33A}', count: 419, y: 715 }
];

// Tres etapas en fila horizontal limpia, todas al mismo y. Cards iguales.
const STAGES = [
  { label: 'tblastn',   sublabel: 'e < 1e-05',           accent: '#4ec5e0', count: '53.396', countLabel: 'hits' },
  { label: 'BEDtools',  sublabel: 'merge 3kb · slop 1kb', accent: '#7dd3a8', count: '3.168',  countLabel: 'locus' },
  { label: 'Exonerate', sublabel: 'protein2genome',      accent: '#f0b65a', count: '4.984',  countLabel: 'modelos' }
];

const STAGE_Y    = 380;   // centro vertical de las cartas
const STAGE_X0   = 560;   // x del centro de la primera carta
const STAGE_GAP  = 180;   // separacion entre centros
const STAGE_W    = 160;
const STAGE_H    = 130;

// Strip inferior (count-up). Ahora termina en "129 candidatas" — las 121
// funcionales son territorio del Bloque II, no de esta escena.
// 4 genomas (Arabidopsis, Oryza, Malus, Hevea) + 1 conjunto Rosaceae (419) = 5 fuentes de query
const STATS = [
  { value: 5,      label: 'FUENTES',     color: '#8b96b3', decimals: 0 },
  { value: 53396,  label: 'HITS',        color: '#4ec5e0', decimals: 0 },
  { value: 3168,   label: 'LOCUS',       color: '#7dd3a8', decimals: 0 },
  { value: 4984,   label: 'MODELOS',     color: '#f0b65a', decimals: 0 },
  { value: 129,    label: 'CANDIDATAS',  color: '#ff5577', decimals: 0 }
];

// Grid de cromosomas (columna por subgenoma, fila por cromosoma)
const CHROM_COLS    = ['A', 'B', 'C', 'D'];
const CHROM_ROWS    = [1, 2, 3, 4, 5, 6, 7];
const GRID_X0       = 1100;   // borde izquierdo del area del grid
const GRID_X1       = 1500;   // borde derecho
const GRID_Y0       = 150;    // borde superior (header esta encima)
const GRID_Y1       = 760;    // borde inferior del grid (bajamos un poco para
                              // hacer hueco al bloque de legenda + disclaimer)
const CHROM_W       = 38;     // ancho de cada mini-cromosoma
const CHROM_H       = 60;     // alto de cada mini-cromosoma

// Fallback con los totales del JSON real (por si fetch falla)
const FALLBACK_DATA = {
  total: 129,
  matriz: {
    '1': { A: 2, B: 1, C: 2, D: 2, total: 7 },
    '2': { A: 4, B: 4, C: 1, D: 4, total: 13 },
    '3': { A: 5, B: 4, C: 4, D: 4, total: 17 },
    '4': { A: 2, B: 2, C: 2, D: 2, total: 8 },
    '5': { A: 5, B: 5, C: 5, D: 3, total: 18 },
    '6': { A: 15, B: 10, C: 9, D: 9, total: 43 },
    '7': { A: 6, B: 6, C: 6, D: 5, total: 23 }
  },
  totales_por_subgenoma: { A: 39, B: 32, C: 29, D: 29 }
};

async function loadLociDistribution() {
  const res = await fetch('assets/loci_129_distribucion.json', { cache: 'force-cache' });
  if (!res.ok) throw new Error('HTTP ' + res.status);
  return await res.json();
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------
const SVG_NS = 'http://www.w3.org/2000/svg';
const TXT = (parent, attrs, text) => {
  const t = document.createElementNS(SVG_NS, 'text');
  for (const k in attrs) {
    if (k === '_data' && attrs._data) {
      for (const dk in attrs._data) t.dataset[dk] = attrs._data[dk];
    } else {
      t.setAttribute(k, attrs[k]);
    }
  }
  if (text != null) t.textContent = text;
  parent.appendChild(t);
  return t;
};
const EL = (tag, attrs, parent) => {
  const e = document.createElementNS(SVG_NS, tag);
  if (attrs) for (const k in attrs) {
    if (k === '_data' && attrs._data) {
      for (const dk in attrs._data) e.dataset[dk] = attrs._data[dk];
    } else {
      e.setAttribute(k, attrs[k]);
    }
  }
  if (parent) parent.appendChild(e);
  return e;
};

// PRNG deterministico para colocar los puntos en cromosomas con apariencia
// "natural" pero estable entre re-renders.
function mulberry32(seed) {
  return function() {
    let t = (seed += 0x6D2B79F5) | 0;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

// ---------------------------------------------------------------------------
// SVG raiz + secciones
// ---------------------------------------------------------------------------
function createPipelineSVG() {
  const svg = document.createElementNS(SVG_NS, 'svg');
  svg.setAttribute('viewBox', '0 0 1600 900');
  svg.setAttribute('width', '96%');
  svg.setAttribute('height', '96%');

  // --- Defs ---------------------------------------------------------------
  const defs = EL('defs', null, svg);
  const fBlur = EL('filter', { id: 'gp-blur', x: '-10%', y: '-10%', width: '120%', height: '120%' }, defs);
  EL('feGaussianBlur', { stdDeviation: '6' }, fBlur);
  const fCardGlow = EL('filter', { id: 'gp-cardglow', x: '-30%', y: '-30%', width: '160%', height: '160%' }, defs);
  EL('feGaussianBlur', { stdDeviation: '4', result: 'b' }, fCardGlow);
  const merge = EL('feMerge', null, fCardGlow);
  EL('feMergeNode', { in: 'b' }, merge);
  EL('feMergeNode', { in: 'SourceGraphic' }, merge);

  // Gradiente vertical para cromosomas (bandeo G-style claro)
  const chromGrad = EL('linearGradient', { id: 'gp-chrom-grad', x1: '0', y1: '0', x2: '0', y2: '1' }, defs);
  const bands = [
    [0.00, 'rgba(232,240,255,0.04)'],
    [0.18, 'rgba(232,240,255,0.22)'],
    [0.30, 'rgba(232,240,255,0.08)'],
    [0.45, 'rgba(232,240,255,0.04)'], // centromero
    [0.55, 'rgba(232,240,255,0.04)'],
    [0.70, 'rgba(232,240,255,0.20)'],
    [0.85, 'rgba(232,240,255,0.08)'],
    [1.00, 'rgba(232,240,255,0.05)']
  ];
  bands.forEach(([o, c]) => EL('stop', { offset: o, 'stop-color': c }, chromGrad));

  // Radial para halo de los locus dots
  const haloGrad = EL('radialGradient', { id: 'gp-halo' }, defs);
  EL('stop', { offset: '0%',   'stop-color': '#ff5577', 'stop-opacity': '0.85' }, haloGrad);
  EL('stop', { offset: '70%',  'stop-color': '#ff5577', 'stop-opacity': '0.15' }, haloGrad);
  EL('stop', { offset: '100%', 'stop-color': '#ff5577', 'stop-opacity': '0' }, haloGrad);

  // ====================== Columna izquierda · cartas de especies ============
  buildSpeciesCards(svg);

  // ====================== Centro · etapas en fila horizontal limpia =========
  buildStages(svg);

  // ====================== Strip inferior ===================================
  buildStatStrip(svg);

  // El grid de cromosomas y el abanico de salida se construyen despues, cuando
  // la JSON haya cargado (buildChromosomeGrid + buildExonerateFan).

  return svg;
}

// ---------------------------------------------------------------------------
// Cartas de especies (izquierda)
// ---------------------------------------------------------------------------
function buildSpeciesCards(svg) {
  const cardW = 270, cardH = 90;
  const cardX = 60;

  SPECIES.forEach((sp, i) => {
    const g = EL('g', {
      opacity: '0',
      _data: { role: 'species-card', spIndex: i }
    }, svg);

    EL('rect', {
      x: cardX, y: sp.y - cardH / 2,
      width: cardW, height: cardH,
      rx: 8, ry: 8,
      class: 'gp-card gp-card-glow'
    }, g);

    EL('rect', {
      x: cardX + 1, y: sp.y - cardH / 2 + 1,
      width: cardW - 2, height: cardH - 2,
      rx: 7, ry: 7, fill: 'none',
      stroke: 'rgba(78,197,224,0.18)', 'stroke-width': '0.5'
    }, g);

    TXT(g, {
      x: cardX + 28, y: sp.y + 12,
      'font-size': 34,
      'text-anchor': 'middle'
    }, sp.glyph);

    TXT(g, {
      x: cardX + 62, y: sp.y - 6,
      class: 'gp-mono', 'font-size': 17, 'font-style': 'italic',
      fill: '#e8f0ff'
    }, sp.short);

    TXT(g, {
      x: cardX + 62, y: sp.y + 20,
      class: 'gp-mono', 'font-size': 13,
      fill: '#8b96b3', 'letter-spacing': '1.5'
    }, `${sp.count} acuaporinas`);

    // Punto de salida derecho
    EL('circle', {
      cx: cardX + cardW + 6, cy: sp.y,
      r: 4, fill: '#4ec5e0', opacity: '0',
      _data: { role: 'source-dot', spIndex: i }
    }, svg);
  });
}

// ---------------------------------------------------------------------------
// Centro: 3 cartas de etapas en una fila horizontal LIMPIA
// ---------------------------------------------------------------------------
function buildStages(svg) {
  const SRC_X = 60 + 270 + 6;   // borde derecho de las cartas + 6
  const stage1CX = STAGE_X0;
  const stage1CY = STAGE_Y;

  // Paths de cada especie hacia la primera etapa (centro de carta tblastn)
  SPECIES.forEach((sp, i) => {
    const d = `M ${SRC_X} ${sp.y} C ${SRC_X + 100} ${sp.y}, ${stage1CX - STAGE_W / 2 - 80} ${stage1CY}, ${stage1CX - STAGE_W / 2} ${stage1CY}`;
    EL('path', {
      d, fill: 'none',
      stroke: '#4ec5e0', 'stroke-opacity': '0.28', 'stroke-width': '1.3',
      'stroke-dasharray': '900', 'stroke-dashoffset': '900',
      id: `gp-flow-src-${i}`,
      _data: { role: 'flow-src', spIndex: i }
    }, svg);

    // Particulas SMIL en el camino especie → tblastn. Empiezan despues de que
    // se dibujen las paths (timing: cartas 0-0.6, paths 0.4-1.2, stages 1.0-1.8;
    // las particulas se hacen visibles a ~1.5 s).
    [0, 0.55].forEach((delayFrac, p) => {
      const dot = EL('circle', {
        r: 2.5, fill: '#4ec5e0', opacity: '0',
        _data: { role: 'flow-particle', src: i }
      }, svg);
      const anim = EL('animateMotion', {
        dur: '3.5s', repeatCount: 'indefinite',
        begin: `${1.6 + delayFrac * 1.75}s`,
        rotate: 'auto'
      }, dot);
      const mp = EL('mpath', {}, anim);
      mp.setAttributeNS('http://www.w3.org/1999/xlink', 'href', `#gp-flow-src-${i}`);
      EL('animate', {
        attributeName: 'opacity',
        from: '0', to: '0.85',
        dur: '0.5s', begin: `${1.6 + delayFrac * 1.75}s`,
        fill: 'freeze'
      }, dot);
    });
  });

  // 3 cartas + tramos entre cartas
  STAGES.forEach((st, i) => {
    const cx = STAGE_X0 + i * STAGE_GAP;
    const left = cx - STAGE_W / 2;
    const top = STAGE_Y - STAGE_H / 2;

    // Carta (rect)
    const cardG = EL('g', {
      opacity: '0',
      _data: { role: 'stage-card', stageIndex: i }
    }, svg);
    EL('rect', {
      x: left, y: top,
      width: STAGE_W, height: STAGE_H,
      rx: 10, ry: 10,
      fill: '#1a2138', stroke: st.accent, 'stroke-width': '1.2',
      'stroke-opacity': '0.65'
    }, cardG);
    EL('rect', {
      x: left + 1, y: top + 1,
      width: STAGE_W - 2, height: STAGE_H - 2,
      rx: 9, ry: 9, fill: 'none',
      stroke: st.accent, 'stroke-opacity': '0.18', 'stroke-width': '0.5'
    }, cardG);

    // Icono (circulo) a la izquierda del titulo, dentro de la carta
    EL('circle', {
      cx: cx - 56, cy: top + 28,
      r: 6, fill: 'none',
      stroke: st.accent, 'stroke-width': '1.4'
    }, cardG);
    EL('circle', {
      cx: cx - 56, cy: top + 28,
      r: 2.4, fill: st.accent
    }, cardG);

    // Titulo (etiqueta de la herramienta)
    TXT(cardG, {
      x: cx - 40, y: top + 32,
      class: 'gp-mono', 'font-size': 20, 'font-weight': '600',
      fill: st.accent, 'text-anchor': 'start'
    }, st.label);

    // Sublabel (parametros) — letter-spacing reducido para que entre el sublabel
    // mas largo (BEDtools "merge 3kb · slop 1kb") sin desbordar la carta de 160px.
    TXT(cardG, {
      x: cx, y: top + 56,
      class: 'gp-mono', 'font-size': 11,
      'letter-spacing': '0.4', fill: '#8b96b3',
      'text-anchor': 'middle'
    }, st.sublabel);

    // Separador horizontal sutil dentro de la carta
    EL('line', {
      x1: left + 16, x2: left + STAGE_W - 16,
      y1: top + 72, y2: top + 72,
      stroke: 'rgba(232,240,255,0.12)', 'stroke-width': '1'
    }, cardG);

    // Conteo grande (inicial = 0; count-up GSAP lleva al valor final al entrar
    // la carta). snapToFinal restaura el valor target si GSAP falla.
    TXT(cardG, {
      x: cx, y: top + 100,
      class: 'gp-serif', 'font-size': 30,
      fill: st.accent, 'text-anchor': 'middle',
      _data: { role: 'stage-count-text', stageIndex: i, target: st.count }
    }, '0');

    // Etiqueta de la unidad (hits/locus/modelos)
    TXT(cardG, {
      x: cx, y: top + 118,
      class: 'gp-mono', 'font-size': 10,
      'letter-spacing': '2', fill: '#8b96b3',
      'text-anchor': 'middle'
    }, st.countLabel.toUpperCase());

    // Conector entre cartas (excepto despues de la ultima)
    if (i < STAGES.length - 1) {
      const ax = cx + STAGE_W / 2;
      const bx = cx + STAGE_GAP - STAGE_W / 2;
      const d = `M ${ax} ${STAGE_Y} L ${bx} ${STAGE_Y}`;
      EL('path', {
        d, fill: 'none',
        stroke: 'rgba(232,240,255,0.20)', 'stroke-width': '1',
        opacity: '0',
        id: `gp-flow-stage-${i}`,
        _data: { role: 'flow-stage', stageIndex: i }
      }, svg);

      // Pequena flechita en el medio (entre cartas)
      const midX = (ax + bx) / 2;
      EL('path', {
        d: `M ${midX - 4} ${STAGE_Y - 4} L ${midX + 4} ${STAGE_Y} L ${midX - 4} ${STAGE_Y + 4}`,
        fill: 'none',
        stroke: st.accent, 'stroke-opacity': '0.55',
        'stroke-width': '1.2',
        'stroke-linejoin': 'round',
        opacity: '0',
        _data: { role: 'stage-arrow', stageIndex: i }
      }, svg);

      // Particulas que viajan entre las cartas — arrancan tras stage count-up
      const color = i === 0 ? '#7dd3a8' : '#f0b65a';
      [0, 0.5].forEach((df, p) => {
        const dot = EL('circle', {
          r: 2.5, fill: color, opacity: '0',
          _data: { role: 'stage-particle', stage: i }
        }, svg);
        const anim = EL('animateMotion', {
          dur: '2.2s', repeatCount: 'indefinite',
          begin: `${2.4 + df * 1.1}s`,
          rotate: 'auto'
        }, dot);
        const mp = EL('mpath', {}, anim);
        mp.setAttributeNS('http://www.w3.org/1999/xlink', 'href', `#gp-flow-stage-${i}`);
        EL('animate', {
          attributeName: 'opacity',
          from: '0', to: '0.9',
          dur: '0.4s', begin: `${2.4 + df * 1.1}s`,
          fill: 'freeze'
        }, dot);
      });
    }
  });
}

// ---------------------------------------------------------------------------
// Flecha Exonerate → grid de cromosomas
// ---------------------------------------------------------------------------
// Flecha simple, recta-curvada, que apunta del borde derecho de la carta
// Exonerate al borde izquierdo del grid (centrado verticalmente). Reemplaza
// al abanico de V4 que no quedaba claro visualmente.
function buildExonerateFan(svg) {
  const lastCX = STAGE_X0 + (STAGES.length - 1) * STAGE_GAP;
  const startX = lastCX + STAGE_W / 2 + 2;
  const startY = STAGE_Y;

  // Borde izquierdo del area util del grid, centro vertical
  const LABEL_LEFT_W = 50;
  const HEADER_TOP_H = 28;
  const usableX0 = GRID_X0 + LABEL_LEFT_W;
  const gridCenterY = GRID_Y0 + HEADER_TOP_H + (GRID_Y1 - GRID_Y0 - HEADER_TOP_H) / 2;
  const endX = usableX0 - 18;
  const endY = gridCenterY;

  // Grupo para controlar opacidad/animacion
  const arrowG = EL('g', {
    opacity: '0',
    _data: { role: 'exonerate-fan' }
  }, svg);

  // Trazo principal: bezier sutil
  const cp1x = startX + 80;
  const cp1y = startY;
  const cp2x = endX - 80;
  const cp2y = endY;
  const d = `M ${startX} ${startY} C ${cp1x} ${cp1y}, ${cp2x} ${cp2y}, ${endX} ${endY}`;

  EL('path', {
    d, fill: 'none',
    stroke: '#ff5577',
    'stroke-opacity': '0.85',
    'stroke-width': '2',
    'stroke-linecap': 'round',
    id: 'gp-fan-line-0',
    _data: { role: 'fan-line', idx: 0 }
  }, arrowG);

  // Punta de flecha (triangle) en endX,endY apuntando a la derecha
  const headSize = 9;
  EL('polygon', {
    points: `${endX},${endY - headSize / 2} ${endX + headSize},${endY} ${endX},${endY + headSize / 2}`,
    fill: '#ff5577',
    'fill-opacity': '0.95',
    _data: { role: 'fan-dot', idx: 0 }
  }, arrowG);

  // Particula viajando sobre la linea para sugerir transferencia de datos
  const part = EL('circle', {
    r: 3, fill: '#ff5577', opacity: '0',
    _data: { role: 'fan-particle', idx: 0 }
  }, arrowG);
  const anim = EL('animateMotion', {
    dur: '2.2s', repeatCount: 'indefinite',
    begin: '4.4s',
    rotate: 'auto'
  }, part);
  const mp = EL('mpath', {}, anim);
  mp.setAttributeNS('http://www.w3.org/1999/xlink', 'href', '#gp-fan-line-0');
  EL('animate', {
    attributeName: 'opacity',
    from: '0', to: '0.95',
    dur: '0.4s', begin: '4.4s',
    fill: 'freeze'
  }, part);

  // Etiqueta encima de la flecha: "129 candidatas"
  const labelX = (startX + endX) / 2;
  const labelY = (startY + endY) / 2 - 14;
  TXT(arrowG, {
    x: labelX, y: labelY,
    class: 'gp-mono', 'font-size': 12,
    'letter-spacing': '1.5',
    'text-anchor': 'middle',
    fill: '#ff5577', opacity: '0',
    _data: { role: 'fan-label' }
  }, '129 candidatas');
}

// ---------------------------------------------------------------------------
// Grid de cromosomas (4 columnas A/B/C/D x 7 filas Chr 1..7) = 28 capsulas
// ---------------------------------------------------------------------------
function buildChromosomeGrid(svg, data) {
  const gridG = EL('g', { _data: { role: 'chrom-grid' } }, svg);

  const totalCols = CHROM_COLS.length;
  const totalRows = CHROM_ROWS.length;

  const LABEL_LEFT_W = 50;
  const HEADER_TOP_H = 28;

  const usableX0 = GRID_X0 + LABEL_LEFT_W;
  const usableX1 = GRID_X1;
  const usableY0 = GRID_Y0 + HEADER_TOP_H;
  const usableY1 = GRID_Y1 - 24;   // dejamos espacio para la leyenda + disclaimer abajo

  const colW = (usableX1 - usableX0) / totalCols;
  const rowH = (usableY1 - usableY0) / totalRows;

  // Headers de columnas (subgenomas) — solo la letra A/B/C/D
  // (Noé verbaliza "subgenoma" en la defensa). Opacity inicial 0 → animacion.
  CHROM_COLS.forEach((sg, c) => {
    const cx = usableX0 + colW * (c + 0.5);
    TXT(gridG, {
      x: cx, y: GRID_Y0 + 22,
      'text-anchor': 'middle',
      class: 'gp-mono', 'font-size': 22, 'letter-spacing': '1',
      'font-weight': 'bold',
      fill: '#e8f0ff', opacity: '0',
      _data: { role: 'grid-col-header', col: c }
    }, sg);
  });

  // Headers de filas (cromosomas) — opacity 0 hasta el GSAP del scaffold (2,4 s)
  CHROM_ROWS.forEach((chr, r) => {
    const cy = usableY0 + rowH * (r + 0.5);
    TXT(gridG, {
      x: GRID_X0 + 14, y: cy + 6,
      class: 'gp-mono', 'font-size': 17, 'letter-spacing': '1.5',
      'font-weight': 'bold',
      fill: '#c7cee2', opacity: '0',
      _data: { role: 'grid-row-label', row: r }
    }, `Chr ${chr}`);
  });

  // Cada celda = mini-cromosoma con N puntos
  CHROM_ROWS.forEach((chr, r) => {
    const rowData = data.matriz[String(chr)] || { A: 0, B: 0, C: 0, D: 0 };
    CHROM_COLS.forEach((sg, c) => {
      const cx = usableX0 + colW * (c + 0.5);
      const cy = usableY0 + rowH * (r + 0.5);
      const count = rowData[sg] || 0;
      drawMiniChromosome(gridG, cx, cy, count, r, c);
    });
  });

  // ----------- Leyenda + disclaimer DEBAJO del grid ------------------------
  const subgTotals = data.totales_por_subgenoma || { A: 0, B: 0, C: 0, D: 0 };
  const chrSums = {};
  CHROM_ROWS.forEach(chr => { chrSums[chr] = (data.matriz[String(chr)] || {}).total || 0; });
  const chr6 = chrSums[6];
  const total = data.total;

  const sumG = EL('g', {
    opacity: '0',
    _data: { role: 'chrom-sum' }
  }, gridG);

  const sumCenterX = (GRID_X0 + GRID_X1) / 2;
  const yBase = GRID_Y1 + 10;

  // 1 · titular único: total + cariotipo
  TXT(sumG, {
    x: sumCenterX, y: yBase,
    'text-anchor': 'middle',
    class: 'gp-serif', 'font-size': 24, fill: '#ff5577'
  }, `${total} candidatas en 28 cromosomas (2n=8x=56)`);

  // 2 · línea de contexto: Chr 6 dominante
  TXT(sumG, {
    x: sumCenterX, y: yBase + 22,
    'text-anchor': 'middle',
    class: 'gp-mono', 'font-size': 11, 'letter-spacing': '1.5',
    fill: '#8b96b3'
  }, `Chr 6: ${chr6} loci (${Math.round(chr6 / total * 100)}%) · posición intracromosómica esquemática`);
}

// ---------------------------------------------------------------------------
// Mini-cromosoma con puntos
// ---------------------------------------------------------------------------
function drawMiniChromosome(parent, cx, cy, count, row, col) {
  // Capsula vertical con centromero
  const halfW = CHROM_W / 2;
  const halfH = CHROM_H / 2;
  const yTop = cy - halfH;
  const yBot = cy + halfH;
  const r = halfW * 0.85;
  const pinchHalfW = halfW - 4;
  const yMid = cy;

  const path = [
    `M ${cx - halfW} ${yTop + r}`,
    `Q ${cx - halfW} ${yTop} ${cx} ${yTop}`,
    `Q ${cx + halfW} ${yTop} ${cx + halfW} ${yTop + r}`,
    `L ${cx + halfW} ${yMid - 6}`,
    `Q ${cx + halfW - 3} ${yMid} ${cx + pinchHalfW} ${yMid}`,
    `Q ${cx + halfW - 3} ${yMid} ${cx + halfW} ${yMid + 6}`,
    `L ${cx + halfW} ${yBot - r}`,
    `Q ${cx + halfW} ${yBot} ${cx} ${yBot}`,
    `Q ${cx - halfW} ${yBot} ${cx - halfW} ${yBot - r}`,
    `L ${cx - halfW} ${yMid + 6}`,
    `Q ${cx - halfW + 3} ${yMid} ${cx - pinchHalfW} ${yMid}`,
    `Q ${cx - halfW + 3} ${yMid} ${cx - halfW} ${yMid - 6}`,
    `Z`
  ].join(' ');

  EL('path', {
    d: path,
    fill: 'url(#gp-chrom-grad)',
    stroke: 'rgba(232,240,255,0.32)',
    'stroke-width': '0.8',
    opacity: '0',
    _data: { role: 'chrom-body', row, col }
  }, parent);

  // Centromero (linea sutil)
  EL('line', {
    x1: cx - pinchHalfW - 1, y1: yMid,
    x2: cx + pinchHalfW + 1, y2: yMid,
    stroke: 'rgba(232,240,255,0.45)',
    'stroke-width': '0.7',
    opacity: '0',
    _data: { role: 'chrom-cent', row, col }
  }, parent);

  if (count <= 0) return;

  // Distribuir 'count' puntos a lo largo del cromosoma con PRNG estable.
  const rng = mulberry32(row * 100 + col * 7 + count * 31);
  for (let i = 0; i < count; i++) {
    const tFrac = 0.05 + 0.90 * ((i + rng() * 0.6 + 0.2) / Math.max(count, 1));
    let yDot = yTop + tFrac * CHROM_H;
    if (Math.abs(yDot - yMid) < 8) {
      yDot += (yDot < yMid ? -10 : 10);
    }
    const offX = (rng() - 0.5) * (CHROM_W * 0.45);

    EL('circle', {
      cx: cx + offX, cy: yDot,
      r: 5, fill: 'url(#gp-halo)',
      opacity: '0',
      _data: { role: 'locus-halo', row, col }
    }, parent);

    EL('circle', {
      cx: cx + offX, cy: yDot,
      r: 1.8, fill: '#ff5577', opacity: '0',
      _data: { role: 'locus-dot', row, col }
    }, parent);
  }
}

// ---------------------------------------------------------------------------
// Strip inferior
// ---------------------------------------------------------------------------
function buildStatStrip(svg) {
  const stripY = 855;
  const x0 = 130, x1 = 1480;
  const sep = (x1 - x0) / (STATS.length - 1);

  EL('line', {
    x1: x0, x2: x1, y1: stripY, y2: stripY,
    stroke: 'rgba(232,240,255,0.10)',
    'stroke-width': '1',
    opacity: '0',
    _data: { role: 'strip-line' }
  }, svg);

  STATS.forEach((st, i) => {
    const cx = x0 + sep * i;

    EL('line', {
      x1: cx, x2: cx, y1: stripY - 5, y2: stripY + 5,
      stroke: 'rgba(232,240,255,0.25)',
      'stroke-width': '1',
      opacity: '0',
      _data: { role: 'strip-tick', statIndex: i }
    }, svg);

    TXT(svg, {
      x: cx, y: stripY - 12,
      'text-anchor': 'middle',
      class: 'gp-serif', 'font-size': 26, fill: st.color, opacity: '0',
      _data: { role: 'stat-value', statIndex: i, target: st.value, decimals: st.decimals }
    }, '0');

    TXT(svg, {
      x: cx, y: stripY + 22,
      'text-anchor': 'middle',
      class: 'gp-mono', 'font-size': 10,
      'letter-spacing': '2.5',
      fill: '#8b96b3', opacity: '0',
      _data: { role: 'stat-label', statIndex: i }
    }, st.label);

    if (i < STATS.length - 1) {
      const midX = cx + sep / 2;
      TXT(svg, {
        x: midX, y: stripY + 4,
        'text-anchor': 'middle',
        class: 'gp-mono', 'font-size': 14,
        fill: 'rgba(232,240,255,0.18)', opacity: '0',
        _data: { role: 'strip-arrow', statIndex: i }
      }, '→');
    }
  });
}

// ---------------------------------------------------------------------------
// Animacion (timeline GSAP) — V4: continuidad y orden estricto
// ---------------------------------------------------------------------------
//
// Plan:
//   0.0 – 0.6  cartas de especies entran (fade + slide desde la izquierda)
//   0.4 – 1.2  flow lines especies → tblastn se dibujan (despues de que las
//              cartas hayan empezado a aterrizar)
//   1.0 – 1.8  cartas de etapas tblastn/BEDtools/Exonerate aparecen
//   1.6 – 2.4  count-up de sub-conteos en cada carta de etapa
//   1.8 – 2.2  conectores rectos entre etapas + flechitas
//   2.4 – ∞    particulas entre etapas (loop continuo)
//   2.4 – 3.4  scaffold del grid: column headers, row labels, capsulas
//   3.0 – 4.2  cascada chromosome-by-chromosome: cuerpo + centromero + halos +
//              puntos (chr 1 → 7, ~0.18 s entre filas)
//   4.0 – 4.7  abanico Exonerate → grid se dibuja + etiqueta inline
//   4.7 – 5.3  bloque de leyenda/disclaimer debajo del grid
//   5.0 – 6.5  strip inferior aparece y los count-up resuelven
// ---------------------------------------------------------------------------
function animatePipeline(svg) {
  const gsap = window.gsap;
  if (!gsap) {
    console.warn('[genome-pipeline] GSAP no cargado, mostrando estado final');
    snapToFinal(svg);
    return;
  }

  const tl = gsap.timeline();

  // 0.0–0.6 · cartas de especies (con stagger pequeno)
  const cards = svg.querySelectorAll('[data-role="species-card"]');
  tl.fromTo(cards,
    { opacity: 0, x: -40, scale: 0.95 },
    { opacity: 1, x: 0, scale: 1, duration: 0.55, stagger: 0.10, ease: 'power3.out' },
    0
  );

  // Source dots aparecen tras la carta
  const sourceDots = svg.querySelectorAll('[data-role="source-dot"]');
  sourceDots.forEach(d => {
    d.style.transformOrigin = `${d.getAttribute('cx')}px ${d.getAttribute('cy')}px`;
    d.style.transformBox = 'fill-box';
  });
  tl.fromTo(sourceDots,
    { scale: 0, opacity: 0 },
    { scale: 1, opacity: 1, duration: 0.35, stagger: 0.08, ease: 'back.out(2)' },
    0.4
  );

  // 0.4–1.2 · flow lines especie → tblastn (despues de que las cartas
  // arranquen su entrada; estan listas a 0.5)
  const flows = svg.querySelectorAll('[data-role="flow-src"]');
  tl.to(flows,
    { strokeDashoffset: 0, duration: 0.8, stagger: 0.06, ease: 'power2.inOut' },
    0.5
  );

  // 1.0–1.8 · cartas de etapas en fila
  const stageCards = svg.querySelectorAll('[data-role="stage-card"]');
  stageCards.forEach(c => { c.style.transformOrigin = 'center center'; c.style.transformBox = 'fill-box'; });
  tl.fromTo(stageCards,
    { opacity: 0, y: 16, scale: 0.95 },
    { opacity: 1, y: 0, scale: 1, duration: 0.55, stagger: 0.14, ease: 'power2.out' },
    1.0
  );

  // 1.1–1.9 · count-up de sub-conteos · arranca con el fade-in de cada carta
  // (carta i empieza a aparecer a 1.0 + i * 0.14; el count-up entra 0,1 s
  // despues para que el numero "0" no aparezca ya dibujado antes que la carta)
  const stageCountTexts = svg.querySelectorAll('[data-role="stage-count-text"]');
  stageCountTexts.forEach((el, i) => {
    const target = parseFloat(el.dataset.target.replace(/\./g, ''));
    const obj = { v: 0 };
    tl.to(obj, {
      v: target,
      duration: 0.75,
      ease: 'power2.out',
      onUpdate: () => {
        el.textContent = obj.v.toLocaleString('es-ES', { maximumFractionDigits: 0, useGrouping: 'always' });
      }
    }, 1.1 + i * 0.14);
  });

  // 1.8–2.2 · conectores rectos entre etapas + flechitas
  tl.to('[data-role="flow-stage"]', { opacity: 1, duration: 0.4, stagger: 0.10, ease: 'power2.out' }, 1.8);
  tl.to('[data-role="stage-arrow"]', { opacity: 1, duration: 0.35, stagger: 0.10, ease: 'power2.out' }, 1.9);

  // 2.4–3.4 · scaffold del grid (columnas, filas, cuerpos vacios)
  const colHeaders = svg.querySelectorAll('[data-role="grid-col-header"]');
  const rowLabels = svg.querySelectorAll('[data-role="grid-row-label"]');
  const chromBodies = svg.querySelectorAll('[data-role="chrom-body"]');

  tl.to(colHeaders, { opacity: 1, duration: 0.35, stagger: 0.06, ease: 'power2.out' }, 2.4);
  tl.to(rowLabels,  { opacity: 1, duration: 0.3, stagger: 0.04, ease: 'power2.out' }, 2.5);

  chromBodies.forEach(b => {
    b.style.transformOrigin = 'center center';
    b.style.transformBox = 'fill-box';
  });

  // 3.0–4.2 · cascada chr-by-chr: para cada fila r, todas sus 4 capsulas +
  // centromero + puntos + halos entran juntos. Stagger 0.16 entre filas.
  CHROM_ROWS.forEach((chr, r) => {
    const at = 3.0 + r * 0.16;
    const rowBodies = svg.querySelectorAll(`[data-role="chrom-body"][data-row="${r}"]`);
    const rowCents = svg.querySelectorAll(`[data-role="chrom-cent"][data-row="${r}"]`);
    const rowDots = svg.querySelectorAll(`[data-role="locus-dot"][data-row="${r}"]`);
    const rowHalos = svg.querySelectorAll(`[data-role="locus-halo"][data-row="${r}"]`);

    tl.fromTo(rowBodies,
      { opacity: 0, scaleY: 0.5 },
      { opacity: 1, scaleY: 1, duration: 0.4, ease: 'power2.out' },
      at
    );
    tl.to(rowCents, { opacity: 1, duration: 0.25 }, at + 0.22);

    rowDots.forEach(d => {
      d.style.transformOrigin = `${d.getAttribute('cx')}px ${d.getAttribute('cy')}px`;
      d.style.transformBox = 'fill-box';
    });
    rowHalos.forEach(h => {
      h.style.transformOrigin = `${h.getAttribute('cx')}px ${h.getAttribute('cy')}px`;
      h.style.transformBox = 'fill-box';
    });

    tl.fromTo(rowDots,
      { scale: 0, opacity: 0 },
      { scale: 1, opacity: 1, duration: 0.3, stagger: 0.008, ease: 'back.out(1.8)' },
      at + 0.28
    );
    tl.fromTo(rowHalos,
      { scale: 0.6, opacity: 0 },
      { scale: 2.0, opacity: 0, duration: 0.9, stagger: 0.008, ease: 'power2.out' },
      at + 0.32
    );
  });

  // 4.0–4.7 · abanico Exonerate → grid (lineas dibujan + dots + label)
  tl.to('[data-role="exonerate-fan"]', { opacity: 1, duration: 0.4 }, 4.0);
  tl.to('[data-role="fan-line"]',
    { strokeDashoffset: 0, duration: 0.7, stagger: 0.06, ease: 'power2.inOut' },
    4.05
  );
  tl.to('[data-role="fan-dot"]',
    { opacity: 0.9, duration: 0.3, stagger: 0.06, ease: 'back.out(2)' },
    4.5
  );
  tl.to('[data-role="fan-label"]',
    { opacity: 1, duration: 0.4 },
    4.5
  );

  // 4.7–5.3 · bloque de leyenda + disclaimer bajo el grid
  tl.fromTo('[data-role="chrom-sum"]',
    { opacity: 0, y: 8 },
    { opacity: 1, y: 0, duration: 0.55, ease: 'power2.out' },
    4.7
  );

  // 5.0–6.5 · strip inferior con count-up
  tl.to('[data-role="strip-line"]',  { opacity: 1, duration: 0.4 }, 5.0);
  tl.to('[data-role="strip-tick"]',  { opacity: 1, duration: 0.25, stagger: 0.06 }, 5.1);
  tl.to('[data-role="strip-arrow"]', { opacity: 1, duration: 0.25, stagger: 0.06 }, 5.15);
  tl.to('[data-role="stat-label"]',  { opacity: 1, duration: 0.35, stagger: 0.08 }, 5.3);

  const statValues = svg.querySelectorAll('[data-role="stat-value"]');
  statValues.forEach((el, i) => {
    const target = parseFloat(el.dataset.target);
    const decimals = parseInt(el.dataset.decimals) || 0;
    const obj = { v: 0 };
    tl.to(el, { opacity: 1, duration: 0.35 }, 5.3 + i * 0.10);
    tl.to(obj, {
      v: target,
      duration: 0.9,
      ease: 'power2.out',
      onUpdate: () => {
        el.textContent = obj.v.toLocaleString('es-ES', {
          minimumFractionDigits: decimals,
          maximumFractionDigits: decimals,
          useGrouping: 'always'
        });
      }
    }, 5.3 + i * 0.10);
  });
}

function snapToFinal(svg) {
  svg.querySelectorAll('[opacity="0"], [data-role]').forEach(el => {
    if (el.tagName === 'g' || el.dataset.role) el.setAttribute('opacity', 1);
  });
  svg.querySelectorAll('[data-role="flow-src"]').forEach(el => el.setAttribute('stroke-dashoffset', 0));
  svg.querySelectorAll('[data-role="fan-line"]').forEach(el => el.setAttribute('stroke-dashoffset', 0));
  svg.querySelectorAll('[data-role="stage-count-text"]').forEach(el => {
    el.textContent = el.dataset.target || el.textContent;
  });
  svg.querySelectorAll('[data-role="stat-value"]').forEach(el => {
    const target = parseFloat(el.dataset.target);
    const decimals = parseInt(el.dataset.decimals) || 0;
    el.textContent = target.toLocaleString('es-ES', {
      minimumFractionDigits: decimals, maximumFractionDigits: decimals,
      useGrouping: 'always'
    });
  });
}
