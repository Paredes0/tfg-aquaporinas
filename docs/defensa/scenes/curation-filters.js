// scenes/curation-filters.js — Bloque II separator (Curado).  v2 (2026-06-08)
//
// "La criba de cuatro filtros": un grupo de candidatas entra por la izquierda y
// atraviesa cuatro filtros independientes —Topología, Motivos, Parentesco y
// Físico-química—. Las que NO superan un filtro quedan CONGELADAS en rojo dentro
// de ese filtro (no han podido pasar); las viables salen por la derecha en verde.
// Cada filtro lleva un icono breve de qué comprueba.
//
// SVG + GSAP (con snap fallback). Animación pausada para que se siga con la vista;
// el estado final es autoexplicativo aunque no se vea el movimiento.

const SVG_NS = 'http://www.w3.org/2000/svg';

const GATES = [
  { n: '01', name: 'TOPOLOGÍA',      tool: 'DeepTMHMM', ask: '¿6 hélices?',        icon: 'helix', color: '#4ec5e0' },
  { n: '02', name: 'MOTIVOS',        tool: 'MEME',      ask: '¿dos motivos NPA?',  icon: 'npa',   color: '#7dd3a8' },
  { n: '03', name: 'PARENTESCO',     tool: 'IQ-TREE',   ask: '¿posición en el árbol?', icon: 'tree', color: '#f0b65a' },
  { n: '04', name: 'FÍSICO-QUÍMICA', tool: 'PCA',       ask: '¿química de su grupo?',  icon: 'pca',  color: '#c98bff' }
];

const GATE_X  = [470, 740, 1010, 1280];
const AXIS_Y  = 500;
const GATE_H  = 300;
const GATE_W  = 64;
const CAP_TOP = AXIS_Y - GATE_H / 2;          // 350
const CAP_BOT = AXIS_Y + GATE_H / 2;          // 650
const ENTRY_X = 215;
const EXIT_X  = 1430;
const N_DOTS  = 30;
// Índice de dot -> filtro en el que se queda atascado (2,2,1,1 = 6 descartes ⇒ 24 pasan).
// Las cifras reales (147 → 121) van en las etiquetas.
const FAIL_AT = { 3: 0, 9: 0, 14: 1, 19: 1, 23: 2, 27: 3 };

const C_CAND = '#aab9d6';   // candidata sin clasificar
const C_PASS = '#7dd3a8';   // viable
const C_FAIL = '#ff5a72';   // descartada (atascada en el filtro)

function mulberry32(seed) {
  return function () {
    let t = (seed += 0x6D2B79F5) | 0;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

const EL = (tag, attrs, parent) => {
  const e = document.createElementNS(SVG_NS, tag);
  if (attrs) for (const k in attrs) e.setAttribute(k, attrs[k]);
  if (parent) parent.appendChild(e);
  return e;
};
const TXT = (parent, attrs, text) => {
  const t = EL('text', attrs, parent);
  if (text != null) t.textContent = text;
  return t;
};

export function init(slide) {
  if (!window.gsap) {
    const s = document.createElement('script');
    s.src = 'vendor/gsap.min.js';
    document.head.appendChild(s);
  }

  const container = document.createElement('div');
  container.style.cssText =
    'position:absolute;inset:0;width:100%;height:100%;z-index:0;' +
    'display:flex;align-items:center;justify-content:center;background:#14182a;';
  slide.insertBefore(container, slide.firstChild);

  const svg = buildSVG();
  container.appendChild(svg);

  requestAnimationFrame(() => requestAnimationFrame(() => animate(svg)));

  return function dispose() { container.remove(); };
}

// --- iconos breves de cada filtro -------------------------------------------
function drawIcon(parent, type, cx, cy, color) {
  const g = EL('g', {}, parent);
  if (type === 'helix') {
    for (let j = 0; j < 6; j++) {
      const x = cx - 16 + j * 6.2;
      const h = (j % 2 === 0) ? 18 : 13;
      EL('rect', { x: x, y: cy - h / 2, width: 3, height: h, rx: 1.5, fill: color }, g);
    }
  } else if (type === 'npa') {
    EL('rect', { x: cx - 17, y: cy - 8, width: 14, height: 16, rx: 3, fill: 'none', stroke: color, 'stroke-width': 1.6 }, g);
    EL('rect', { x: cx + 3, y: cy - 8, width: 14, height: 16, rx: 3, fill: 'none', stroke: color, 'stroke-width': 1.6 }, g);
    EL('line', { x1: cx - 3, y1: cy, x2: cx + 3, y2: cy, stroke: color, 'stroke-width': 1.6 }, g);
    EL('circle', { cx: cx - 10, cy: cy, r: 1.8, fill: color }, g);
    EL('circle', { cx: cx + 10, cy: cy, r: 1.8, fill: color }, g);
  } else if (type === 'tree') {
    EL('path', {
      d: `M ${cx} ${cy + 9} L ${cx} ${cy - 2} M ${cx} ${cy - 2} L ${cx - 12} ${cy - 2} L ${cx - 12} ${cy - 9} M ${cx} ${cy - 2} L ${cx + 12} ${cy - 2} L ${cx + 12} ${cy - 9}`,
      fill: 'none', stroke: color, 'stroke-width': 1.6, 'stroke-linecap': 'round', 'stroke-linejoin': 'round'
    }, g);
    [cx - 12, cx + 12].forEach(tx => EL('circle', { cx: tx, cy: cy - 10, r: 2.2, fill: color }, g));
  } else if (type === 'pca') {
    EL('ellipse', { cx: cx, cy: cy, rx: 17, ry: 10, fill: 'none', stroke: color, 'stroke-width': 1.4, 'stroke-dasharray': '3 3' }, g);
    [[-7, -2], [3, -4], [6, 3], [-4, 4], [0, 1]].forEach(([dx, dy]) =>
      EL('circle', { cx: cx + dx, cy: cy + dy, r: 2, fill: color }, g));
  }
  return g;
}

function buildSVG() {
  const svg = document.createElementNS(SVG_NS, 'svg');
  svg.setAttribute('viewBox', '0 0 1600 900');
  svg.setAttribute('width', '95%');
  svg.setAttribute('height', '95%');

  const defs = EL('defs', null, svg);
  // glow para los puntos
  const glow = EL('radialGradient', { id: 'cf-glow' }, defs);
  EL('stop', { offset: '0%',   'stop-color': '#dfe9ff', 'stop-opacity': '0.85' }, glow);
  EL('stop', { offset: '100%', 'stop-color': '#dfe9ff', 'stop-opacity': '0' }, glow);
  // sombra suave para las cápsulas
  const sh = EL('filter', { id: 'cf-shadow', x: '-40%', y: '-20%', width: '180%', height: '140%' }, defs);
  EL('feDropShadow', { dx: 0, dy: 6, stdDeviation: 10, 'flood-color': '#000', 'flood-opacity': '0.35' }, sh);

  // guía horizontal sutil
  EL('line', {
    x1: ENTRY_X, x2: EXIT_X, y1: AXIS_Y, y2: AXIS_Y,
    stroke: 'rgba(232,240,255,0.07)', 'stroke-width': 1, 'data-role': 'axis', opacity: 0
  }, svg);

  // --- filtros ---
  GATES.forEach((gt, i) => {
    const gx = GATE_X[i];
    const grp = EL('g', { 'data-role': 'gate', opacity: 0 }, svg);

    // cápsula
    const capGrad = EL('linearGradient', { id: `cf-cap-${i}`, x1: 0, y1: 0, x2: 0, y2: 1 }, svg.querySelector('defs'));
    EL('stop', { offset: '0%',   'stop-color': gt.color, 'stop-opacity': 0.16 }, capGrad);
    EL('stop', { offset: '100%', 'stop-color': gt.color, 'stop-opacity': 0.04 }, capGrad);

    EL('rect', {
      x: gx - GATE_W / 2, y: CAP_TOP, width: GATE_W, height: GATE_H, rx: GATE_W / 2, ry: GATE_W / 2,
      fill: `url(#cf-cap-${i})`, stroke: gt.color, 'stroke-opacity': 0.8, 'stroke-width': 1.8,
      filter: 'url(#cf-shadow)'
    }, grp);
    // brillo interior superior
    EL('rect', {
      x: gx - GATE_W / 2 + 3, y: CAP_TOP + 3, width: GATE_W - 6, height: GATE_H - 6, rx: GATE_W / 2 - 3,
      fill: 'none', stroke: gt.color, 'stroke-opacity': 0.18, 'stroke-width': 1
    }, grp);

    // número arriba
    TXT(grp, { x: gx, y: CAP_TOP - 56, 'text-anchor': 'middle', 'font-family': 'IBM Plex Mono, monospace', 'font-size': 21, 'letter-spacing': 2, fill: gt.color }, gt.n);
    // icono
    drawIcon(grp, gt.icon, gx, CAP_TOP - 30, gt.color);

    // nombre / pregunta / herramienta (debajo)
    TXT(grp, { x: gx, y: CAP_BOT + 36, 'text-anchor': 'middle', 'font-family': 'IBM Plex Mono, monospace', 'font-size': 20, 'letter-spacing': 1.4, 'font-weight': 600, fill: gt.color }, gt.name);
    TXT(grp, { x: gx, y: CAP_BOT + 62, 'text-anchor': 'middle', 'font-family': 'IBM Plex Sans, sans-serif', 'font-size': 16, fill: '#c7cee2' }, gt.ask);
    TXT(grp, { x: gx, y: CAP_BOT + 84, 'text-anchor': 'middle', 'font-family': 'IBM Plex Mono, monospace', 'font-size': 13, 'letter-spacing': 1, fill: '#8b96b3' }, gt.tool);
  });

  // --- etiqueta de entrada ---
  const entryG = EL('g', { 'data-role': 'entry-label', opacity: 0 }, svg);
  TXT(entryG, { x: ENTRY_X, y: CAP_TOP - 44, 'text-anchor': 'middle', 'font-family': 'Cormorant Garamond, Georgia, serif', 'font-size': 60, fill: '#e8f0ff' }, '147');
  TXT(entryG, { x: ENTRY_X, y: CAP_TOP - 18, 'text-anchor': 'middle', 'font-family': 'IBM Plex Mono, monospace', 'font-size': 16, 'letter-spacing': 2, fill: '#8b96b3' }, 'CANDIDATAS');

  // --- etiqueta de salida ---
  const exitG = EL('g', { 'data-role': 'exit-label', opacity: 0 }, svg);
  TXT(exitG, { x: EXIT_X, y: CAP_TOP - 44, 'text-anchor': 'middle', 'font-family': 'Cormorant Garamond, Georgia, serif', 'font-size': 68, fill: C_PASS }, '121');
  TXT(exitG, { x: EXIT_X, y: CAP_TOP - 18, 'text-anchor': 'middle', 'font-family': 'IBM Plex Mono, monospace', 'font-size': 16, 'letter-spacing': 2, fill: '#8b96b3' }, 'VIABLES ✓');

  // --- candidatas (puntos) ---
  const rng = mulberry32(20260608);
  let surv = 0;
  for (let i = 0; i < N_DOTS; i++) {
    const ang = rng() * Math.PI * 2;
    const rad = 22 + rng() * 92;
    const sx = ENTRY_X + Math.cos(ang) * rad * 0.85;
    const sy = AXIS_Y + Math.sin(ang) * rad;

    const failGate = (i in FAIL_AT) ? FAIL_AT[i] : -1;
    let ex, ey;
    if (failGate >= 0) {
      // queda atascado dentro de la cápsula del filtro
      ex = GATE_X[failGate] + (rng() - 0.5) * (GATE_W - 26);
      ey = AXIS_Y + (rng() - 0.5) * (GATE_H - 70);
    } else {
      // columna ordenada a la salida (2 columnas)
      const col = surv % 2, row = Math.floor(surv / 2);
      ex = EXIT_X - 11 + col * 22;
      ey = (AXIS_Y - 138) + row * 24;
      surv++;
    }

    const halo = EL('circle', { cx: sx, cy: sy, r: 12, fill: 'url(#cf-glow)', opacity: 0, 'data-role': 'halo' }, svg);
    const dot = EL('circle', {
      cx: sx, cy: sy, r: 6, fill: C_CAND, opacity: 0, 'data-role': 'dot',
      'data-ex': ex, 'data-ey': ey, 'data-fail': failGate
    }, svg);
    dot._halo = halo;
  }

  return svg;
}

function animate(svg) {
  const gsap = window.gsap;
  const dots = Array.from(svg.querySelectorAll('[data-role="dot"]'));

  if (!gsap) {
    svg.querySelectorAll('[data-role="axis"],[data-role="gate"],[data-role="entry-label"],[data-role="exit-label"]')
      .forEach(el => el.setAttribute('opacity', 1));
    dots.forEach(d => {
      const fail = parseInt(d.dataset.fail);
      d.setAttribute('cx', d.dataset.ex);
      d.setAttribute('cy', d.dataset.ey);
      d.setAttribute('fill', fail >= 0 ? C_FAIL : C_PASS);
      d.setAttribute('opacity', fail >= 0 ? 0.92 : 1);
    });
    return;
  }

  const tl = gsap.timeline();
  tl.to('[data-role="axis"]', { opacity: 1, duration: 0.6 }, 0);
  tl.to('[data-role="gate"]', { opacity: 1, duration: 0.7, stagger: 0.18, ease: 'power2.out' }, 0.1);
  tl.to('[data-role="entry-label"]', { opacity: 1, duration: 0.6 }, 0.4);

  // entrada de los puntos (nube)
  dots.forEach((d, i) => {
    const at = 1.0 + i * 0.05;
    tl.to(d, { opacity: 1, duration: 0.4 }, at);
    tl.to(d._halo, { opacity: 0.45, duration: 0.4 }, at);
  });

  // viaje por la criba — pausado
  dots.forEach((d, i) => {
    const fail = parseInt(d.dataset.fail);
    const ex = parseFloat(d.dataset.ex);
    const ey = parseFloat(d.dataset.ey);
    const at = 2.0 + i * 0.085;

    if (fail >= 0) {
      // avanza hasta su filtro y se queda atascado, en rojo (congelado)
      const dur = 1.0 + fail * 0.55;
      tl.to(d,       { attr: { cx: ex, cy: ey }, duration: dur, ease: 'power2.out' }, at);
      tl.to(d._halo, { attr: { cx: ex, cy: ey }, duration: dur, ease: 'power2.out' }, at);
      tl.to(d,       { fill: C_FAIL, duration: 0.4 }, at + dur - 0.2);
      tl.to(d._halo, { opacity: 0, duration: 0.5 }, at + dur);
    } else {
      // atraviesa y se ordena a la salida; al llegar, se valida (verde)
      tl.to(d,       { attr: { cx: ex, cy: ey }, duration: 2.8, ease: 'power1.inOut' }, at);
      tl.to(d._halo, { attr: { cx: ex, cy: ey }, duration: 2.8, ease: 'power1.inOut' }, at);
      tl.to(d,       { fill: C_PASS, duration: 0.5 }, at + 2.4);
    }
  });

  tl.to('[data-role="exit-label"]', { opacity: 1, duration: 0.7, ease: 'power2.out' }, '>-0.6');
}
