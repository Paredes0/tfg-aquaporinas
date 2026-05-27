// scenes/tree-static.js — arbol filogenetico REAL renderizado en la slide 14
// (Bloque III, "Arbol filogenetico final"). Es la version estatica que
// sustituye al SVG estilizado a mano que habia en el HTML.
//
// Renderiza un cladograma rectangular (izquierda → derecha) con:
//   - Esqueleto en gris oscuro sobre fondo blanco (slide 14 es light)
//   - Ramas coloreadas por subfamilia (PIP rojo, TIP azul, NIP verde,
//     SIP naranja, XIP violeta)
//   - Etiquetas de subfamilia + cuenta a la derecha
//   - Tres controles literarios marcados como puntos con etiqueta:
//     FaPIP1;1, FaPIP2;1, FaNIP1;1
//
// La escena se inserta dentro de un contenedor del slide marcado con
// `data-tree-target`. Si no encuentra ese contenedor, hace append al slide.
//
// Animacion minimal: fade-in suave de 1 s sobre todo el grupo del arbol.
import {
  parseNewick, fetchNewick, assignSubfamilies,
  postOrderLeaves,
  SUBFAMILIES, SF_BY_NAME, MUTED_LIGHT
} from './shared/newick.js';

const VIEW_W = 800;
const VIEW_H = 580;

// Margenes para el dibujo del arbol dentro del viewBox
const X_MIN = 40;
const X_MAX = 640;     // dejamos espacio a la derecha para etiquetas/controles
const Y_MIN = 30;
const Y_MAX = 540;     // dejamos espacio abajo para barra de escala

// Color de las ramas mezcladas/profundas sobre fondo claro
const SKELETON = 'rgba(20,24,42,0.32)';
const SKELETON_LIGHT = 'rgba(20,24,42,0.20)';

// Nombres exactos de los controles tal como aparecen en el .treefile
// (basado en convenciones del proyecto: FaPIP1_1, FaPIP2_1, FaNIP1_1 con
// posibles variantes; usamos detectores tolerantes).
const CONTROL_MATCHERS = [
  { label: 'FaPIP1;1', match: name => /^FaPIP1[_;.]1$/.test(name) || /FaPIP1_1\b/.test(name) },
  { label: 'FaPIP2;1', match: name => /^FaPIP2[_;.]1$/.test(name) || /FaPIP2_1\b/.test(name) },
  { label: 'FaNIP1;1', match: name => /^FaNIP1[_;.]1$/.test(name) || /FaNIP1_1\b/.test(name) }
];

export function init(slide) {
  // Localizar contenedor objetivo dentro del slide
  const target = slide.querySelector('[data-tree-target]') || slide;

  const container = document.createElement('div');
  container.style.cssText = 'width:100%;height:100%;display:flex;align-items:center;justify-content:center;';
  target.appendChild(container);

  const placeholder = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
  placeholder.setAttribute('viewBox', `0 0 ${VIEW_W} ${VIEW_H}`);
  placeholder.setAttribute('width', '100%');
  placeholder.style.maxWidth = '100%';
  placeholder.style.height = 'auto';
  container.appendChild(placeholder);

  let disposed = false;

  loadAndRender(container, placeholder).catch(err => {
    console.warn('[tree-static] No se pudo cargar el arbol real, usando placeholder', err);
    if (disposed) return;
    if (placeholder.parentNode === container) container.removeChild(placeholder);
    const fallback = createPlaceholderSVG();
    container.appendChild(fallback);
  });

  return function dispose() {
    disposed = true;
    container.remove();
  };
}

async function loadAndRender(container, placeholder) {
  const newick = await fetchNewick();
  const root = parseNewick(newick);
  assignSubfamilies(root, { muted: SKELETON });

  const layout = layoutRectangular(root);
  const svg = createTreeSVG(root, layout);

  if (placeholder.parentNode === container) container.removeChild(placeholder);
  container.appendChild(svg);

  // Fade-in muy sutil de todo el grupo
  if (window.gsap) {
    window.gsap.fromTo(svg, { opacity: 0 }, { opacity: 1, duration: 1.0, ease: 'power2.out' });
  } else {
    svg.style.opacity = '1';
  }
}

// =============================================================================
// Layout rectangular (clasico cladograma left → right)
// =============================================================================
function layoutRectangular(root) {
  const leaves = [];
  postOrderLeaves(root, leaves);
  const n = leaves.length;

  // y por orden post-order
  leaves.forEach((leaf, idx) => {
    leaf._idx = idx;
    leaf.y = Y_MIN + (Y_MAX - Y_MIN) * (n === 1 ? 0.5 : idx / (n - 1));
  });

  // profundidad acumulada
  function computeDepth(node, depth) {
    node.depth = depth;
    for (const c of node.children) computeDepth(c, depth + (c.length || 0));
  }
  computeDepth(root, 0);

  let maxDepth = 0;
  (function findMax(node) {
    if (node.depth > maxDepth) maxDepth = node.depth;
    for (const c of node.children) findMax(c);
  })(root);
  if (maxDepth <= 0) maxDepth = 1;

  const xScale = (X_MAX - X_MIN) / maxDepth;
  (function assignX(node) {
    node.x = X_MIN + node.depth * xScale;
    for (const c of node.children) assignX(c);
  })(root);

  (function assignY(node) {
    if (node.children.length === 0) return node.y;
    let sum = 0;
    for (const c of node.children) sum += assignY(c);
    node.y = sum / node.children.length;
    return node.y;
  })(root);

  return { leaves, root, maxDepth, n };
}

// =============================================================================
// SVG render
// =============================================================================
function createTreeSVG(root, layout) {
  const svgNS = 'http://www.w3.org/2000/svg';
  const svg = document.createElementNS(svgNS, 'svg');
  svg.setAttribute('viewBox', `0 0 ${VIEW_W} ${VIEW_H}`);
  svg.setAttribute('width', '100%');
  svg.style.maxWidth = '100%';
  svg.style.height = 'auto';
  svg.style.opacity = '0';

  // Grupos por subfamilia + uno para el esqueleto
  const groups = {};
  SUBFAMILIES.forEach(sf => {
    const g = document.createElementNS(svgNS, 'g');
    g.dataset.role = 'sf-group';
    g.dataset.sfName = sf.name;
    svg.appendChild(g);
    groups[sf.name] = g;
  });
  const muteG = document.createElementNS(svgNS, 'g');
  muteG.dataset.role = 'mute-group';
  svg.appendChild(muteG);

  // Dibujar aristas
  function drawEdges(node) {
    if (node.children.length === 0) return;
    // Linea vertical entre los hijos extremos
    const ys = node.children.map(c => c.y);
    const yMin = Math.min(...ys), yMax = Math.max(...ys);
    if (yMax - yMin > 0.5) {
      const v = document.createElementNS(svgNS, 'line');
      v.setAttribute('x1', node.x); v.setAttribute('x2', node.x);
      v.setAttribute('y1', yMin);   v.setAttribute('y2', yMax);
      v.setAttribute('stroke', node.branchColor);
      v.setAttribute('stroke-width', node.branchColor === SKELETON ? 1.0 : 1.5);
      v.setAttribute('stroke-linecap', 'round');
      const tgt = node.branchSubfamily ? groups[node.branchSubfamily] : muteG;
      tgt.appendChild(v);
    }
    for (const c of node.children) {
      const h = document.createElementNS(svgNS, 'line');
      h.setAttribute('x1', node.x); h.setAttribute('x2', c.x);
      h.setAttribute('y1', c.y);    h.setAttribute('y2', c.y);
      h.setAttribute('stroke', c.branchColor);
      h.setAttribute('stroke-width', c.branchColor === SKELETON ? 1.0 : 1.5);
      h.setAttribute('stroke-linecap', 'round');
      const tgt = c.branchSubfamily ? groups[c.branchSubfamily] : muteG;
      tgt.appendChild(h);
      drawEdges(c);
    }
  }
  drawEdges(root);

  // ---- Etiquetas por subfamilia ------------------------------------------
  // Calcular el rango y por subfamilia + el x maximo (para alinear)
  const rangesBySf = {};
  for (const leaf of layout.leaves) {
    const sf = leaf.branchSubfamily;
    if (!sf) continue;
    if (!rangesBySf[sf]) rangesBySf[sf] = { yMin: leaf.y, yMax: leaf.y, xMax: leaf.x };
    rangesBySf[sf].yMin = Math.min(rangesBySf[sf].yMin, leaf.y);
    rangesBySf[sf].yMax = Math.max(rangesBySf[sf].yMax, leaf.y);
    rangesBySf[sf].xMax = Math.max(rangesBySf[sf].xMax, leaf.x);
  }

  const labelsLayer = document.createElementNS(svgNS, 'g');
  svg.appendChild(labelsLayer);
  SUBFAMILIES.forEach(sf => {
    const r = rangesBySf[sf.name];
    if (!r) return;
    // Linea guia vertical
    const guide = document.createElementNS(svgNS, 'line');
    guide.setAttribute('x1', X_MAX + 18); guide.setAttribute('x2', X_MAX + 18);
    guide.setAttribute('y1', r.yMin);     guide.setAttribute('y2', r.yMax);
    guide.setAttribute('stroke', sf.color);
    guide.setAttribute('stroke-width', '2.5');
    guide.setAttribute('stroke-linecap', 'round');
    guide.setAttribute('opacity', '0.85');
    labelsLayer.appendChild(guide);

    const label = document.createElementNS(svgNS, 'text');
    label.setAttribute('x', X_MAX + 32);
    label.setAttribute('y', (r.yMin + r.yMax) / 2 + 5);
    label.setAttribute('font-family', 'IBM Plex Mono, monospace');
    label.setAttribute('font-size', 16);
    label.setAttribute('font-weight', '600');
    label.setAttribute('letter-spacing', '2');
    label.setAttribute('fill', sf.color);
    label.textContent = `${sf.name} · ${sf.count}`;
    labelsLayer.appendChild(label);
  });

  // ---- Controles literarios -----------------------------------------------
  const controlsLayer = document.createElementNS(svgNS, 'g');
  svg.appendChild(controlsLayer);
  CONTROL_MATCHERS.forEach((ctrl, i) => {
    const leaf = layout.leaves.find(l => ctrl.match(l.name));
    if (!leaf) return;
    // Punto destacado (anillo)
    const ring = document.createElementNS(svgNS, 'circle');
    ring.setAttribute('cx', leaf.x);
    ring.setAttribute('cy', leaf.y);
    ring.setAttribute('r', 4.5);
    ring.setAttribute('fill', 'none');
    ring.setAttribute('stroke', '#14182a');
    ring.setAttribute('stroke-width', '1');
    controlsLayer.appendChild(ring);

    const dot = document.createElementNS(svgNS, 'circle');
    dot.setAttribute('cx', leaf.x);
    dot.setAttribute('cy', leaf.y);
    dot.setAttribute('r', 2.5);
    dot.setAttribute('fill', '#14182a');
    controlsLayer.appendChild(dot);

    // Label (a la izquierda de la hoja para no chocar con la etiqueta SF)
    const lbl = document.createElementNS(svgNS, 'text');
    lbl.setAttribute('x', leaf.x - 8);
    lbl.setAttribute('y', leaf.y - 6);
    lbl.setAttribute('text-anchor', 'end');
    lbl.setAttribute('font-family', 'IBM Plex Mono, monospace');
    lbl.setAttribute('font-size', 10);
    lbl.setAttribute('letter-spacing', '1');
    lbl.setAttribute('fill', '#3a4055');
    lbl.textContent = ctrl.label;
    controlsLayer.appendChild(lbl);
  });

  // ---- Encabezado de cuenta total ----------------------------------------
  const header = document.createElementNS(svgNS, 'text');
  header.setAttribute('x', X_MIN);
  header.setAttribute('y', 18);
  header.setAttribute('font-family', 'IBM Plex Mono, monospace');
  header.setAttribute('font-size', 11);
  header.setAttribute('letter-spacing', '2.5');
  header.setAttribute('fill', '#5a6485');
  header.textContent = `${layout.n} SECUENCIAS · Q.PLANT+R6 · MAFFT · ClipKIT · IQ-TREE 3`;
  svg.appendChild(header);

  // ---- Barra de escala (esquina inferior izquierda) ----------------------
  const scaleG = document.createElementNS(svgNS, 'g');
  scaleG.setAttribute('transform', `translate(${X_MIN}, ${Y_MAX + 20})`);
  const scaleLine = document.createElementNS(svgNS, 'line');
  scaleLine.setAttribute('x1', 0); scaleLine.setAttribute('x2', 60);
  scaleLine.setAttribute('y1', 0); scaleLine.setAttribute('y2', 0);
  scaleLine.setAttribute('stroke', '#5a6485');
  scaleLine.setAttribute('stroke-width', '1');
  scaleG.appendChild(scaleLine);
  const t1 = document.createElementNS(svgNS, 'line');
  t1.setAttribute('x1', 0); t1.setAttribute('x2', 0);
  t1.setAttribute('y1', -4); t1.setAttribute('y2', 4);
  t1.setAttribute('stroke', '#5a6485'); t1.setAttribute('stroke-width', '1');
  scaleG.appendChild(t1);
  const t2 = document.createElementNS(svgNS, 'line');
  t2.setAttribute('x1', 60); t2.setAttribute('x2', 60);
  t2.setAttribute('y1', -4); t2.setAttribute('y2', 4);
  t2.setAttribute('stroke', '#5a6485'); t2.setAttribute('stroke-width', '1');
  scaleG.appendChild(t2);
  const scaleLabel = document.createElementNS(svgNS, 'text');
  scaleLabel.setAttribute('x', 68);
  scaleLabel.setAttribute('y', 4);
  scaleLabel.setAttribute('font-family', 'IBM Plex Mono, monospace');
  scaleLabel.setAttribute('font-size', 10);
  scaleLabel.setAttribute('fill', '#5a6485');
  scaleLabel.textContent = '0,1 sub/sitio';
  scaleG.appendChild(scaleLabel);
  svg.appendChild(scaleG);

  return svg;
}

// =============================================================================
// FALLBACK · placeholder esquematico si la carga del .treefile falla
// =============================================================================
function createPlaceholderSVG() {
  const svgNS = 'http://www.w3.org/2000/svg';
  const svg = document.createElementNS(svgNS, 'svg');
  svg.setAttribute('viewBox', `0 0 ${VIEW_W} ${VIEW_H}`);
  svg.setAttribute('width', '100%');
  svg.style.maxWidth = '100%';
  svg.style.height = 'auto';

  // Un cladograma mini-estilizado (similar al SVG original que sustituimos):
  // 5 ramas horizontales + trunk + labels. Sobre blanco.
  const baseY = 290;
  const trunkX = 100;
  const branchEndX = 600;

  const trunk = document.createElementNS(svgNS, 'line');
  trunk.setAttribute('x1', 50); trunk.setAttribute('y1', baseY);
  trunk.setAttribute('x2', trunkX); trunk.setAttribute('y2', baseY);
  trunk.setAttribute('stroke', SKELETON); trunk.setAttribute('stroke-width', '1.5');
  svg.appendChild(trunk);

  const ys = [80, 180, 290, 400, 490];
  SUBFAMILIES.forEach((sf, i) => {
    const y = ys[i];
    const conn = document.createElementNS(svgNS, 'line');
    conn.setAttribute('x1', trunkX); conn.setAttribute('y1', baseY);
    conn.setAttribute('x2', trunkX); conn.setAttribute('y2', y);
    conn.setAttribute('stroke', SKELETON); conn.setAttribute('stroke-width', '1.5');
    svg.appendChild(conn);

    const horiz = document.createElementNS(svgNS, 'line');
    horiz.setAttribute('x1', trunkX); horiz.setAttribute('y1', y);
    horiz.setAttribute('x2', branchEndX); horiz.setAttribute('y2', y);
    horiz.setAttribute('stroke', sf.color); horiz.setAttribute('stroke-width', '2');
    svg.appendChild(horiz);

    const label = document.createElementNS(svgNS, 'text');
    label.setAttribute('x', branchEndX + 18);
    label.setAttribute('y', y + 6);
    label.setAttribute('font-family', 'IBM Plex Mono, monospace');
    label.setAttribute('font-size', 16);
    label.setAttribute('font-weight', '600');
    label.setAttribute('letter-spacing', '2');
    label.setAttribute('fill', sf.color);
    label.textContent = `${sf.name} · ${sf.count}`;
    svg.appendChild(label);
  });

  return svg;
}
