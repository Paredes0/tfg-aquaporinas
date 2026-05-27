// scenes/expression-tissues.js — el acuaporinoma en marcha (Bloque IV)
// 22 librerías RNA-seq, 6 tejidos, 32 grupos homeólogos con dominancia desigual
//
// La planta se importa del SVG vectorial fxa_vectorizado.svg (mismo que usa el
// visor eFP del Bloque III: capa1 con g4=fruto rojo, g48=hoja, g133=raíz,
// g144=corona, g230=fruto verde). La yema axilar no está en ese SVG y se
// dibuja como pequeña elipse anclada cerca de la corona.
export function init(slide) {
  if (!window.gsap) {
    const script = document.createElement('script');
    script.src = 'vendor/gsap.min.js';
    document.head.appendChild(script);
  }

  const container = document.createElement('div');
  container.style.cssText = 'position:absolute;inset:0;width:100%;height:100%;z-index:0;display:flex;align-items:center;justify-content:center;';
  slide.insertBefore(container, slide.firstChild);

  const svg = createExpressionSVG();
  container.appendChild(svg);

  // Importar la planta vectorial de forma asíncrona y, cuando esté lista,
  // lanzar la animación. Si fallase la carga del SVG, animamos igualmente
  // sin la planta (los heat-spots, panel y barras quedan reconocibles).
  let started = false;
  const startOnce = () => {
    if (started) return;
    started = true;
    setTimeout(() => animateExpression(svg), 100);
  };

  loadPlantSvg(svg)
    .then(() => startOnce())
    .catch(err => {
      console.warn('[expression-tissues] No se pudo cargar fxa_vectorizado.svg', err);
      startOnce();
    });

  return function dispose() {
    container.remove();
  };
}

// ---- Carga e inserción del SVG vectorial -----------------------------------
// El SVG original tiene viewBox 0 0 210 297 (mm). El bloque útil de la planta
// ocupa aproximadamente x:59-145 (ancho ~86) y:58-180 (alto ~122).
// Lo escalamos x5 y trasladamos (40, -145) para centrarlo en (550, 450) del
// viewBox 0 0 1600 900 de la escena.
const PLANT_TRANSFORM = 'translate(40, -145) scale(5)';

async function loadPlantSvg(sceneSvg) {
  const res = await fetch('assets/fxa_vectorizado.svg', { cache: 'force-cache' });
  if (!res.ok) throw new Error('HTTP ' + res.status);
  const text = await res.text();
  const parser = new DOMParser();
  const doc = parser.parseFromString(text, 'image/svg+xml');
  const errNode = doc.querySelector('parsererror');
  if (errNode) throw new Error('SVG parse error');

  // Tomamos solamente el contenido interno (g id="layer1") para evitar anidar
  // un <svg> dentro de otro con su propio viewBox.
  const layer = doc.getElementById('layer1') || doc.documentElement;
  const plantGroup = document.createElementNS('http://www.w3.org/2000/svg', 'g');
  plantGroup.dataset.role = 'plant-fragaria';
  plantGroup.setAttribute('transform', PLANT_TRANSFORM);
  // Sutil filtro para realzar sobre fondo oscuro (#14182a). Los rellenos
  // originales ya son verdes/marrones oscuros y rojo intenso, así que un
  // brightness moderado + saturación basta.
  plantGroup.setAttribute('style', 'filter: brightness(1.18) saturate(1.05); opacity: 0;');

  // Importamos los hijos directos de layer1 dentro del group de la escena.
  Array.from(layer.children).forEach(child => {
    plantGroup.appendChild(document.importNode(child, true));
  });

  // Insertar la planta DETRÁS de los heat-spots pero por encima del fondo.
  // Como heat-spots se añaden después del placeholder dejado en createExpressionSVG,
  // colocamos la planta antes del primer heat-spot.
  const firstHeat = sceneSvg.querySelector('[data-role="heat-spot"]');
  if (firstHeat) {
    sceneSvg.insertBefore(plantGroup, firstHeat);
  } else {
    sceneSvg.appendChild(plantGroup);
  }

  // Pequeña yema axilar: el SVG eFP no la tiene en la planta principal (se
  // muestra como inset aparte). La dibujamos a mano cerca de la base de la
  // corona como elipse pálida con destello para que el heat-spot tenga algo
  // a lo que anclar visualmente.
  const svgNS = 'http://www.w3.org/2000/svg';
  const budGroup = document.createElementNS(svgNS, 'g');
  budGroup.dataset.role = 'plant-fragaria';
  budGroup.setAttribute('style', 'opacity: 0;');

  const budBody = document.createElementNS(svgNS, 'ellipse');
  budBody.setAttribute('cx', 610);
  budBody.setAttribute('cy', 470);
  budBody.setAttribute('rx', 18);
  budBody.setAttribute('ry', 12);
  budBody.setAttribute('fill', '#c9c3b2');
  budBody.setAttribute('opacity', '0.85');
  budGroup.appendChild(budBody);

  const budStem = document.createElementNS(svgNS, 'path');
  budStem.setAttribute('d', 'M 560 500 Q 585 488 605 472');
  budStem.setAttribute('fill', 'none');
  budStem.setAttribute('stroke', '#7cb342');
  budStem.setAttribute('stroke-width', '3');
  budStem.setAttribute('stroke-linecap', 'round');
  budStem.setAttribute('opacity', '0.75');
  budGroup.appendChild(budStem);

  const budLeaf1 = document.createElementNS(svgNS, 'path');
  budLeaf1.setAttribute('d', 'M 598 462 Q 610 446 622 462');
  budLeaf1.setAttribute('fill', '#6b9e3a');
  budLeaf1.setAttribute('opacity', '0.75');
  budGroup.appendChild(budLeaf1);

  const budLeaf2 = document.createElementNS(svgNS, 'ellipse');
  budLeaf2.setAttribute('cx', 628);
  budLeaf2.setAttribute('cy', 464);
  budLeaf2.setAttribute('rx', 8);
  budLeaf2.setAttribute('ry', 4);
  budLeaf2.setAttribute('fill', '#8bc34a');
  budLeaf2.setAttribute('opacity', '0.65');
  budLeaf2.setAttribute('transform', 'rotate(-18 628 464)');
  budGroup.appendChild(budLeaf2);

  if (firstHeat) {
    sceneSvg.insertBefore(budGroup, firstHeat);
  } else {
    sceneSvg.appendChild(budGroup);
  }
}

function createExpressionSVG() {
  const svgNS = 'http://www.w3.org/2000/svg';
  const svg = document.createElementNS(svgNS, 'svg');
  svg.setAttribute('viewBox', '0 0 1600 900');
  svg.setAttribute('width', '90%');
  svg.setAttribute('height', '90%');

  // ---------- Radial gradients for heat spots (defs) ----------
  const defs = document.createElementNS(svgNS, 'defs');
  const heatColors = {
    ambar:   '#f0b65a',
    crimson: '#ff5577',
    gold:    '#e8a548',
    lightAmbar: '#f5cf94'
  };
  Object.entries(heatColors).forEach(([key, color]) => {
    const grad = document.createElementNS(svgNS, 'radialGradient');
    grad.setAttribute('id', `heat-${key}`);
    grad.setAttribute('cx', '50%');
    grad.setAttribute('cy', '50%');
    grad.setAttribute('r', '50%');
    const s1 = document.createElementNS(svgNS, 'stop');
    s1.setAttribute('offset', '0%');
    s1.setAttribute('stop-color', color);
    s1.setAttribute('stop-opacity', '0.9');
    const s2 = document.createElementNS(svgNS, 'stop');
    s2.setAttribute('offset', '60%');
    s2.setAttribute('stop-color', color);
    s2.setAttribute('stop-opacity', '0.35');
    const s3 = document.createElementNS(svgNS, 'stop');
    s3.setAttribute('offset', '100%');
    s3.setAttribute('stop-color', color);
    s3.setAttribute('stop-opacity', '0');
    grad.appendChild(s1); grad.appendChild(s2); grad.appendChild(s3);
    defs.appendChild(grad);
  });
  svg.appendChild(defs);

  // ---------- 6 tissue heat-spots overlaid on plant ----------
  // Coordenadas re-ancladas a la anatomía del SVG vectorial (transform
  // translate(40,-145) scale(5)). Mapeo mm → pixel: sx = mm_x*5 + 40,
  // sy = mm_y*5 - 145.
  //   leaf       (102.1, 94.8 mm) → (550, 330)
  //   green_fruit (92.7, 81.9 mm) → (504, 265)
  //   red_fruit  (120.2, 110.3 mm) → (640, 405)
  //   crown      (102.2, 130.5 mm) → (550, 510)
  //   roots       (97.7, 159.2 mm) → (530, 650)
  //   axillary_bud (no en SVG, dibujado a mano en (610, 470))
  //
  // Orden de aparición narrativo: HOJA → CORONA → FRUTO ROJO → FRUTO VERDE → YEMA → RAÍZ
  const tissues = [
    { id: 'HOJA',        cx: 550, cy: 330, r: 90, gradient: 'heat-ambar',      heatPct: 47,  tpm: 1240, color: '#f0b65a' },
    { id: 'CORONA',      cx: 550, cy: 510, r: 50, gradient: 'heat-ambar',      heatPct: 35,  tpm: 920,  color: '#f0b65a' },
    { id: 'FRUTO ROJO',  cx: 640, cy: 405, r: 70, gradient: 'heat-crimson',    heatPct: 100, tpm: 3748, color: '#ff5577' },
    { id: 'FRUTO VERDE', cx: 504, cy: 265, r: 38, gradient: 'heat-lightAmbar', heatPct: 30,  tpm: 740,  color: '#f5cf94' },
    { id: 'YEMA',        cx: 610, cy: 470, r: 32, gradient: 'heat-gold',       heatPct: 65,  tpm: 1750, color: '#e8a548' },
    { id: 'RAÍZ',        cx: 530, cy: 650, r: 75, gradient: 'heat-crimson',    heatPct: 89,  tpm: 2380, color: '#ff5577' }
  ];

  tissues.forEach((t, i) => {
    const spot = document.createElementNS(svgNS, 'circle');
    spot.setAttribute('cx', t.cx);
    spot.setAttribute('cy', t.cy);
    spot.setAttribute('r', t.r);
    spot.setAttribute('fill', `url(#${t.gradient})`);
    spot.setAttribute('opacity', 0);
    spot.dataset.role = 'heat-spot';
    spot.dataset.tissueIndex = i;
    svg.appendChild(spot);
  });

  // ---------- Right side info panel ----------
  // Title at x=900 y=200
  const panelTitle = document.createElementNS(svgNS, 'text');
  panelTitle.setAttribute('x', 900);
  panelTitle.setAttribute('y', 200);
  panelTitle.setAttribute('font-family', 'IBM Plex Mono');
  panelTitle.setAttribute('font-size', 14);
  panelTitle.setAttribute('fill', '#8b96b3');
  panelTitle.setAttribute('letter-spacing', '3');
  panelTitle.textContent = 'EXPRESIÓN POR TEJIDO';
  panelTitle.setAttribute('opacity', 0);
  panelTitle.dataset.role = 'panel-title';
  svg.appendChild(panelTitle);

  // Max TPM for normalising the bar width
  const maxTpm = Math.max(...tissues.map(t => t.tpm));
  const barXStart = 1100;
  const barXEnd = 1300;
  const barMaxW = barXEnd - barXStart;

  tissues.forEach((t, i) => {
    const rowY = 240 + i * 65;
    const rowGroup = document.createElementNS(svgNS, 'g');
    rowGroup.dataset.role = 'info-row';
    rowGroup.dataset.tissueIndex = i;
    rowGroup.setAttribute('opacity', 0);

    // Tissue name (left)
    const name = document.createElementNS(svgNS, 'text');
    name.setAttribute('x', 900);
    name.setAttribute('y', rowY);
    name.setAttribute('font-family', 'IBM Plex Mono');
    name.setAttribute('font-size', 16);
    name.setAttribute('fill', t.color);
    name.setAttribute('letter-spacing', '1.5');
    name.textContent = t.id;
    rowGroup.appendChild(name);

    // Bar background (rail)
    const railH = 8;
    const rail = document.createElementNS(svgNS, 'rect');
    rail.setAttribute('x', barXStart);
    rail.setAttribute('y', rowY - railH);
    rail.setAttribute('width', barMaxW);
    rail.setAttribute('height', railH);
    rail.setAttribute('fill', 'rgba(232,240,255,0.06)');
    rail.setAttribute('stroke', 'rgba(232,240,255,0.15)');
    rail.setAttribute('stroke-width', '1');
    rail.setAttribute('rx', '2');
    rowGroup.appendChild(rail);

    // Bar fill (proportional)
    const w = (t.tpm / maxTpm) * barMaxW;
    const bar = document.createElementNS(svgNS, 'rect');
    bar.setAttribute('x', barXStart);
    bar.setAttribute('y', rowY - railH);
    bar.setAttribute('width', w);
    bar.setAttribute('height', railH);
    bar.setAttribute('fill', t.color);
    bar.setAttribute('opacity', '0.85');
    bar.setAttribute('rx', '2');
    rowGroup.appendChild(bar);

    // TPM value (right)
    const tpmText = document.createElementNS(svgNS, 'text');
    tpmText.setAttribute('x', 1320);
    tpmText.setAttribute('y', rowY + 4);
    tpmText.setAttribute('font-family', 'Cormorant Garamond, Georgia, serif');
    tpmText.setAttribute('font-size', 24);
    tpmText.setAttribute('fill', '#e8f0ff');
    tpmText.textContent = t.tpm.toLocaleString('es-ES');
    rowGroup.appendChild(tpmText);

    svg.appendChild(rowGroup);
  });

  // ---------- Bottom: subgenome dominance bars ----------
  const subgenomes = [
    { name: 'A', count: 9,  color: '#a78b6f' },
    { name: 'B', count: 5,  color: '#d8567e' },
    { name: 'C', count: 5,  color: '#4ec5c0' },
    { name: 'D', count: 12, color: '#ff5577' }
  ];
  const sgMaxCount = 12;
  const sgMaxHeight = 80;
  const sgBarW = 60;
  const sgBaseline = 830;

  // Distribute 4 columns evenly from x=300 to x=1300
  const sgX0 = 300, sgX1 = 1300;
  subgenomes.forEach((sg, i) => {
    const cx = sgX0 + (sgX1 - sgX0) * (i / (subgenomes.length - 1));
    const h = (sg.count / sgMaxCount) * sgMaxHeight;

    const bar = document.createElementNS(svgNS, 'rect');
    bar.setAttribute('x', cx - sgBarW / 2);
    bar.setAttribute('y', sgBaseline);
    bar.setAttribute('width', sgBarW);
    bar.setAttribute('height', 0);
    bar.setAttribute('fill', sg.color);
    bar.setAttribute('opacity', '0.85');
    bar.setAttribute('rx', '2');
    bar.dataset.role = 'sg-bar';
    bar.dataset.sgIndex = i;
    bar.dataset.finalY = sgBaseline - h;
    bar.dataset.finalH = h;
    svg.appendChild(bar);

    const lbl = document.createElementNS(svgNS, 'text');
    lbl.setAttribute('x', cx);
    lbl.setAttribute('y', sgBaseline + 30);
    lbl.setAttribute('text-anchor', 'middle');
    lbl.setAttribute('font-family', 'IBM Plex Mono');
    lbl.setAttribute('font-size', 16);
    lbl.setAttribute('fill', sg.color);
    lbl.setAttribute('letter-spacing', '2');
    lbl.textContent = `${sg.name} · ${sg.count}`;
    lbl.setAttribute('opacity', 0);
    lbl.dataset.role = 'sg-label';
    lbl.dataset.sgIndex = i;
    svg.appendChild(lbl);
  });

  return svg;
}

function animateExpression(svg) {
  const gsap = window.gsap;
  if (!gsap) {
    console.warn('[expression-tissues] GSAP no cargado, mostrando estado final sin animación');
    svg.querySelectorAll('[data-role="plant-fragaria"]').forEach(el => el.style.opacity = 1);
    svg.querySelectorAll(
      '[data-role="heat-spot"], [data-role="info-row"], ' +
      '[data-role="panel-title"], [data-role="sg-label"]'
    ).forEach(el => el.setAttribute('opacity', 1));
    // Snap subgenome bars to final height
    svg.querySelectorAll('[data-role="sg-bar"]').forEach(el => {
      el.setAttribute('y', el.dataset.finalY);
      el.setAttribute('height', el.dataset.finalH);
    });
    return;
  }

  const tl = gsap.timeline();

  // 0–0.7s: planta vectorial fade-in suave
  const plantGroups = svg.querySelectorAll('[data-role="plant-fragaria"]');
  tl.to(plantGroups, { opacity: 1, duration: 0.9, ease: 'power2.out' }, 0);

  // Panel title appears early
  tl.to(svg.querySelector('[data-role="panel-title"]'), { opacity: 1, duration: 0.4 }, 0.3);

  // 0.5–3s: 6 heat-spots animate in cascade in narrative order
  // (HOJA → CORONA → FRUTO ROJO → FRUTO VERDE → YEMA → RAÍZ — that's the order they were appended)
  const heatSpots = svg.querySelectorAll('[data-role="heat-spot"]');
  heatSpots.forEach(s => {
    s.style.transformOrigin = `${s.getAttribute('cx')}px ${s.getAttribute('cy')}px`;
    s.style.transformBox = 'fill-box';
  });
  heatSpots.forEach((spot, i) => {
    const at = 0.7 + i * 0.4;
    tl.fromTo(spot,
      { scale: 0.3, opacity: 0 },
      { scale: 1, opacity: 1, duration: 0.7, ease: 'power2.out' },
      at
    );
  });

  // 1–3.5s: 6 info-panel rows fade in, lagging 0.4s behind each spot's start
  const infoRows = svg.querySelectorAll('[data-role="info-row"]');
  infoRows.forEach((row, i) => {
    const at = 0.7 + i * 0.4 + 0.4;
    tl.to(row, { opacity: 1, duration: 0.45, ease: 'power2.out' }, at);
  });

  // 3.5–5s: subgenome bars grow (cascade A → B → C → D)
  const sgBars = svg.querySelectorAll('[data-role="sg-bar"]');
  const sgLabels = svg.querySelectorAll('[data-role="sg-label"]');
  sgBars.forEach((bar, i) => {
    const finalY = parseFloat(bar.dataset.finalY);
    const finalH = parseFloat(bar.dataset.finalH);
    tl.to(bar, {
      attr: { y: finalY, height: finalH },
      duration: 0.6,
      ease: 'power2.out'
    }, 3.7 + i * 0.2);
  });
  tl.to(sgLabels, { opacity: 1, duration: 0.4, stagger: 0.2, ease: 'power2.out' }, 3.9);
}
