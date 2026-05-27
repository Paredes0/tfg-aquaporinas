// scenes/forest.js — Bloque III separator (slide 13).
//
// V5 (2026-05-26): muestra la imagen REAL del árbol filogenético circular
// generado externamente (arbol_con_heatmap_2.png, cropeado a
// arbol_circular_sin_leyendas.png sin la leyenda del autor original).
// La animación es un fade-in suave; la información detallada está en slide 14
// (tree-static.js renderiza el .treefile en rectangular). Aquí solo sirve
// como "billboard" cinematográfico al inicio del bloque.

const SUBFAMILY_LEGEND = [
  { name: 'PIP', color: '#E74C3C', count: 32 },
  { name: 'TIP', color: '#3498DB', count: 34 },
  { name: 'NIP', color: '#2ECC71', count: 37 },
  { name: 'SIP', color: '#F39C12', count: 12 },
  { name: 'XIP', color: '#9B59B6', count: 6 }
];

export function init(slide) {
  const container = document.createElement('div');
  // Padding: top deja sitio al scene-overlay-bar (~250px), right a la leyenda
  // de subfamilias (~280px), bottom al corner-note (~120px), left mínimo.
  container.style.cssText = `
    position: absolute;
    inset: 0;
    width: 100%;
    height: 100%;
    z-index: 0;
    display: flex;
    align-items: center;
    justify-content: center;
    padding: 180px 240px 80px 60px;
    box-sizing: border-box;
    background: #14182a;
  `;
  slide.insertBefore(container, slide.firstChild);

  // Imagen del árbol — fit al espacio disponible (ya descontado padding)
  const img = document.createElement('img');
  img.src = 'assets/arbol_circular_sin_leyendas.png';
  img.alt = 'Árbol filogenético circular · 281 secuencias · Q.PLANT+R6';
  img.style.cssText = `
    max-width: 100%;
    max-height: 100%;
    object-fit: contain;
    opacity: 0;
    transform: scale(0.94);
    transition: opacity 1.2s ease, transform 1.4s cubic-bezier(.2,.7,.2,1);
  `;
  container.appendChild(img);

  // Leyenda de subfamilias (columna derecha) — la crop quitó la leyenda original,
  // así que la añadimos como overlay limpio del deck.
  const legend = document.createElement('div');
  legend.style.cssText = `
    position: absolute;
    right: 56px;
    top: 50%;
    transform: translateY(-50%);
    display: flex;
    flex-direction: column;
    gap: 18px;
    z-index: 2;
    opacity: 0;
    transition: opacity 0.8s ease 0.7s;
  `;
  SUBFAMILY_LEGEND.forEach(sf => {
    const row = document.createElement('div');
    row.style.cssText = 'display:flex; align-items:center; gap:14px;';

    const dot = document.createElement('span');
    dot.style.cssText = `
      width: 18px; height: 18px; border-radius: 3px;
      background: ${sf.color};
      box-shadow: 0 0 14px ${sf.color}88;
    `;

    const label = document.createElement('span');
    label.style.cssText = `
      font-family: 'IBM Plex Mono', monospace;
      font-size: 18px;
      color: #e8f0ff;
      letter-spacing: 0.08em;
    `;
    label.innerHTML =
      '<strong style="color:' + sf.color + '">' + sf.name + '</strong>' +
      ' <span style="color:#8b96b3; font-size:14px;">· ' + sf.count + '</span>';

    row.appendChild(dot);
    row.appendChild(label);
    legend.appendChild(row);
  });
  container.appendChild(legend);

  // Anotación inferior izquierda — explica los anillos para el tribunal.
  // bottom: 100px para no chocar con el .footer-rule del slide (bottom: 28px + altura).
  const corner = document.createElement('div');
  corner.style.cssText = `
    position: absolute;
    left: 56px;
    bottom: 100px;
    z-index: 2;
    opacity: 0;
    transition: opacity 0.8s ease 0.9s;
    font-family: 'IBM Plex Mono', monospace;
    font-size: 13px;
    color: #8b96b3;
    letter-spacing: 0.1em;
    max-width: 360px;
    line-height: 1.5;
  `;
  corner.innerHTML =
    '<div style="margin-bottom: 8px;">281 secuencias · 430 sitios · Q.PLANT+R6</div>' +
    '<div style="color: #5a6485; font-size: 12px;">log L = −45.149,26</div>' +
    '<div style="color: #5a6485; font-size: 11px; margin-top: 10px;">' +
    'Anillo interno: subfamilia. Anillo externo: TPM por tejido (heatmap).' +
    '</div>';
  container.appendChild(corner);

  // Disparar fade-in tras dos rAFs (el navegador necesita registrar el initial state)
  requestAnimationFrame(() => {
    requestAnimationFrame(() => {
      img.style.opacity = '1';
      img.style.transform = 'scale(1)';
      legend.style.opacity = '1';
      corner.style.opacity = '1';
    });
  });

  // Fallback si la imagen no carga
  img.onerror = () => {
    img.style.display = 'none';
    const error = document.createElement('div');
    error.style.cssText = `
      color: #8b96b3;
      font-family: 'IBM Plex Mono', monospace;
      font-size: 15px;
      opacity: 0.7;
      text-align: center;
      padding: 40px;
      max-width: 600px;
    `;
    error.innerHTML =
      '[árbol filogenético no disponible]<br/>' +
      '<span style="font-size:12px; color:#5a6485;">' +
      'Buscado en assets/arbol_circular_sin_leyendas.png' +
      '</span>';
    container.appendChild(error);
  };

  return function dispose() {
    container.remove();
  };
}
