// scenes/shared/molstar-wrapper.js — envuelve la API de Mol* viewer
// Soporta múltiples viewers en la misma página con backgrounds distintos
// usando CSS scoped al container vía clase única (molstar-bg-XXXXXX).

function injectMolstarBgCss(hex) {
  const cleanHex = hex.replace('#', '').toLowerCase();
  const cssId = `molstar-bg-${cleanHex}`;
  if (document.getElementById(cssId)) return `molstar-bg-${cleanHex}`;
  const style = document.createElement('style');
  style.id = cssId;
  style.textContent = `
    .${cssId} .msp-plugin,
    .${cssId} .msp-plugin-content,
    .${cssId} .msp-viewport,
    .${cssId} .msp-viewport-controls,
    .${cssId} .msp-canvas3d,
    .${cssId} canvas.msp-viewport-controls-info,
    .${cssId} .msp-layout-expanded,
    .${cssId} .msp-layout-standard { background: ${hex} !important; }
  `;
  document.head.appendChild(style);
  return cssId;
}

// Mol* viewer.js intenta consultar un endpoint PDBe local hardcoded
// (http://localhost:9000/v2/list_entries/...) que no servimos. Suprime esos
// 6 errores rojos en consola sin afectar el rendering del PDB.
function installPdbeFetchSilencer() {
  if (window.__molstarFetchSilenced) return;
  const origFetch = window.fetch.bind(window);
  const matches = (s) => s.includes('localhost:9000') || s.includes('/v2/list_entries') || s.includes('pdbe-kb');
  window.fetch = function patched(input, init) {
    let href = '';
    try {
      if (typeof input === 'string') href = input;
      else if (input instanceof URL) href = input.href;
      else if (input && typeof input.url === 'string') href = input.url;
      else href = String(input);
    } catch (e) {}
    if (matches(href)) {
      return Promise.resolve(new Response('[]', { status: 200, headers: { 'Content-Type': 'application/json' } }));
    }
    return origFetch(input, init);
  };
  // Mol* may also use XMLHttpRequest for some endpoints — patch open() too.
  const OrigXHROpen = window.XMLHttpRequest.prototype.open;
  window.XMLHttpRequest.prototype.open = function (method, url, ...rest) {
    if (typeof url === 'string' && matches(url)) {
      this.__silencedUrl = url;
      return OrigXHROpen.call(this, method, 'data:application/json,[]', ...rest);
    }
    return OrigXHROpen.call(this, method, url, ...rest);
  };
  window.__molstarFetchSilenced = true;
}

export async function createViewer(container, { pdbUrl, backgroundColor = 0x14182a, autoRotate = true } = {}) {
  // 1. CSS scoped al container con clase única por color
  const hex = '#' + backgroundColor.toString(16).padStart(6, '0');
  const bgClass = injectMolstarBgCss(hex);
  container.classList.add(bgClass);
  container.style.background = hex;
  installPdbeFetchSilencer();

  // 2. Cargar Mol* viewer.js si no está ya
  if (!window.molstar) {
    await new Promise((resolve, reject) => {
      const s = document.createElement('script');
      s.src = 'vendor/molstar/molstar.js';
      s.onload = resolve;
      s.onerror = reject;
      document.head.appendChild(s);
    });
  }

  // 3. Constructor de Viewer
  const viewer = await window.molstar.Viewer.create(container, {
    layoutIsExpanded: false,
    layoutShowControls: false,
    layoutShowRemoteState: false,
    layoutShowSequence: false,
    layoutShowLog: false,
    layoutShowLeftPanel: false,
    viewportShowExpand: false,
    viewportShowSelectionMode: false,
    viewportShowAnimation: false,
    backgroundColor
  });

  // 4. Forzar fondo transparente del canvas Mol* → deja ver el bg del container/slide
  if (viewer.plugin?.canvas3d) {
    try {
      viewer.plugin.canvas3d.setProps({
        renderer: {
          backgroundColor: backgroundColor,
          transparentBackground: true
        }
      });
    } catch (e) {
      console.warn('[molstar-wrapper] setProps renderer fallo:', e);
    }
  }

  // pdbUrl puede ser null/undefined → no cargar nada (MVS u otro proceso lo hará)
  if (pdbUrl) {
    await viewer.loadStructureFromUrl(pdbUrl, 'pdb', false);
  }

  // 5. Tras cargar la estructura, re-aplicar para asegurar (algunos builds de Mol* resetean al load)
  if (viewer.plugin?.canvas3d) {
    try {
      viewer.plugin.canvas3d.setProps({
        renderer: {
          backgroundColor: backgroundColor,
          transparentBackground: true
        }
      });
    } catch (e) {}
  }

  if (autoRotate) {
    viewer.plugin.canvas3d?.setProps({
      trackball: { animate: { name: 'spin', params: { speed: 0.5 } } }
    });
  }

  return viewer;
}
