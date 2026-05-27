// scenes/aquaporin-monomer.js — Mol* viewer del monómero del PDB 1Z98 (chain A)
// Tema LIGHT (fondo blanco), vista lateral con auto-rotación lenta para apreciar las 6 TMH.
// Anotaciones HTML superpuestas: 6 TMH · 2 NPA · ar/R.
import { isWebGLAvailable } from './shared/three-helpers.js';

export function init(slide) {
  if (!isWebGLAvailable()) {
    showFallback(slide);
    return () => {};
  }

  const cssLink = document.createElement('link');
  cssLink.rel = 'stylesheet';
  cssLink.href = 'vendor/molstar/molstar.css';
  document.head.appendChild(cssLink);

  // Container target dentro del slide (donde aparecerá el viewer)
  const target = slide.querySelector('[data-aquaporin-target]');
  if (!target) {
    console.warn('[aquaporin-monomer] no [data-aquaporin-target] in slide');
    return () => {};
  }

  // No borrar overlays HTML (leyenda, pie PDB) — solo añadir el viewer detrás
  const container = document.createElement('div');
  container.style.cssText = 'position:absolute;inset:0;width:100%;height:100%;background:#ffffff;border-radius:4px;overflow:hidden;z-index:1;';
  target.insertBefore(container, target.firstChild);

  let molstarPlugin;

  import('./shared/molstar-wrapper.js?v=20260527-17').then(async ({ createViewer }) => {
    // Crear viewer vacío (sin cargar PDB). MVS se encarga de cargar + colorear.
    molstarPlugin = await createViewer(container, {
      pdbUrl: null,
      backgroundColor: 0xffffff,
      autoRotate: true
    });
    window.__aqpViewer = molstarPlugin;
    try {
      await applyResidueColoring(molstarPlugin);
    } catch (e) {
      console.warn('[aquaporin-monomer] MVS coloreado fallo, recargando PDB sin colorear:', e);
      try {
        await molstarPlugin.loadStructureFromUrl('assets/pdb/1Z98_chainA.pdb', 'pdb', false);
      } catch (e2) {}
    }
    // Rotación lenta sobre el eje "up" de la cámara (eje del poro Z global → vertical en pantalla)
    try {
      molstarPlugin.plugin.canvas3d?.setProps({
        trackball: { animate: { name: 'spin', params: { speed: 0.25 } } }
      });
    } catch (e) {}
  }).catch(err => {
    console.error('[aquaporin-monomer] init fallo', err);
    showFallback(slide);
  });

  return function dispose() {
    if (molstarPlugin) {
      try { molstarPlugin.dispose(); } catch (e) {}
    }
    container.remove();
    cssLink.remove();
  };
}

async function applyResidueColoring(viewer) {
  const plugin = viewer.plugin;
  const mvsExt = window.molstar?.PluginExtensions?.mvs;
  if (!mvsExt) throw new Error('MolViewSpec extension no disponible');

  // MVS carga el PDB y crea todas las representaciones — no usamos plugin.clear() porque
  // el viewer está vacío (createViewer con pdbUrl:null).

  // Construir escena MVS con coloreado por elementos funcionales del PDB 1Z98 (SoPIP2;1):
  //   · TMH (6 hélices transmembrana): azul · cartoon
  //   · Hélices reentrantes HB (101-111) y HE (222-234): verde · cartoon
  //   · 2 NPA (101-103 y 222-224): rojo · ball-and-stick
  //   · Filtro ar/R (Arg 225): ámbar · ball-and-stick
  const b = mvsExt.MVSData.createBuilder();
  const s = b.download({ url: 'assets/pdb/1Z98_chainA.pdb' })
    .parse({ format: 'pdb' })
    .modelStructure({});

  // Función helper: crear array de residuos a partir de rangos [start,end]
  const range = (ranges) => {
    const out = [];
    for (const [start, end] of ranges) {
      for (let i = start; i <= end; i++) out.push({ label_seq_id: i });
    }
    return out;
  };

  // Capa 1 · LOOPS y extremos N/C en cartoon gris claro
  // (residuos que NO son TMH ni reentrantes; SoPIP2;1 chain A va de 24 a 274)
  const loopRanges = [
    [24, 32],   // N-terminal antes de TMH1
    [65, 71],   // loop A (entre TMH1 y TMH2)
    [94, 100],  // loop B inicial (antes de HB)
    [112, 114], // loop B final (después de HB, antes de TMH3)
    [149, 159], // loop C (entre TMH3 y TMH4)
    [183, 197], // loop D (entre TMH4 y TMH5)
    [219, 221], // loop E inicial (antes de HE)
    [235, 242], // loop E final (después de HE, antes de TMH6)
    [264, 274]  // C-terminal después de TMH6
  ];
  s.component({ selector: range(loopRanges) })
    .representation({ type: 'cartoon' })
    .color({ color: '#cdd2d8' });

  // Capa 2 · 6 TMH cartoon azul
  // TMH1: 33-64 · TMH2: 72-93 · TMH3: 115-148 · TMH4: 160-182 · TMH5: 198-218 · TMH6: 243-263
  const tmhRanges = [
    [33, 64], [72, 93], [115, 148], [160, 182], [198, 218], [243, 263]
  ];
  s.component({ selector: range(tmhRanges) })
    .representation({ type: 'cartoon' })
    .color({ color: '#3498DB' });

  // Capa 3 · 2 hélices reentrantes cartoon verde (HB 101-111, HE 222-234)
  s.component({ selector: range([[101, 111], [222, 234]]) })
    .representation({ type: 'cartoon' })
    .color({ color: '#2ECC71' });

  // Capa 4 · 2 NPA · ball-and-stick rojo (residuos 101-103 + 222-224)
  const npa = [101, 102, 103, 222, 223, 224].map(i => ({ label_seq_id: i }));
  s.component({ selector: npa })
    .representation({ type: 'ball_and_stick' })
    .color({ color: '#E74C3C' });

  // Capa 5 · filtro ar/R · ball-and-stick ámbar
  // Phe81 (TMH2), His210 (loop E), Thr213, Arg225 (Tornroth-Horsefield et al. 2005)
  const arR = [81, 210, 213, 225].map(i => ({ label_seq_id: i }));
  s.component({ selector: arR })
    .representation({ type: 'ball_and_stick' })
    .color({ color: '#F39C12' });

  // Capa 6 · residuos del gating · 5 residuos mecanísticamente críticos
  // Asp28 (capping Ca²⁺), Ser115 + Ser188 (fosforilación → abre),
  // His193 (sensor pH del loop D), Leu197 (tapón hidrofóbico que cierra el poro)
  // (Tornroth-Horsefield et al. 2005, Nature · refinado por Nyblom et al. 2009, JMB)
  const gating = [28, 115, 188, 193, 197].map(i => ({ label_seq_id: i }));
  s.component({ selector: gating })
    .representation({ type: 'ball_and_stick' })
    .color({ color: '#9B59B6' });

  // Cámara · vista lateral con el eje del poro (Z global) vertical en pantalla.
  // Convención biológica: extracelular ARRIBA (ar/R), citoplasma ABAJO (compuerta).
  // En el PDB 1Z98 el lado extracelular está en Z negativo, así que up=(0,0,-1)
  // pone Z negativo arriba → extracelular arriba.
  b.camera({
    target: [32, 30, 2],
    position: [115, 30, 2],
    up: [0, 0, -1]
  });

  await mvsExt.loadMVS(plugin, b.getState());
}

function showFallback(slide) {
  const target = slide.querySelector('[data-aquaporin-target]');
  if (!target) return;
  const note = document.createElement('div');
  note.style.cssText = 'display:flex;align-items:center;justify-content:center;width:100%;height:100%;color:#8b96b3;font-family:IBM Plex Mono,monospace;font-size:14px;';
  note.textContent = '[Monómero del PDB 1Z98 no disponible — WebGL requerido]';
  target.appendChild(note);
}
