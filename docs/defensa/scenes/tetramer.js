// scenes/tetramer.js — Mol* viewer con PDB 1Z98 (SoPIP2;1 espinaca)
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

  const container = document.createElement('div');
  container.style.cssText = 'position:absolute;inset:0;width:100%;height:100%;z-index:0;background:#14182a;';
  slide.insertBefore(container, slide.firstChild);

  let molstarPlugin;

  import('./shared/molstar-wrapper.js?v=20260527-02').then(async ({ createViewer }) => {
    molstarPlugin = await createViewer(container, {
      pdbUrl: 'assets/pdb/1Z98.pdb',
      backgroundColor: 0x14182a,
      autoRotate: true
    });
  }).catch(err => {
    console.error('[tetramer] init fallo', err);
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

function showFallback(slide) {
  const img = document.createElement('img');
  img.src = 'assets/biorender/fallback-tetramer.png';
  img.style.cssText = 'position:absolute;inset:0;width:100%;height:100%;object-fit:contain;z-index:0;background:#14182a;';
  img.onerror = () => { img.style.display = 'none'; };
  slide.insertBefore(img, slide.firstChild);
}
