// scenes/shared/three-helpers.js — utilidades comunes a las escenas Three.js
export function isWebGLAvailable() {
  try {
    const canvas = document.createElement('canvas');
    return !!(window.WebGLRenderingContext && (canvas.getContext('webgl2') || canvas.getContext('webgl')));
  } catch (e) {
    return false;
  }
}

export function createCanvas(slide, { fullscreen = true } = {}) {
  const canvas = document.createElement('canvas');
  canvas.style.cssText = fullscreen
    ? 'position:absolute;inset:0;width:100%;height:100%;z-index:0;display:block;'
    : 'width:100%;height:100%;display:block;';
  slide.insertBefore(canvas, slide.firstChild);
  return canvas;
}

export function resizeRenderer(renderer, camera) {
  const parent = renderer.domElement.parentElement;
  const w = parent.clientWidth;
  const h = parent.clientHeight;
  renderer.setSize(w, h, false);
  renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
  if (camera.isPerspectiveCamera) {
    camera.aspect = w / h;
    camera.updateProjectionMatrix();
  }
}

export function disposeScene(scene) {
  scene.traverse((obj) => {
    if (obj.geometry) obj.geometry.dispose();
    if (obj.material) {
      if (Array.isArray(obj.material)) obj.material.forEach(m => m.dispose());
      else obj.material.dispose();
    }
  });
}
