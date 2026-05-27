// deck-orchestrator.js — listens to slidechange and orchestrates per-slide init/dispose
const ORCHESTRATOR_VERSION = '20260527-34';  // cache-bust dynamic imports + iframes
const deckStage = document.querySelector('deck-stage');
if (!deckStage) {
  console.warn('[orchestrator] No <deck-stage> encontrado.');
}

const activeContext = { dispose: null };

function disposeActive() {
  if (typeof activeContext.dispose === 'function') {
    try { activeContext.dispose(); } catch (e) { console.warn('[orchestrator] dispose error', e); }
  }
  activeContext.dispose = null;
}

async function initSlide(slide) {
  if (!slide) return;
  const scene = slide.dataset.scene;
  const embed = slide.dataset.embed;
  const video = slide.dataset.video;
  const introVideo = slide.dataset.introVideo;

  let sceneDispose = null;

  if (scene && scene.startsWith('threejs:')) {
    const sceneId = scene.split(':')[1];
    try {
      const module = await import(`./scenes/${sceneId}.js?v=${ORCHESTRATOR_VERSION}`);
      sceneDispose = module.init(slide);
    } catch (e) {
      console.error(`[orchestrator] Escena ${sceneId} fallo al iniciar`, e);
      showSceneFallback(slide, sceneId);
    }
  } else if (scene === 'omni' && video) {
    initOmniVideo(slide, video);
  } else if (embed) {
    initEmbed(slide, embed);
  }

  // Intro video: se reproduce una vez sobre la escena y hace fade-out
  if (introVideo) {
    const introDispose = playIntroVideo(slide, introVideo);
    const prev = sceneDispose || activeContext.dispose;
    activeContext.dispose = () => {
      try { introDispose && introDispose(); } catch (e) {}
      try { prev && prev(); } catch (e) {}
    };
  } else if (sceneDispose) {
    activeContext.dispose = sceneDispose;
  }
}

function playIntroVideo(slide, videoSrc) {
  const video = document.createElement('video');
  video.src = videoSrc;
  video.autoplay = true;
  video.muted = true;
  video.playsInline = true;
  video.preload = 'auto';
  video.style.cssText = 'position:absolute;inset:0;width:100%;height:100%;object-fit:cover;z-index:50;opacity:1;transition:opacity 1.2s ease;';

  let removed = false;
  const fadeOut = () => {
    if (removed) return;
    video.style.opacity = '0';
    setTimeout(() => { if (!removed) { video.remove(); removed = true; } }, 1300);
  };

  video.onended = fadeOut;
  video.onerror = () => { video.remove(); removed = true; };

  // Fallback temporal: si el vídeo dura más de 9 s o el ended no dispara, fade igualmente
  const safetyTimer = setTimeout(fadeOut, 9000);

  slide.appendChild(video);

  return function disposeIntroVideo() {
    clearTimeout(safetyTimer);
    if (!removed) {
      video.pause();
      video.remove();
      removed = true;
    }
  };
}

function showSceneFallback(slide, sceneId) {
  const img = document.createElement('img');
  img.src = `assets/biorender/fallback-${sceneId}.png`;
  img.style.cssText = 'position:absolute;inset:0;width:100%;height:100%;object-fit:cover;z-index:0;';
  img.onerror = () => { img.style.display = 'none'; };
  slide.appendChild(img);
  activeContext.dispose = () => img.remove();
}

function initOmniVideo(slide, videoSrc) {
  let video = slide.querySelector('video[data-omni]');
  if (!video) {
    video = document.createElement('video');
    video.dataset.omni = '';
    video.src = videoSrc;
    video.autoplay = true;
    video.muted = true;
    video.playsInline = true;
    video.loop = true;
    video.preload = 'auto';
    video.style.cssText = 'position:absolute;inset:0;width:100%;height:100%;object-fit:cover;z-index:0;';
    video.onerror = () => showVideoFallback(slide);
    slide.insertBefore(video, slide.firstChild);
  } else {
    video.currentTime = 0;
    video.play().catch(() => showVideoFallback(slide));
  }
  activeContext.dispose = () => { video.pause(); };
}

function showVideoFallback(slide) {
  const gradient = document.createElement('div');
  gradient.style.cssText = 'position:absolute;inset:0;z-index:0;background: radial-gradient(ellipse at 50% 50%, rgba(78,197,224,0.3) 0%, transparent 60%), #06090f;';
  slide.insertBefore(gradient, slide.firstChild);
}

function initEmbed(slide, embedId) {
  const iframe = document.createElement('iframe');
  iframe.src = `embeds/${embedId}_clean.html?v=${ORCHESTRATOR_VERSION}`;
  iframe.style.cssText = 'width:100%;height:100%;border:0;background:#ffffff;';
  iframe.setAttribute('sandbox', 'allow-scripts allow-same-origin');
  const container = slide.querySelector('[data-embed-target]');
  if (container) {
    container.replaceChildren(iframe);
  }
  activeContext.dispose = () => { iframe.remove(); };

  const timer = setTimeout(() => {
    if (!iframe.contentDocument || iframe.contentDocument.readyState === 'uninitialized') {
      console.warn(`[orchestrator] iframe ${embedId} timeout, mostrando fallback`);
      showEmbedFallback(slide, embedId);
    }
  }, 3000);
  iframe.onload = () => clearTimeout(timer);
}

function showEmbedFallback(slide, embedId) {
  const fallbacks = {
    'pca': 'assets/figuras/Figura_06_PCA_fisicoquimico.png',
    'efp': null
  };
  const src = fallbacks[embedId];
  if (!src) return;
  const img = document.createElement('img');
  img.src = src;
  img.style.cssText = 'width:100%;height:auto;display:block;';
  const container = slide.querySelector('[data-embed-target]');
  if (container) {
    container.replaceChildren(img);
  }
}

if (deckStage) {
  deckStage.addEventListener('slidechange', (e) => {
    disposeActive();
    initSlide(e.detail.slide);
  });

  // Si el slidechange inicial ya disparó antes de registrar el listener (carrera
  // por type="module" diferido), inicializar manualmente el slide activo.
  const tryInitActive = () => {
    const active = deckStage.querySelector('section[data-deck-active]')
                || deckStage.querySelector('section:not([hidden])');
    if (active) {
      disposeActive();
      initSlide(active);
    }
  };
  if (document.readyState === 'complete') {
    tryInitActive();
  } else {
    window.addEventListener('load', tryInitActive, { once: true });
  }
}
