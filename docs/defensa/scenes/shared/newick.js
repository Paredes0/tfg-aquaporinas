// scenes/shared/newick.js
// Newick parser + subfamily propagation utilities, extraidos de la primera
// version de forest.js para reutilizarse desde forest.js (radial) y
// tree-static.js (rectangular estatico en slide 14).
//
// API publica:
//   parseNewick(str) -> rootNode
//   assignSubfamilies(rootNode)       (muta el arbol; coloca subfamily,
//                                       sfSet, branchColor, branchSubfamily)
//   postOrderLeaves(node, out)
//   fetchNewick(url, cache)            -> Promise<string>
//   SUBFAMILIES                        constantes oficiales del TFG
//   SF_BY_NAME                         indice por nombre
//   MUTED                              color gris para clados mezclados
//
// El objetivo es no duplicar 200 LOC en dos escenas con la misma logica.

export const SUBFAMILIES = [
  { name: 'PIP', color: '#E74C3C', count: 32 },
  { name: 'TIP', color: '#3498DB', count: 34 },
  { name: 'NIP', color: '#2ECC71', count: 37 },
  { name: 'SIP', color: '#F39C12', count: 12 },
  { name: 'XIP', color: '#9B59B6', count: 6 }
];
export const SF_BY_NAME = Object.fromEntries(SUBFAMILIES.map(sf => [sf.name, sf]));
export const MUTED = 'rgba(232,240,255,0.22)';
export const MUTED_LIGHT = 'rgba(20,24,42,0.20)'; // para fondos blancos (slide 14)

// ---------------------------------------------------------------------------
// 1 · Parser
// ---------------------------------------------------------------------------
// Newick de IQ-TREE: (...)<support>:<branchLen>. <support> puede ser una
// etiqueta compuesta tipo "98.9/1/100".
export function parseNewick(str) {
  let i = 0;
  const s = str;

  function readToken() {
    const start = i;
    while (i < s.length && !'(),:;'.includes(s[i])) i++;
    return s.slice(start, i);
  }

  function parseNode() {
    const node = { name: '', length: 0, children: [], support: null };

    if (s[i] === '(') {
      i++;
      while (true) {
        node.children.push(parseNode());
        if (s[i] === ',') { i++; continue; }
        if (s[i] === ')') { i++; break; }
        if (i >= s.length) break;
      }
      const tok = readToken();
      if (tok.length) {
        if (/[\/]/.test(tok) || /^[0-9.]+$/.test(tok)) {
          node.support = tok;
        } else {
          node.name = tok;
        }
      }
    } else {
      node.name = readToken();
    }

    if (s[i] === ':') {
      i++;
      const lenTok = readToken();
      const v = parseFloat(lenTok);
      node.length = Number.isFinite(v) ? v : 0;
    }
    return node;
  }

  const root = parseNode();
  if (s[i] === ';') i++;
  return root;
}

// ---------------------------------------------------------------------------
// 2 · Fetch helper (cacheado en modulo)
// ---------------------------------------------------------------------------
let _cached = null;
export async function fetchNewick(url = 'assets/arbol_acuaporinas.treefile') {
  if (_cached) return _cached;
  const res = await fetch(url, { cache: 'force-cache' });
  if (!res.ok) throw new Error('HTTP ' + res.status);
  _cached = (await res.text()).trim();
  return _cached;
}

// ---------------------------------------------------------------------------
// 3 · Detectar subfamilia por nombre de hoja
// ---------------------------------------------------------------------------
export function detectSubfamilyFromName(name) {
  for (const sf of SUBFAMILIES) {
    if (name.includes(sf.name)) return sf.name;
  }
  return null;
}

// ---------------------------------------------------------------------------
// 4 · Post-order leaves
// ---------------------------------------------------------------------------
export function postOrderLeaves(node, out) {
  if (node.children.length === 0) {
    out.push(node);
    return;
  }
  for (const c of node.children) postOrderLeaves(c, out);
}

// ---------------------------------------------------------------------------
// 5 · Propagacion de subfamilias por clado monofiletico + marcado de ramas
// ---------------------------------------------------------------------------
export function assignSubfamilies(root, opts = {}) {
  const mutedColor = opts.muted || MUTED;
  const leaves = [];
  postOrderLeaves(root, leaves);
  leaves.forEach(l => { l.subfamily = detectSubfamilyFromName(l.name); });

  function setSF(node) {
    if (node.children.length === 0) {
      node.sfSet = new Set(node.subfamily ? [node.subfamily] : []);
      return node.sfSet;
    }
    const combined = new Set();
    for (const c of node.children) {
      const childSet = setSF(c);
      childSet.forEach(x => combined.add(x));
    }
    node.sfSet = combined;
    return combined;
  }
  setSF(root);

  function fillFromAncestor(node, currentSF) {
    let inherited = currentSF;
    if (node.sfSet && node.sfSet.size === 1) {
      inherited = node.sfSet.values().next().value;
    }
    if (node.children.length === 0) {
      if (!node.subfamily && inherited) node.subfamily = inherited;
      return;
    }
    for (const c of node.children) fillFromAncestor(c, inherited);
  }
  fillFromAncestor(root, null);

  function markBranchColor(node) {
    if (node.sfSet && node.sfSet.size === 1) {
      const sf = node.sfSet.values().next().value;
      node.branchColor = (SF_BY_NAME[sf] && SF_BY_NAME[sf].color) || mutedColor;
      node.branchSubfamily = sf;
    } else {
      node.branchColor = mutedColor;
      node.branchSubfamily = null;
    }
    for (const c of node.children) markBranchColor(c);
  }
  markBranchColor(root);
}
