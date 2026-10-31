// CGV app: state, loading, WebGL read rendering, Canvas2D for ruler/coverage, UI.
'use strict';

// ======================= Display profile =======================
const DEFAULT_PROFILE = {
  version: 2,
  readHeight: 10,                // px, smooth 1–24
  colorBy: 'pair',               // none | strand | pair | tag
  colorTag: 'HP',
  groupBy: 'none',               // none | strand | pair | tag | firstOfPair
  groupTag: 'HP',
  sortBy: 'start',               // start | base | mapq | strand | insert
  mapqAlpha: { low: 0, high: 30, minAlpha: 0.18 },
  baseQualityShading: true,
  filters: { minMapq: 0, hideDuplicates: true, hideSecondary: true, hideSupplementary: false, hideQcFail: true },
  showSoftClips: true,
  coverage: { height: 56, alleleThreshold: 0.2 },
  trackMinHeight: 240,
  showGenes: true,               // MANE transcripts under the reference (GRCh38 only)
  startReference: 'hg38',        // opened at start: hg38 | hg19 | chm13 | demo | none           // px: with many samples each keeps at least this height; the view scrolls
  downsample: { enabled: true, maxReads: 100, window: 50 }, // like IGV: reads kept per window; coverage counts all
  maxReads: 300000,
  maxLoadMB: 64,                 // windows with more compressed alignment data ask before loading
  hideToolbar: false,
  viewPairs: false,              // show mates joined in one row
  cardPair: true,                // read details: the read and its mate side by side
  pairs: { maxSpan: 10000 },     // pairs spanning more than this (bp) are shown as single reads
  language: 'en',                // en | ru
};
const PROFILE_KEY = 'cgv.profile.v1';
function clone(o) { return JSON.parse(JSON.stringify(o)); }
function mergeProfile(base, over) {
  const out = clone(base);
  for (const k of Object.keys(over || {})) {
    if (!Object.prototype.hasOwnProperty.call(base, k)) continue; // own keys only: "__proto__" or "constructor" from a pasted profile are ignored
    if (base[k] && typeof base[k] === 'object' && !Array.isArray(base[k])) out[k] = Object.assign({}, base[k], over[k] && typeof over[k] === 'object' ? over[k] : {});
    else if (typeof over[k] === typeof base[k]) out[k] = over[k];
  }
  return sanitizeProfile(out, base);
}
// Allowed values for a profile pasted or stored by the user: numbers are clamped to a sane range,
// choices must be one of the listed ones, tags two characters. Anything else keeps the default.
const PROFILE_RULES = {
  readHeight: [1, 24], 'mapqAlpha.low': [0, 255], 'mapqAlpha.high': [0, 255], 'mapqAlpha.minAlpha': [0, 1],
  'filters.minMapq': [0, 255], 'coverage.height': [14, 400], 'coverage.alleleThreshold': [0, 1],
  'downsample.maxReads': [1, 100000], 'downsample.window': [1, 1000000], maxReads: [1000, 5000000],
  maxLoadMB: [1, 4096], trackMinHeight: [80, 2000], startReference: ['hg38', 'hg19', 'chm13', 'demo', 'none'], 'pairs.maxSpan': [0, 10000000], version: [0, 1000],
  colorBy: ['none', 'strand', 'pair', 'tag'], groupBy: ['none', 'strand', 'pair', 'tag', 'firstOfPair'],
  sortBy: ['start', 'base', 'mapq', 'strand', 'insert'], language: ['en', 'ru'],
  colorTag: /^[A-Za-z][A-Za-z0-9]$/, groupTag: /^[A-Za-z][A-Za-z0-9]$/,
};
function sanitizeProfile(p, base) {
  const fix = (obj, def, path) => {
    for (const k of Object.keys(def)) {
      const key = path ? path + '.' + k : k, rule = PROFILE_RULES[key], v = obj[k];
      if (def[k] && typeof def[k] === 'object') { if (!v || typeof v !== 'object' || Array.isArray(v)) obj[k] = clone(def[k]); else fix(v, def[k], key); continue; }
      if (typeof v !== typeof def[k]) { obj[k] = def[k]; continue; }
      if (typeof v === 'number') { if (!Number.isFinite(v)) obj[k] = def[k]; else if (Array.isArray(rule)) obj[k] = Math.max(rule[0], Math.min(rule[1], v)); }
      else if (typeof v === 'string' && rule && (Array.isArray(rule) ? !rule.includes(v) : !rule.test(v))) obj[k] = def[k];
    }
    for (const k of Object.keys(obj)) if (!Object.prototype.hasOwnProperty.call(def, k)) delete obj[k];
  };
  fix(p, base, '');
  return p;
}
let P = DEFAULT_PROFILE;
try { const s = localStorage.getItem(PROFILE_KEY); if (s) P = mergeProfile(DEFAULT_PROFILE, JSON.parse(s)); } catch (e) { P = clone(DEFAULT_PROFILE); }
if (P === DEFAULT_PROFILE) P = clone(DEFAULT_PROFILE);
// profile v1 kept the old 8 px default read height; v2 defaults to 10 px
if (P.version < 2) { if (P.readHeight === 8) P.readHeight = 10; P.version = 2; try { localStorage.setItem(PROFILE_KEY, JSON.stringify(P)); } catch (e) { /* storage unavailable */ } }
function saveProfile() { try { localStorage.setItem(PROFILE_KEY, JSON.stringify(P)); } catch (e) { /* storage unavailable: the profile lives until reload */ } }

// ======================= State =======================
// Samples are tracks, stacked under one ruler and reference. All tracks show the same window.
// Track: { id, bam, name, demo, map (view contig → BAM ref), back (BAM ref → view contig),
//          data { reads, cov, insert, total, sampled }, shown, layout, scroll, msg, vao, fbuf, cbuf, n }
const S = {
  tracks: [], trackSeq: 0, fasta: null,
  selected: null,      // selected: { T, rec, row, gi, pn }
  hl: null,            // { name, id }: the highlighted read of track id, outlined with its mate in every panel
  stats: { reads: 0, loadMs: 0, buildMs: 0, drawMs: 0, bytes: 0 },
};

// ======================= Panels (split screen) =======================
// The view can be split into panels side by side, as in IGV: each panel shows its own region of the
// genome, with every track, its own ruler, loaded window, layout and scrolling. The code below the
// panel machinery works on "the current panel": S.view, S.data, S.msg ... and a track's data,
// layout, scroll and GPU buffers (T.data, T.shown, T.layout, T.scroll, T.vao ...) are the current
// panel's. usePanel(pn) makes a panel current; every entry point (an event, a frame, a load that
// resumes after an await) chooses its panel first.
const PANEL_KEYS = ['view', 'data', 'msg', 'loading', 'loadToken', 'forceLoad', 'sortPos', 'rubber'];
const TRACK_PANEL_KEYS = ['data', 'shown', 'layout', 'scroll', 'vao', 'fbuf', 'cbuf', 'n', 'msg'];
const PANEL_GAP = 6, PANEL_HDR = 20, PANEL_MIN_W = 60;
let panelSeq = 0;
function newPanel(view) {
  return { id: ++panelSeq, view: view || { refId: 0, start: 0, bpPerPx: 1 }, data: null, msg: null, loading: false, loadToken: 0, forceLoad: false,
    sortPos: null, rubber: null, x0: 0, w: 0, ts: new Map(), focusName: null };
}
let CUR = newPanel();
S.panels = [CUR]; S.active = CUR;
for (const k of PANEL_KEYS) Object.defineProperty(S, k, { get: () => CUR[k], set: v => { CUR[k] = v; }, enumerable: true });
function usePanel(pn) { CUR = pn; W = pn.w || W; }
// Runs f(pn) with each panel current in turn, then returns to the panel that was current.
function forEachPanel(f) {
  const prev = CUR;
  for (const pn of S.panels.slice()) { usePanel(pn); f(pn); }
  usePanel(S.panels.includes(prev) ? prev : S.active);
}
// A track's state in the current panel.
function trackState(T) {
  let x = CUR.ts.get(T.id);
  if (!x) CUR.ts.set(T.id, x = { data: null, shown: null, layout: null, scroll: 0, vao: null, fbuf: null, cbuf: null, n: 0, msg: null });
  return x;
}
function makeTrack(o) {
  for (const k of TRACK_PANEL_KEYS) { delete o[k]; Object.defineProperty(o, k, { get() { return trackState(o)[k]; }, set(v) { trackState(o)[k] = v; } }); }
  return o;
}
// As many panels as fit at a useful width, and never fewer than 6.
function maxPanels() { return Math.max(6, Math.floor((FULL_W + PANEL_GAP) / (PANEL_MIN_W + PANEL_GAP))); }
// The chromosome map runs across the top; panel headers (with several panels) come under it.
const IDEO_H = 30;
function ideoH() { return S.contigs ? IDEO_H : 0; }
function panelHdr() { return ideoH() + (S.panels.length > 1 ? PANEL_HDR : 0); }
// Places the panels side by side over the full width. keepSpan: each panel keeps the stretch of
// genome it shows (panels added or closed); otherwise it keeps its scale (the window was resized).
function layoutPanels(keepSpan) {
  const n = S.panels.length, step = (FULL_W + PANEL_GAP) / n;
  const prev = CUR;
  S.panels.forEach((pn, i) => {
    const oldW = pn.w, x0 = Math.round(i * step), x1 = Math.round((i + 1) * step - PANEL_GAP * (n > 1 ? 1 : 0));
    pn.x0 = x0; pn.w = Math.max(1, (i === n - 1 ? FULL_W : x1) - x0);
    if (oldW && S.contigs && oldW !== pn.w) {
      usePanel(pn);
      const v = pn.view, c = v.start + oldW * v.bpPerPx / 2;
      if (keepSpan) v.bpPerPx = oldW * v.bpPerPx / pn.w;
      v.start = c - pn.w * v.bpPerPx / 2; clampView();
    }
  });
  usePanel(S.panels.includes(prev) ? prev : S.active);
}
function panelAtX(x) { for (const pn of S.panels) if (x < pn.x0 + pn.w + PANEL_GAP / 2) return pn; return S.panels[S.panels.length - 1]; }
// The panel under a pointer event, made current, and the event's position inside it.
function local(e) { const pn = panelAtX(e.offsetX); usePanel(pn); return { pn, x: e.offsetX - pn.x0, y: e.offsetY }; }
function freePanel(pn) { const prev = CUR; usePanel(pn); for (const T of S.tracks) freeVao(T); usePanel(S.panels.includes(prev) && prev !== pn ? prev : S.active); }
// Back to one panel (a new genome, a sample opened at its own place).
function resetPanels() {
  for (const pn of S.panels) if (pn !== S.panels[0]) freePanel(pn);
  const keep = S.panels[0];
  S.panels = [keep]; S.active = keep; usePanel(keep); layoutPanels(false);
}
function closePanel(pn) {
  if (S.panels.length < 2) return;
  const i = S.panels.indexOf(pn); if (i < 0) return;
  freePanel(pn); S.panels.splice(i, 1);
  if (S.active === pn) S.active = S.panels[Math.min(i, S.panels.length - 1)];
  if (S.selected && S.selected.pn === pn) S.selected = null;
  if (cardState && cardState.pn === pn) hideCard();
  usePanel(S.active); layoutPanels(true); requestLoad(); draw(); syncLocus();
}
// Panels for a list of loci ({ id, s, e }): one locus keeps the active panel, several replace them.
function setPanels(list) {
  list = list.slice(0, maxPanels());
  const keep = S.active;
  for (const pn of S.panels) if (pn !== keep) freePanel(pn);
  S.panels = [keep];
  for (let i = 1; i < list.length; i++) S.panels.push(newPanel({ refId: list[i].id, start: 0, bpPerPx: 1 }));
  S.active = keep; layoutPanels(false);
  S.panels.forEach((pn, i) => { usePanel(pn); goto(list[i].id, list[i].s, list[i].e); });
  usePanel(keep); syncLocus();
}
// IGV's "View mate region in split screen": a new panel next to this one, centered on the mate, at
// the same scale; the read and its mate are outlined in every panel.
function openMatePanel(rec, T) {
  const src = CUR, v = T.back[rec.mateRef];
  if (!(rec.flag & 1) || rec.flag & 8 || rec.mateRef < 0) return;
  if (v < 0) { toast(t('sample.mateNoContig', { name: T.bam.refs[rec.mateRef].name })); return; }
  if (S.panels.length >= maxPanels()) { toast(t('split.full', { n: S.panels.length })); return; }
  const span = W * S.view.bpPerPx, mc = recTags(rec).MC, mlen = typeof mc === 'string' ? cigarStrRefLen(mc) : rec.end - rec.pos;
  const center = rec.matePos + Math.max(1, mlen) / 2;
  const pn = newPanel({ refId: v, start: 0, bpPerPx: S.view.bpPerPx });
  S.panels.splice(S.panels.indexOf(src) + 1, 0, pn);
  layoutPanels(true);
  usePanel(pn); S.hl = pn.focusName = { name: rec.name, id: T.id };
  goto(v, center - span / 2, center + span / 2);
  S.active = pn; syncLocus();
}

// ======================= DOM =======================
const $ = id => document.getElementById(id);
const stage = $('stage'), glCanvas = $('gl'), ov = $('ov'), ctx = ov.getContext('2d');
let W = 0, H = 0, DPR = 1, FULL_W = 0; // W: width of the current panel; FULL_W: of the whole view

// track geometry (CSS px)
const RULER_H = 26, REF_H = 14, TRACK_HDR = 18;
function covTop() { return panelHdr() + RULER_H + REF_H + geneH() + 4; } // where the tracks begin
// ---- gene track (MANE, GRCh38): a fixed band under the reference sequence ----
const GENE_ROW = 21, GENE_H = 2 * GENE_ROW + 4;
// The bundled genes are GRCh38: shown only when the coordinates are GRCh38 (chr1 of its length).
function genesActive() {
  if (!GENES.ready || !S.contigs) return false;
  if (S._g38 && S._g38[0] === S.contigs) return S._g38[1];
  const i = contigId('chr1'), on = i >= 0 && S.contigs[i].len === 248956422;
  S._g38 = [S.contigs, on]; return on;
}
function geneH() { return P.showGenes && genesActive() ? GENE_H : 0; }
function geneTop() { return panelHdr() + RULER_H + REF_H + 2; }
function groupGap() { return P.groupBy === 'none' ? 0 : 20; }
// Vertical place of track i: the space under the ruler is shared equally. With one track there is
// no header and the layout is the same as before tracks existed.
// When the samples do not fit at their minimum height, they are stacked taller than the view and the
// whole stack scrolls under the fixed ruler (S.vscroll, shared by all panels).
function trackSlot() { const N = Math.max(1, S.tracks.length), avail = H - covTop(); return N > 1 ? Math.max(avail / N, P.trackMinHeight) : avail; }
function tracksHeight() { return trackSlot() * Math.max(1, S.tracks.length); }
function clampVScroll() { S.vscroll = Math.max(0, Math.min(tracksHeight() - (H - covTop()), S.vscroll || 0)); }
function vscrollBy(d) { const v = S.vscroll || 0; S.vscroll = v + d; clampVScroll(); if (S.vscroll !== v) draw(); return S.vscroll - v; }
function trackGeom(i) {
  const N = Math.max(1, S.tracks.length), top = covTop(), slot = trackSlot(), y0 = top + i * slot - (S.vscroll || 0);
  const hdr = N > 1 ? TRACK_HDR : 0;
  const covH = N > 1 ? Math.max(14, Math.min(P.coverage.height, Math.round((slot - hdr) * 0.3))) : P.coverage.height;
  const ct = y0 + hdr, rt = ct + covH + 6;
  return { y0, hdr, covTop: ct, covH, readsTop: rt, bottom: Math.max(rt, y0 + slot - (i < N - 1 ? 4 : 0)), end: y0 + slot };
}
function trackAtY(y) { for (let i = 0; i < S.tracks.length; i++) { const g = trackGeom(i); if (y >= g.y0 && y < g.end) return i; } return -1; }

// ======================= Colors (from CSS tokens) =======================
let COL = {};
function hexToRgb(h) {
  h = h.trim();
  if (h.startsWith('rgb')) { const m = h.match(/[\d.]+/g).map(Number); return [m[0], m[1], m[2]]; }
  if (h[0] === '#') { if (h.length === 4) h = '#' + h[1] + h[1] + h[2] + h[2] + h[3] + h[3]; return [parseInt(h.slice(1, 3), 16), parseInt(h.slice(3, 5), 16), parseInt(h.slice(5, 7), 16)]; }
  return [128, 128, 128];
}
function readColors() {
  const cs = getComputedStyle(document.documentElement);
  const g = n => cs.getPropertyValue(n).trim();
  COL = {
    css: { bg: g('--bg'), fg: g('--fg'), muted: g('--muted'), line: g('--line'), cov: g('--cov'), accent: g('--accent'), panel: g('--panel'), gene: g('--gene'), geneClin: g('--gene-clin') },
    read: hexToRgb(g('--read')), fwd: hexToRgb(g('--read-fwd')), rev: hexToRgb(g('--read-rev')),
    del: hexToRgb(g('--del')), skip: hexToRgb(g('--muted')), ins: hexToRgb(g('--ins')), clip: hexToRgb(g('--clip')),
    big: hexToRgb(g('--pair-large')), small: hexToRgb(g('--pair-small')), rf: hexToRgb(g('--pair-rf')), ff: hexToRgb(g('--pair-ff')), mu: hexToRgb(g('--pair-unmapped')),
    base: { 65: hexToRgb(g('--base-a')), 67: hexToRgb(g('--base-c')), 71: hexToRgb(g('--base-g')), 84: hexToRgb(g('--base-t')), 78: hexToRgb(g('--muted')) },
    cat: [1, 2, 3, 4, 5, 6, 7, 8].map(i => hexToRgb(g('--cat-' + i))),
    link: hexToRgb(g('--muted')), overlap: hexToRgb(g('--fg')),
  };
  COL.baseCss = { 65: g('--base-a'), 67: g('--base-c'), 71: g('--base-g'), 84: g('--base-t'), 78: g('--muted') };
}

// ======================= WebGL =======================
const gl = glCanvas.getContext('webgl2', { premultipliedAlpha: false, antialias: false, alpha: true });
const GL = {};
function initGL() {
  if (!gl) return false;
  const vs = `#version 300 es
  precision highp float;
  layout(location=0) in vec2 a_corner;
  layout(location=1) in vec4 a_x;      // x0, x1, row, grp
  layout(location=2) in vec3 a_y;      // yA, yB, minpx
  layout(location=3) in vec4 a_color;
  uniform float u_viewInt, u_viewFrac, u_bpPerPx, u_rowH, u_gap, u_grpGap, u_scroll, u_top, u_mmFade;
  uniform vec2 u_size;
  out vec4 v_color;
  void main() {
    // minpx < 0: an isolated mismatch that fades when zoomed out
    float minpx = abs(a_y.z);
    float fade = a_y.z < 0.0 ? clamp(u_mmFade / u_bpPerPx, 0.0, 1.0) : 1.0;
    // positions are whole numbers relative to the window start; the view start comes as a whole part
    // and a fraction, so the subtraction stays exact in 32-bit floats even far into a large window
    float px0 = ((a_x.x - u_viewInt) - u_viewFrac) / u_bpPerPx;
    float px1 = ((a_x.y - u_viewInt) - u_viewFrac) / u_bpPerPx;
    float h = max(u_rowH - u_gap, 1.0);
    float y0 = u_top + a_x.z * u_rowH + a_x.w * u_grpGap - u_scroll;
    float ya = y0 + h * a_y.x, yb = y0 + h * a_y.y;
    if (yb - ya < 1.0) yb = ya + 1.0;
    float x, y;
    if (a_y.z >= 499.0) {
      // strand arrow, as in IGV: a corner cut off the read's 3' end in the background color, so the
      // end comes to a point. 500/501: right end (forward read), top/bottom corner; 502/503: left end.
      // x0, x1: the read's extent; the cut is at most half the row high and a quarter of the read long.
      float code = a_y.z - 500.0, right = code < 1.5 ? 1.0 : 0.0, bottom = mod(code, 2.0);
      float w = min(h * 0.45, (px1 - px0) * 0.25);
      float e = right > 0.5 ? px1 : px0, dir = right > 0.5 ? -1.0 : 1.0;
      float edge = bottom > 0.5 ? yb : ya, mid = (ya + yb) * 0.5;
      if (w < 1.5 || h < 4.0) w = 0.0;
      // triangle (e + dir*w, edge), (e, edge), (e, mid); the 4th corner repeats the 3rd
      x = a_corner.y > 0.5 ? e : (a_corner.x < 0.5 ? e + dir * w : e);
      y = a_corner.y > 0.5 ? mid : edge;
    } else {
      if (px1 - px0 < minpx) { float c = (px0 + px1) * 0.5; px0 = c - minpx * 0.5; px1 = c + minpx * 0.5; }
      x = mix(px0, px1, a_corner.x); y = mix(ya, yb, a_corner.y);
    }
    gl_Position = vec4(x / u_size.x * 2.0 - 1.0, 1.0 - y / u_size.y * 2.0, 0.0, 1.0);
    v_color = vec4(a_color.rgb, a_color.a * fade);
  }`;
  const fs = `#version 300 es
  precision mediump float;
  in vec4 v_color; out vec4 o;
  void main() { o = v_color; }`;
  const sh = (t, src) => { const s = gl.createShader(t); gl.shaderSource(s, src); gl.compileShader(s); if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(s)); return s; };
  const prog = gl.createProgram();
  gl.attachShader(prog, sh(gl.VERTEX_SHADER, vs)); gl.attachShader(prog, sh(gl.FRAGMENT_SHADER, fs)); gl.linkProgram(prog);
  if (!gl.getProgramParameter(prog, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(prog));
  GL.prog = prog;
  GL.u = {}; for (const n of ['u_viewInt', 'u_viewFrac', 'u_bpPerPx', 'u_rowH', 'u_gap', 'u_grpGap', 'u_scroll', 'u_top', 'u_size', 'u_mmFade']) GL.u[n] = gl.getUniformLocation(prog, n);
  GL.quad = gl.createBuffer(); gl.bindBuffer(gl.ARRAY_BUFFER, GL.quad);
  gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([0, 0, 1, 0, 0, 1, 1, 1]), gl.STATIC_DRAW);
  Object.assign(GL, makeVao());
  gl.enable(gl.BLEND); gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
  return true;
}
// One vertex array + instance buffers; each track has its own.
function makeVao() {
  const vao = gl.createVertexArray(); gl.bindVertexArray(vao);
  gl.bindBuffer(gl.ARRAY_BUFFER, GL.quad);
  gl.enableVertexAttribArray(0); gl.vertexAttribPointer(0, 2, gl.FLOAT, false, 0, 0);
  const fbuf = gl.createBuffer(), cbuf = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, fbuf);
  gl.enableVertexAttribArray(1); gl.vertexAttribPointer(1, 4, gl.FLOAT, false, 28, 0); gl.vertexAttribDivisor(1, 1);
  gl.enableVertexAttribArray(2); gl.vertexAttribPointer(2, 3, gl.FLOAT, false, 28, 16); gl.vertexAttribDivisor(2, 1);
  gl.bindBuffer(gl.ARRAY_BUFFER, cbuf);
  gl.enableVertexAttribArray(3); gl.vertexAttribPointer(3, 4, gl.UNSIGNED_BYTE, true, 0, 0); gl.vertexAttribDivisor(3, 1);
  gl.bindVertexArray(null);
  return { vao, fbuf, cbuf, n: 0 };
}
function freeVao(o) { if (!GL_OK || !o.vao) return; gl.deleteVertexArray(o.vao); gl.deleteBuffer(o.fbuf); gl.deleteBuffer(o.cbuf); o.vao = null; }
const GL_OK = (() => { try { return initGL(); } catch (e) { console.error(e); return false; } })();
// Which GPU path the browser uses: a software renderer is the usual cause of slow frames.
const GL_RENDERER = (() => { if (!GL_OK) return ''; try { const d = gl.getExtension('WEBGL_debug_renderer_info'); return String(gl.getParameter(d ? d.UNMASKED_RENDERER_WEBGL : gl.RENDERER)); } catch (e) { return ''; } })();
const GL_SOFTWARE = /swiftshader|llvmpipe|softpipe|basic render|software|microsoft basic/i.test(GL_RENDERER);
// Drivers compile shaders lazily on the first draw, which can freeze one frame for hundreds of ms.
// Draw one invisible instance during start-up so that cost is paid before any data arrives.
function warmUpGL() {
  if (!GL_OK) return;
  gl.bindBuffer(gl.ARRAY_BUFFER, GL.fbuf); gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([0, 1, 0, 0, 0, 1, 1]), gl.STATIC_DRAW);
  gl.bindBuffer(gl.ARRAY_BUFFER, GL.cbuf); gl.bufferData(gl.ARRAY_BUFFER, new Uint8Array([0, 0, 0, 0]), gl.STATIC_DRAW);
  gl.useProgram(GL.prog); gl.bindVertexArray(GL.vao);
  gl.uniform2f(GL.u.u_size, 1, 1); gl.uniform1f(GL.u.u_bpPerPx, 1); gl.uniform1f(GL.u.u_rowH, 1);
  gl.drawArraysInstanced(gl.TRIANGLE_STRIP, 0, 4, 1);
  gl.clearColor(0, 0, 0, 0); gl.clear(gl.COLOR_BUFFER_BIT);
}

// ======================= Instance buffer =======================
class InstBuf {
  constructor(cap) { this.f = new Float32Array(cap * 7); this.c = new Uint8Array(cap * 4); this.n = 0; }
  grow() { const f = new Float32Array(this.f.length * 2); f.set(this.f); this.f = f; const c = new Uint8Array(this.c.length * 2); c.set(this.c); this.c = c; }
  add(x0, x1, row, grp, yA, yB, minpx, rgb, a) {
    if ((this.n + 1) * 7 > this.f.length) this.grow();
    const i = this.n * 7, j = this.n * 4;
    this.f[i] = x0; this.f[i + 1] = x1; this.f[i + 2] = row; this.f[i + 3] = grp; this.f[i + 4] = yA; this.f[i + 5] = yB; this.f[i + 6] = minpx;
    this.c[j] = rgb[0]; this.c[j + 1] = rgb[1]; this.c[j + 2] = rgb[2]; this.c[j + 3] = Math.round(a * 255);
    this.n++;
  }
}

// ======================= Derived data =======================
function mapqAlpha(m) {
  const { low, high, minAlpha } = P.mapqAlpha;
  if (high <= low) return m >= high ? 1 : minAlpha;
  const t = Math.min(1, Math.max(0, (m - low) / (high - low)));
  return minAlpha + (1 - minAlpha) * t;
}
function qualAlpha(q) { if (!P.baseQualityShading || q === 255) return 1; if (q < 5) return 0.2; if (q < 20) return 0.2 + 0.8 * (q - 5) / 15; return 1; }
function hashStr(s) { let h = 2166136261; s = String(s); for (let i = 0; i < s.length; i++) { h ^= s.charCodeAt(i); h = Math.imul(h, 16777619); } return h >>> 0; }

function readColor(rec, T) {
  switch (P.colorBy) {
    case 'strand': return rec.flag & F_REV ? COL.rev : COL.fwd;
    case 'tag': { const v = recTags(rec)[P.colorTag]; return v === undefined ? COL.read : COL.cat[hashStr(v) % 8]; }
    case 'pair': {
      const o = pairOrientation(rec), ins = T.data.insert;
      if (o === 'other') return COL.cat[(rec.mateRef >>> 0) % 8];
      if (o === 'mateUnmapped') return COL.mu;
      if (o === 'RF') return COL.rf;
      if (o === 'FF' || o === 'RR') return COL.ff;
      if (o === 'FR') { const t = Math.abs(rec.tlen); if (t > ins.hi) return COL.big; if (t && t < ins.lo) return COL.small; }
      return COL.read;
    }
    default: return COL.read;
  }
}

function groupKey(rec) {
  switch (P.groupBy) {
    case 'strand': return rec.flag & F_REV ? t('g.rev') : t('g.fwd');
    case 'firstOfPair': return rec.flag & F_R1 ? t('g.first') : rec.flag & F_R2 ? t('g.second') : t('g.unpaired');
    case 'pair': { const o = pairOrientation(rec); return { FR: t('g.FR'), RF: 'RF', FF: 'FF', RR: 'RR', other: t('g.other'), mateUnmapped: t('g.mateUnmapped'), unpaired: t('g.unpaired') }[o]; }
    case 'tag': { const v = recTags(rec)[P.groupTag]; return v === undefined ? t('g.tagNone', { tag: P.groupTag }) : P.groupTag + ': ' + v; }
    default: return '';
  }
}

// Reads in a layout item: the item itself, or both mates of a pair unit.
function readsOf(it) { return it.reads || [it]; }
function refAtFactory() { const ref = S.data.ref; return p => ref.at(p); }
// Soft-clipped bases that differ from the reference (all of them where the reference is unknown).
function clipMismatches(rec) {
  if (rec._clipMM) return rec._clipMM;
  const out = [], [L, R] = clipLens(rec), ref = S.data.ref;
  const add = (pos, qi) => { const b = readBaseAt(rec, qi), rb = ref.at(pos); if (b !== 78 && b !== rb) out.push({ pos, b, q: rec.qual[qi] }); };
  for (let i = 0; i < L; i++) add(rec.pos - L + i, i);
  for (let i = 0; i < R; i++) add(rec.end + i, rec.lseq - R + i);
  return (rec._clipMM = out);
}
function feats(rec) { if (!rec.feats) rec.feats = computeFeatures(rec, refAtFactory()); return rec.feats; }

// Rebuilds every panel (a setting changed); rebuildPanel only the current one (its data arrived).
function rebuild() { const t0 = performance.now(); forEachPanel(rebuildPanel); S.stats.buildMs = performance.now() - t0; }
function rebuildPanel() {
  const t0 = performance.now();
  for (const T of S.tracks) rebuildTrack(T);
  S.stats.buildMs = performance.now() - t0;
  S.tracks.forEach(clampTrackScroll);
}
// Rebuild of one track: layout → sort → instances. Filters and downsampling acted while loading.
function rebuildTrack(T) {
  const d = S.data, td = T.data;
  if (!d || !td) { T.shown = null; T.n = 0; return; }
  // reads are already filtered (and downsampled) while loading; coverage came with them
  let reads = td.reads;
  let sampled = reads.length < td.total;
  if (reads.length > P.maxReads) { // safety cap; the choice depends on each read alone, not on its neighbours
    const keep = Math.floor(4294967296 * P.maxReads / reads.length);
    reads = reads.filter(r => hashStr(readKey(r)) < keep); sampled = true;
  }
  const cov = td.cov;
  // Stable layout: rows are kept while nothing that decides placement changes (see ReadLayout).
  const sorting = P.sortBy !== 'start' && S.sortPos != null;
  const sig = JSON.stringify([T.id, S.view.refId, P.groupBy, P.groupTag, P.sortBy, sorting ? S.sortPos : null, P.filters, P.downsample, LANG, P.viewPairs, P.viewPairs ? P.pairs.maxSpan : 0, !!P.showSoftClips]);
  if (!T.layout || T.layout.sig !== sig) T.layout = new ReadLayout(sig, 2);
  // Layout items: reads, or in pair view units holding both mates ({ key, pos, end, reads, mate }).
  // With soft clips shown, items span the clipped bases too, so the layout leaves room for them.
  const clips = !!P.showSoftClips, units = P.viewPairs || clips;
  const items = P.viewPairs ? pairUnits(reads, P.pairs.maxSpan, clips)
    : clips ? reads.map(r => { const [a, b] = readExtent(r, true); return { key: readKey(r), pos: a, end: b, reads: [r], mate: null }; }).sort((x, y) => x.pos - y.pos || (x.key < y.key ? -1 : 1))
    : reads;
  const keyOf = units ? (u => u.key) : readKey, firstRead = units ? (u => u.reads[0]) : (r => r);
  let prio = null;
  if (sorting) {
    const k = sortKeyFactory(S.sortPos), p0 = S.sortPos;
    prio = { pos: p0, key: units ? (u => { const r = u.reads.find(x => x.pos <= p0 && x.end > p0); return r ? k(r) : [7, 0]; }) : k };
  }
  const { groups, totalRows: rowCursor } = T.layout.update(items, it => groupKey(firstRead(it)), prio, keyOf);
  // instances: bodies → lines → features on top
  const ib = new InstBuf(Math.max(1024, reads.length * 3));
  const bodies = [], marks = [], links = [];
  groups.forEach((g, gi) => g.rows.forEach((row, ri) => {
    const R = g.start + ri;
    for (const it of row) {
      for (const rec of readsOf(it)) bodies.push([rec, R, gi]);
      if (it.reads && (it.reads.length === 2 || it.mate)) links.push([it, R, gi]);
    }
  }));
  // pair view: a line between mates (to the mate's place when it is not loaded); overlap of mates shaded
  for (const [u, R, gi] of links) {
    const a = mapqAlpha(Math.min(...u.reads.map(r => r.mapq)));
    if (u.reads.length === 2) {
      const [l, r] = u.reads[0].pos <= u.reads[1].pos ? u.reads : [u.reads[1], u.reads[0]];
      if (r.pos > l.end) ib.add(l.end - d.start, r.pos - d.start, R, gi, 0.47, 0.53, 0, COL.link, a);
    } else if (u.mate) {
      const r = u.reads[0], [m0, m1] = u.mate;
      if (m0 >= r.end) ib.add(r.end - d.start, m0 - d.start, R, gi, 0.47, 0.53, 0, COL.link, a * 0.6);
      else if (m1 <= r.pos) ib.add(m1 - d.start, r.pos - d.start, R, gi, 0.47, 0.53, 0, COL.link, a * 0.6);
    }
  }
  for (const [rec, R, gi] of bodies) {
    const col = readColor(rec, T), a = mapqAlpha(rec.mapq);
    // bodies per aligned block, broken at D/N
    let r = rec.pos, segStart = rec.pos;
    for (const c of rec.cigar) {
      const op = c & 15, n = c >>> 4;
      if (op === 0 || op === 7 || op === 8) r += n;
      else if (op === 2 || op === 3) {
        if (r > segStart) ib.add(segStart - d.start, r - d.start, R, gi, 0, 1, 1, col, a);
        marks.push([op === 2 ? 2 : 3, r, r + n, R, gi, a]);
        r += n; segStart = r;
      }
    }
    if (r > segStart) ib.add(segStart - d.start, r - d.start, R, gi, 0, 1, 1, col, a);
    if (P.showSoftClips) { // clipped bases: the read's own color, a little thinner so the clip boundary shows
      const [L, Rc] = clipLens(rec);
      if (L) ib.add(rec.pos - L - d.start, rec.pos - d.start, R, gi, 0.12, 0.88, 1, col, a);
      if (Rc) ib.add(rec.end - d.start, rec.end + Rc - d.start, R, gi, 0.12, 0.88, 1, col, a);
    }
  }
  for (const [t, x0, x1, R, gi, a] of marks) ib.add(x0 - d.start, x1 - d.start, R, gi, t === 2 ? 0.44 : 0.48, t === 2 ? 0.56 : 0.52, 0, t === 2 ? COL.del : COL.skip, a);
  for (const [rec, R, gi] of bodies) {
    const a = mapqAlpha(rec.mapq);
    for (const f of feats(rec)) {
      if (f.t === 0) ib.add(f.pos - d.start, f.pos + 1 - d.start, R, gi, 0, 1, flaggedPos(cov, d, f.pos) ? 1 : -1, COL.base[f.b] || COL.base[78], a * qualAlpha(f.q));
      else if (f.t === 1) ib.add(f.pos - d.start, f.pos - d.start, R, gi, -0.08, 1.08, 2, COL.ins, Math.max(a, 0.6));
      else if ((f.t === 4 || f.t === 5) && !P.showSoftClips && f.len >= 1) ib.add(f.pos - d.start, f.pos - d.start, R, gi, 0, 1, 1.5, COL.clip, Math.max(a, 0.7));
    }
    if (P.showSoftClips) for (const m of clipMismatches(rec)) ib.add(m.pos - d.start, m.pos + 1 - d.start, R, gi, 0.1, 0.9, -1, COL.base[m.b] || COL.base[78], a * qualAlpha(m.q));
  }
  for (const [u, R, gi] of links) {
    if (u.reads.length !== 2) continue;
    const [p, q] = u.reads, o0 = Math.max(p.pos, q.pos), o1 = Math.min(p.end, q.end);
    if (o1 > o0) ib.add(o0 - d.start, o1 - d.start, R, gi, 0, 1, 1, COL.overlap, 0.22);
  }
  // strand arrows last: the 3' end of every read is cut to a point (seen when the row is tall enough)
  const bgc = hexToRgb(COL.css.bg);
  for (const [rec, R, gi] of bodies) {
    const [e0, e1] = readExtent(rec, P.showSoftClips), c = rec.flag & F_REV ? 502 : 500;
    ib.add(e0 - d.start, e1 - d.start, R, gi, 0, 1, c, bgc, 1);
    ib.add(e0 - d.start, e1 - d.start, R, gi, 0, 1, c + 1, bgc, 1);
  }
  T.shown = { reads, cov, groups, totalRows: rowCursor, sampled, total: td.total };
  if (GL_OK) {
    if (!T.vao) Object.assign(T, makeVao());
    gl.bindBuffer(gl.ARRAY_BUFFER, T.fbuf); gl.bufferData(gl.ARRAY_BUFFER, ib.f.subarray(0, ib.n * 7), gl.STATIC_DRAW);
    gl.bindBuffer(gl.ARRAY_BUFFER, T.cbuf); gl.bufferData(gl.ARRAY_BUFFER, ib.c.subarray(0, ib.n * 4), gl.STATIC_DRAW);
    T.n = ib.n;
  }
}

// A position where non-reference alleles pass the threshold (same rule as the coverage highlight).
const CNT = new Uint32Array(6), CNT2 = new Uint32Array(6);
function flaggedPos(cov, d, pos) {
  if (cov.depthOnly) return false;
  const refc = d.ref.at(pos); if (!refc) return false;
  const c = cov.counts(pos, CNT), t = c[0] + c[1] + c[2] + c[3] + c[4]; if (t < 2) return false;
  const ri = { 65: 0, 67: 1, 71: 2, 84: 3 }[refc]; if (ri === undefined) return false;
  return (t - c[ri]) / t >= P.coverage.alleleThreshold;
}

// Sort key of a read at a position; smaller goes higher.
function sortKeyFactory(pos) {
  const rb = refAtFactory()(pos);
  return rec => {
    switch (P.sortBy) {
      case 'base': { const b = baseAtRefPos(rec, pos); if (!b) return [8, 0]; if (b === 45) return [2, 0]; return b === rb ? [3, 0] : [1, b]; }
      case 'mapq': return [1, -rec.mapq];
      case 'strand': return [1, rec.flag & F_REV ? 1 : 0];
      case 'insert': return [1, -Math.abs(rec.tlen)];
    }
    return [1, 0];
  };
}

// ======================= View =======================
// Coordinates come from the reference when one is loaded, otherwise from the BAM header.
// S.contigs: [{ name, len }]; per track T.map[viewContig] = BAM ref id or -1; T.back[bamRef] = view contig or -1.
function refLen() { return S.contigs ? S.contigs[S.view.refId].len : 0; }
function contigId(name) {
  if (!S.contigs) return -1;
  if (!S._ck) { S._ck = new Map(); S.contigs.forEach((c, i) => { for (const k of contigKeys(c.name)) if (!S._ck.has(k)) S._ck.set(k, i); }); }
  for (const k of contigKeys(name)) if (S._ck.has(k)) return S._ck.get(k);
  return -1;
}
function viewEnd() { return S.view.start + W * S.view.bpPerPx; }
function viewCenter() { return S.view.start + W * S.view.bpPerPx / 2; }
function clampView() {
  const L = refLen(); if (!L || !W) return;
  const minBpp = 1 / 48, maxBpp = Math.max(minBpp, L / W);
  S.view.bpPerPx = Math.min(maxBpp, Math.max(minBpp, S.view.bpPerPx));
  const span = W * S.view.bpPerPx;
  S.view.start = Math.min(L - span * 0.5, Math.max(-span * 0.5, S.view.start));
}
// Without groups the scroll range follows the rows the visible part uses, so a deep spot off screen
// does not leave empty space under a coverage dip. With groups their offsets stay fixed.
function contentHeight(T) {
  const sh = T.shown; if (!sh) return 0;
  if (sh.groups.length === 1) return rowsInUse(sh.groups[0], S.view.start, viewEnd()) * P.readHeight;
  return sh.totalRows * P.readHeight + sh.groups.length * groupGap();
}
function clampTrackScroll(T, i) { const g = trackGeom(i), maxS = Math.max(0, contentHeight(T) - (g.bottom - g.readsTop) + 24); T.scroll = Math.min(maxS, Math.max(0, T.scroll || 0)); }
function clampScroll() { forEachPanel(() => S.tracks.forEach(clampTrackScroll)); }
function goto(refId, start, end) {
  S.view.refId = refId;
  const span = Math.max(20, end - start);
  S.view.bpPerPx = span / Math.max(1, W);
  S.view.start = start;
  clampView(); for (const T of S.tracks) T.scroll = 0; requestLoad(); draw(); syncLocus();
}
function zoomAt(px, factor) {
  const bp = S.view.start + px * S.view.bpPerPx;
  S.view.bpPerPx *= factor; clampView();
  S.view.start = bp - px * S.view.bpPerPx; clampView();
  requestLoad(); draw(); syncLocus();
}
function pan(dpx) { S.view.start += dpx * S.view.bpPerPx; clampView(); requestLoad(); draw(); syncLocus(); }

// ======================= Workers =======================
// Samples are read in Web Workers made from the core script's own text (see the end of core.js):
// several samples inflate in parallel and the page stays responsive. Where workers cannot start
// (a page policy, an old browser) everything runs on the page as before.
const POOL = { workers: [], next: 0, seq: 0, wait: new Map(), ok: null, url: null };
// A worker that does not answer in time is stopped and replaced: a stuck decoder never hangs CGV,
// and the page does not retry the same data itself (that would freeze it instead).
const WORKER_TIMEOUT = { open: 30000, fetch: 60000 };
function makeWorker() {
  const w = new Worker(POOL.url);
  w.onmessage = e => {
    if (e.data.op === 'ref') { answerRef(w, e.data); return; } // a CRAM in the worker needs reference bases
    const x = POOL.wait.get(e.data.id); if (x) { POOL.wait.delete(e.data.id); clearTimeout(x.timer); x.cb(e.data); }
  };
  w.onerror = () => failWorker(w, 'worker');
  return w;
}
// Ends every call waiting on worker w with an error, and puts a fresh worker in its place.
function failWorker(w, code, ms) {
  for (const [id, x] of POOL.wait) if (x.w === w) { POOL.wait.delete(id); clearTimeout(x.timer); x.cb({ ok: false, err: { code, vars: { s: Math.round((ms || 0) / 1000) }, message: 'worker ' + code } }); }
  const i = POOL.workers.indexOf(w); if (i < 0) return;
  try { w.terminate(); } catch (e) { /* already gone */ }
  try { POOL.workers[i] = makeWorker(); } catch (e) { POOL.workers.splice(i, 1); }
  // tracks that used it open their file again in a new worker before their next read
  for (const T of S.tracks) if (T.worker === w) { T.worker = null; T.reattach = true; }
}
function workerPool() {
  if (POOL.ok !== null) return POOL.ok && POOL.workers.length > 0;
  try {
    const code = $('cgv-core') && $('cgv-core').textContent;
    if (!code || typeof Worker === 'undefined') throw new Error('no worker');
    POOL.url = URL.createObjectURL(new Blob([code], { type: 'text/javascript' }));
    const n = Math.max(1, Math.min(6, navigator.hardwareConcurrency || 4));
    for (let i = 0; i < n; i++) POOL.workers.push(makeWorker());
    POOL.ok = true;
  } catch (e) { POOL.ok = false; }
  return POOL.ok;
}
function callWorker(w, msg) {
  return new Promise((resolve, reject) => {
    const id = ++POOL.seq;
    const cb = m => m.ok ? resolve(m) : reject(Object.assign(new CgvError(m.err.code || 'worker', m.err.vars), m.err.code ? {} : { message: m.err.message }));
    const ms = WORKER_TIMEOUT[msg.op] || WORKER_TIMEOUT.fetch, timer = setTimeout(() => failWorker(w, 'workerTimeout', ms), ms);
    POOL.wait.set(id, { cb, w, timer });
    w.postMessage(Object.assign({ id }, msg));
  });
}
// Reference bases for a CRAM decoded in a worker, from the page's reference and contig matching.
async function answerRef(w, m) {
  const T = S.tracks.find(x => x.id === m.key);
  let bases = null;
  try { bases = T && T.bam.refProvider ? await T.bam.refProvider(m.refId, m.a, m.b) : null; } catch (e) { bases = null; }
  w.postMessage({ op: 'refResult', rid: m.rid, bases }, bases && bases.buffer && bases.byteOffset === 0 && bases.byteLength === bases.buffer.byteLength ? [bases.buffer] : []);
}
const srcDesc = s => s instanceof UrlSource ? { url: s.url, headers: s.headers } : { blob: s.blob, name: s.name };
// Gives a track its worker; on any failure the track simply reads on the page.
async function attachWorker(T) {
  if (!workerPool()) return;
  const w = POOL.workers[POOL.next++ % POOL.workers.length];
  try { await callWorker(w, { op: 'open', key: T.id, cram: !!T.bam.isCram, bam: srcDesc(T.bam.src), bai: srcDesc(T.bam.baiSrc) }); T.worker = w; }
  catch (e) { T.worker = null; }
}
async function trackFetch(T, br, ls, le, opts) {
  if (T.reattach) { T.reattach = false; await attachWorker(T); }
  if (T.worker) {
    try {
      const m = await callWorker(T.worker, { op: 'fetch', key: T.id, args: [br, ls, le, opts] });
      T.bam.maxSpan = Math.max(T.bam.maxSpan, m.maxSpan || 0);
      const c = m.res.cov, cov = new SparseCoverage(c.start, c.end);
      for (const [bi, blk] of c.blocks) cov.blocks.set(bi, blk);
      if (c.depthOnly) cov.depthOnly = true;
      return { reads: m.res.reads, total: m.res.total, sampled: m.res.sampled, cov };
    } catch (e) { if (e.code !== 'worker') throw e; T.worker = null; } // the worker died: read on the page
  }
  return T.bam.fetch(br, ls, le, opts);
}

// ======================= Loading =======================
let loadTimer = 0;
function requestLoad() { clearTimeout(loadTimer); loadTimer = setTimeout(loadAll, 120); }
// Every panel checks its window; panels whose data still covers their view return at once.
function loadAll() { forEachPanel(() => { ensureData(); }); }
async function ensureData() {
  if (!S.contigs || !W) return;
  const pn = CUR, back = () => { if (!S.panels.includes(pn)) return false; usePanel(pn); return true; }; // after each await
  const vs = Math.max(0, Math.floor(S.view.start)), ve = Math.min(refLen(), Math.ceil(viewEnd()));
  const d = S.data;
  const span0 = ve - vs;
  // a window loaded with depth-only coverage is loaded again in detail when zoomed in below COV_DETAIL_BP
  if (d && d.refId === S.view.refId && vs >= d.start && ve <= d.end && !(d.depthOnly && span0 <= COV_DETAIL_BP)) return;
  const span = span0, margin = Math.round(Math.min(span, 20000) / 2);
  let ls = Math.max(0, vs - margin); const le = Math.min(refLen(), ve + margin);
  const tracks = S.tracks.slice(), refId = S.view.refId;
  // with downsampling, reads starting before the query start are left out (see BamFile.fetch), so the
  // query starts one read length earlier, on a window boundary; reads kept never depend on the query
  if (tracks.length && P.downsample.enabled) {
    const w = P.downsample.window, back = Math.max(1000, ...tracks.map(T => T.bam.maxSpan));
    ls = Math.max(0, Math.floor((ls - back) / w) * w);
  }
  const clear = () => { S.data = null; for (const T of tracks) { T.data = null; T.shown = null; T.n = 0; } };
  if (!tracks.length && span > REF_VIEW_BP) { clear(); S.msg = null; draw(); return; } // reference only: bases are drawn when zoomed in
  // decide by the amount of data the index points to, not by window length: 30x WGS and a
  // 3000x panel differ a hundredfold for the same window
  const est = tracks.reduce((a, T) => a + (T.map[refId] >= 0 ? T.bam.estimateBytes(T.map[refId], ls, le) : 0), 0);
  if (est > P.maxLoadMB * 1048576 && !S.forceLoad) {
    clear();
    const size = fmtBytes(est);
    S.msg = { text: t('wide', { span: fmtBp(span), size }), action: t('wide.btn', { size }) };
    draw(); return;
  }
  S.forceLoad = false;
  const token = ++S.loadToken; S.loading = true; S.msg = null; draw();
  const t0 = performance.now();
  try {
    const depthOnly = span > COV_DETAIL_BP;
    const opts = { depthOnly, filters: P.filters, downsample: P.downsample.enabled ? { max: P.downsample.maxReads, window: P.downsample.window } : null };
    // samples load side by side; one that fails shows its error in its own track
    const results = await Promise.all(tracks.map(async T => {
      const br = T.map[refId];
      if (br < 0) return { T, res: null, msg: t('sample.noContig', { name: S.contigs[refId].name }) };
      try { return { T, res: await trackFetch(T, br, ls, le, opts), msg: null }; }
      catch (e) { return { T, res: null, msg: errText(e), error: true }; }
    }));
    if (!back() || token !== S.loadToken) return;
    const all = results.flatMap(r => r.res ? r.res.reads : []);
    // Reference bases only where they are used: under reads (mismatches, allele highlight) and the
    // view itself when it is narrow enough to show bases. A 46 Mb window with a few reads reads a few kb.
    const want = all.map(r => { const [a, b] = readExtent(r, P.showSoftClips); return [Math.max(ls, a), Math.min(le, b)]; }).filter(x => x[1] > x[0]);
    if (span <= REF_VIEW_BP) want.push([ls, le]);
    const parts = mergeIntervals(want, 2000), ref = new SparseRef(), name = S.contigs[refId].name;
    let fromMD = false, ri = 0, sorted = null;
    for (const [a, b] of parts) {
      let seq = null;
      if (S.fasta) { try { seq = await S.fasta.fetch(name, a, b); } catch (e) { seq = null; } if (!back()) return; }
      if (!seq || seq.length !== b - a) { // no reference: rebuild bases from the reads' MD tags
        seq = new Uint8Array(b - a); fromMD = true;
        if (!sorted) sorted = results.length > 1 ? all.slice().sort((x, y) => x.pos - y.pos) : all;
        while (ri < sorted.length && sorted[ri].end <= a) ri++;
        for (let k = ri; k < sorted.length && sorted[k].pos < b; k++) fillRefFromMD(sorted[k], seq, a);
      }
      ref.add(a, seq);
    }
    if (token !== S.loadToken) return;
    S.refFromMD = fromMD && tracks.length > 0;
    S.data = { refId, start: ls, end: le, ref, depthOnly };
    for (const { T, res, msg, error } of results) {
      T.data = { reads: res ? res.reads : [], insert: insertSizeBounds(res ? res.reads : []),
        cov: res ? res.cov : new SparseCoverage(ls, le), total: res ? res.total : 0, sampled: res ? res.sampled : [] };
      T.msg = msg ? { text: msg, error: !!error } : null;
    }
    S.stats.loadMs = performance.now() - t0;
    S.loading = false;
    // the records were replaced: a selection in this panel no longer points at them
    if (S.selected && S.selected.pn === pn) { S.selected = null; if (cardState && cardState.pn === pn) hideCard(); }
    rebuildPanel();
    if (pn.focusName) { focusRead(pn.focusName); pn.focusName = null; }
    draw();
  } catch (e) {
    if (!back() || token !== S.loadToken) return;
    S.loading = false; S.msg = { text: errText(e), error: true }; draw();
  }
}

// Filters and downsampling act while reading, so changing them re-reads the window.
const REF_VIEW_BP = 200000; // windows up to this size also load the reference between reads
const COV_DETAIL_BP = 1000000; // wider views count coverage depth only (no allele colors), which loads several times faster
function reloadData() { if (!S.contigs) return; forEachPanel(() => { S.data = null; for (const T of S.tracks) { T.data = null; T.shown = null; T.n = 0; } }); S.selected = null; hideCard(); loadAll(); }
// Scrolls each track of the current panel so that a read with this name is in sight (a mate just
// opened in its own panel).
function focusRead(f) {
  const name = f.name;
  S.tracks.forEach((T, i) => {
    const sh = T.shown; if (!sh || T.id !== f.id) return;
    const rh = P.readHeight, gg = groupGap(), g = trackGeom(i);
    for (let gi = 0; gi < sh.groups.length; gi++) {
      const G = sh.groups[gi], ri = G.rows.findIndex(row => row.some(it => readsOf(it).some(r => r.name === name)));
      if (ri < 0) continue;
      const y = (G.start + ri) * rh + gi * gg + gg;
      T.scroll = Math.max(0, y - (g.bottom - g.readsTop) / 3); clampTrackScroll(T, i);
      return;
    }
  });
}

// ======================= Rendering =======================
let raf = 0;
function draw() { if (!raf) raf = requestAnimationFrame(() => { raf = 0; render(); }); }
function fmtBp(n) {
  const L = NUM_LOCALE[LANG];
  if (n >= 1e6) return (n / 1e6).toLocaleString(L, { maximumFractionDigits: 2 }) + ' ' + t('unit.mb');
  if (n >= 1e3) return (n / 1e3).toLocaleString(L, { maximumFractionDigits: 1 }) + ' ' + t('unit.kb');
  return Math.round(n) + ' ' + t('unit.bp');
}
function fmtBytes(n) { return n >= 1 << 30 ? (n / (1 << 30)).toLocaleString(NUM_LOCALE[LANG], { maximumFractionDigits: 1 }) + ' GB' : Math.max(1, Math.round(n / (1 << 20))) + ' MB'; }
function fmtInt(n) { return Math.round(n).toLocaleString(NUM_LOCALE[LANG]); }

const frameTimes = [];
function render() {
  const t0 = performance.now(), tp = {};
  ctx.setTransform(DPR, 0, 0, DPR, 0, 0);
  ctx.clearRect(0, 0, FULL_W, H);
  clampVScroll();
  if (GL_OK) { gl.viewport(0, 0, glCanvas.width, glCanvas.height); gl.clearColor(0, 0, 0, 0); gl.clear(gl.COLOR_BUFFER_BIT); }
  const prev = CUR;
  for (const pn of S.panels) { usePanel(pn); renderPanel(pn, tp); }
  usePanel(S.panels.includes(prev) ? prev : S.active);
  drawPanelFrames();
  drawIdeogram();
  drawMessage(); updateVBar();
  S.stats.drawMs = performance.now() - t0;
  S.stats.parts = tp;
  frameTimes.push(S.stats.drawMs); if (frameTimes.length > 15) frameTimes.shift();
  updateStatus();
}
// One panel: WebGL reads inside the panel's columns, then Canvas2D shifted to the panel and clipped.
function renderPanel(pn, tp) {
  let tl = performance.now();
  const add = (k, v) => { tp[k] = (tp[k] || 0) + v; };
  if (!drag) S.tracks.forEach(clampTrackScroll);
  if (GL_OK && S.data && S.tracks.some(T => T.shown && T.n)) {
    const vx = Math.round(pn.x0 * DPR), vw = Math.max(1, Math.round((pn.x0 + pn.w) * DPR) - vx);
    gl.viewport(vx, 0, vw, glCanvas.height);
    gl.useProgram(GL.prog);
    const vrel = S.view.start - S.data.start, vint = Math.floor(vrel);
    gl.uniform1f(GL.u.u_viewInt, vint); gl.uniform1f(GL.u.u_viewFrac, vrel - vint);
    gl.uniform1f(GL.u.u_bpPerPx, S.view.bpPerPx);
    gl.uniform1f(GL.u.u_rowH, P.readHeight);
    gl.uniform1f(GL.u.u_gap, P.readHeight >= 5 ? 1 : 0);
    gl.uniform1f(GL.u.u_grpGap, groupGap());
    gl.uniform2f(GL.u.u_size, vw / DPR, H);
    gl.uniform1f(GL.u.u_mmFade, 1.5);
    gl.enable(gl.SCISSOR_TEST);
    const clipTop = covTop();
    S.tracks.forEach((T, i) => {
      if (!T.shown || !T.n) return;
      const g = trackGeom(i), y0 = Math.max(g.readsTop, clipTop), y1 = Math.min(g.bottom, H);
      if (y1 <= y0) return; // scrolled out of sight
      gl.scissor(vx, Math.max(0, Math.floor((H - y1) * DPR)), vw, Math.max(0, Math.floor((y1 - y0) * DPR)));
      gl.bindVertexArray(T.vao);
      gl.uniform1f(GL.u.u_scroll, T.scroll);
      gl.uniform1f(GL.u.u_top, g.readsTop + groupGap());
      gl.drawArraysInstanced(gl.TRIANGLE_STRIP, 0, 4, T.n);
    });
    gl.disable(gl.SCISSOR_TEST);
  }
  let n = performance.now(); add('reads', n - tl); tl = n;
  ctx.save(); ctx.translate(pn.x0, 0); ctx.beginPath(); ctx.rect(0, 0, pn.w, H); ctx.clip();
  const part = (k, f) => { f(); const m = performance.now(); add(k, m - tl); tl = m; };
  part('ruler', () => { ctx.save(); ctx.translate(0, panelHdr()); drawRuler(); drawRefTrack(); ctx.restore(); drawGenes(pn); });
  drawRubber();
  // tracks: below the ruler only, and only those in sight
  ctx.save(); ctx.beginPath(); ctx.rect(0, covTop(), pn.w, H - covTop()); ctx.clip();
  const inSight = i => { const g = trackGeom(i); return g.end > covTop() && g.y0 < H; };
  part('coverage', () => S.tracks.forEach((T, i) => { if (inSight(i)) drawCoverage(T, trackGeom(i)); }));
  part('labels', () => S.tracks.forEach((T, i) => { if (!inSight(i)) return; const g = trackGeom(i); drawTrackHeader(T, i, g); drawReadsOverlay(T, g); }));
  ctx.restore();
  ctx.restore();
}
// With several panels: a header over each (its locus and ×, the active one underlined) and a
// line between panels.
function drawPanelFrames() {
  if (S.panels.length < 2 || !S.contigs) return;
  ctx.save();
  const prev = CUR;
  for (const pn of S.panels) {
    usePanel(pn);
    const x = pn.x0, w = pn.w, act = pn === S.active, y0 = ideoH();
    ctx.save(); ctx.translate(0, y0);
    ctx.fillStyle = COL.css.panel; ctx.fillRect(x, 0, w, PANEL_HDR - 1);
    if (act) { ctx.fillStyle = COL.css.accent; ctx.fillRect(x, PANEL_HDR - 3, w, 2); }
    ctx.font = (act ? '600 ' : '') + '11px ' + SANS; ctx.textAlign = 'center'; ctx.textBaseline = 'middle'; ctx.fillStyle = act ? COL.css.fg : COL.css.muted;
    const num = S.panels.indexOf(pn) + 1 + '  '; // the same number marks the panel on the chromosome map
    let label = num + locusText(); const maxW = w - 30;
    if (ctx.measureText(label).width > maxW) label = num + S.contigs[S.view.refId].name + ':' + fmtInt(Math.round(viewCenter()) + 1);
    while (label.length > 2 && ctx.measureText(label).width > maxW) label = label.slice(0, -2);
    ctx.fillText(label, x + w / 2 - 6, PANEL_HDR / 2);
    ctx.font = '13px ' + SANS; ctx.fillStyle = COL.css.muted; ctx.fillText('×', x + w - 10, PANEL_HDR / 2);
    if (pn !== S.panels[0]) { ctx.fillStyle = COL.css.fg; ctx.globalAlpha = 0.3; ctx.fillRect(x - PANEL_GAP / 2 - 0.5, 0, 1, H); ctx.globalAlpha = 1; }
    ctx.restore();
  }
  usePanel(S.panels.includes(prev) ? prev : S.active);
  ctx.textBaseline = 'alphabetic';
  ctx.restore();
}

function niceStep(raw) { const p = Math.pow(10, Math.floor(Math.log10(raw))); for (const m of [1, 2, 5, 10]) if (m * p >= raw) return m * p; return 10 * p; }

function drawRuler() {
  if (!S.contigs) return;
  const bpp = S.view.bpPerPx, step = niceStep(bpp * 110);
  ctx.fillStyle = COL.css.muted; ctx.strokeStyle = COL.css.line; ctx.lineWidth = 1;
  ctx.font = '11px ' + MONO; ctx.textBaseline = 'alphabetic'; ctx.textAlign = 'center';
  const first = Math.ceil(S.view.start / step) * step;
  ctx.beginPath();
  for (let b = first; b < viewEnd(); b += step) {
    const x = Math.round((b - S.view.start) / bpp) + 0.5;
    ctx.moveTo(x, RULER_H - 7); ctx.lineTo(x, RULER_H);
    const lbl = fmtInt(b + 1), tw = ctx.measureText(lbl).width / 2;
    ctx.fillText(lbl, Math.min(W - tw - 2, Math.max(tw + 2, x)), RULER_H - 10);
  }
  ctx.moveTo(0, RULER_H + 0.5); ctx.lineTo(W, RULER_H + 0.5);
  ctx.stroke();
  // sort position
  if (S.sortPos != null && P.sortBy !== 'start') {
    const x0 = (S.sortPos - S.view.start) / bpp, x1 = (S.sortPos + 1 - S.view.start) / bpp;
    ctx.strokeStyle = COL.css.accent; ctx.setLineDash([3, 3]);
    ctx.beginPath(); ctx.moveTo(Math.round(x0) + 0.5, RULER_H); ctx.lineTo(Math.round(x0) + 0.5, H);
    if (x1 - x0 > 3) { ctx.moveTo(Math.round(x1) - 0.5, RULER_H); ctx.lineTo(Math.round(x1) - 0.5, H); }
    ctx.stroke(); ctx.setLineDash([]);
  }
}

// Region being selected on the ruler: shaded band with its length.
function drawRubber() {
  const r = S.rubber; if (!r || !r.moved) return;
  const x0 = Math.min(r.x0, r.x1), x1 = Math.max(r.x0, r.x1);
  ctx.save();
  ctx.fillStyle = COL.css.accent; ctx.globalAlpha = 0.14; ctx.fillRect(x0, 0, x1 - x0, H);
  ctx.globalAlpha = 1; ctx.strokeStyle = COL.css.accent; ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(Math.round(x0) + 0.5, 0); ctx.lineTo(Math.round(x0) + 0.5, H); ctx.moveTo(Math.round(x1) - 0.5, 0); ctx.lineTo(Math.round(x1) - 0.5, H); ctx.stroke();
  const a = Math.floor(S.view.start + x0 * S.view.bpPerPx), b = Math.ceil(S.view.start + x1 * S.view.bpPerPx);
  const label = `${fmtInt(Math.max(1, a + 1))}–${fmtInt(b)} · ${fmtBp(b - a)}`;
  ctx.font = '600 11px ' + SANS; const tw = ctx.measureText(label).width + 12;
  const lx = Math.max(2, Math.min(W - tw - 2, (x0 + x1) / 2 - tw / 2));
  ctx.fillStyle = COL.css.accent; ctx.fillRect(lx, covTop() + 2, tw, 18);
  ctx.fillStyle = COL.css.bg; ctx.textAlign = 'left'; ctx.textBaseline = 'middle'; ctx.fillText(label, lx + 6, covTop() + 11.5);
  ctx.textBaseline = 'alphabetic';
  ctx.restore();
}

function drawRefTrack() {
  const d = S.data; if (!d) return;
  const bpp = S.view.bpPerPx, ppb = 1 / bpp, y = RULER_H + 2;
  if (ppb < 1.2) return;
  const b0 = Math.max(d.start, Math.floor(S.view.start)), b1 = Math.min(d.end, Math.ceil(viewEnd()));
  ctx.textAlign = 'center'; ctx.font = '600 11px ' + MONO;
  for (let b = b0; b < b1; b++) {
    const c = d.ref.at(b); if (!c) continue;
    const x = (b - S.view.start) * ppb;
    if (ppb >= 8) { ctx.fillStyle = COL.baseCss[c] || COL.css.muted; ctx.fillText(String.fromCharCode(c), x + ppb / 2, y + 10); }
    else { ctx.fillStyle = COL.baseCss[c] || COL.css.muted; ctx.fillRect(x, y + 2, Math.max(1, ppb - (ppb > 3 ? 1 : 0)), REF_H - 6); }
  }
}

function drawCoverage(T, g) {
  const sh = T.shown, d = S.data; if (!sh || !d) return;
  const top = g.covTop, h = g.covH, bpp = S.view.bpPerPx, ppb = 1 / bpp;
  const cov = sh.cov;
  const b0 = Math.max(d.start, Math.floor(S.view.start)), b1 = Math.min(d.end, Math.ceil(viewEnd()));
  let maxD = 1;
  cov.forEach(b0, b1, (pos, blk, o) => { const t = blk[o] + blk[o + 1] + blk[o + 2] + blk[o + 3] + blk[o + 4] + blk[o + 5]; if (t > maxD) maxD = t; });
  const scale = Math.max(1, h - 10) / maxD, base = top + h;
  const drawStack = (x, w, pos) => {
    const c = cov.counts(pos, CNT2); let y = base;
    for (const [k, ch] of [[0, 65], [1, 67], [2, 71], [3, 84]]) { const hh = c[k] * scale; if (hh > 0) { ctx.fillStyle = COL.baseCss[ch]; ctx.fillRect(x, y - hh, w, hh); y -= hh; } }
  };
  if (ppb >= 1) {
    cov.forEach(b0, b1, (pos, blk, o) => {
      const x = (pos - S.view.start) * ppb, w = Math.max(1, ppb - (ppb > 4 ? 1 : 0));
      if (flaggedPos(cov, d, pos)) drawStack(x, w, pos);
      else { const t = blk[o] + blk[o + 1] + blk[o + 2] + blk[o + 3] + blk[o + 4]; if (t) { ctx.fillStyle = COL.css.cov; ctx.fillRect(x, base - t * scale, w, t * scale); } }
    });
  } else {
    for (let px = 0; px < W; px++) {
      const s0 = Math.floor(S.view.start + px * bpp), e0 = Math.max(Math.floor(S.view.start + (px + 1) * bpp), s0 + 1);
      if (!cov.hasData(s0, e0)) continue; // empty stretches of a wide window cost nothing
      let mx = 0, fpos = -1;
      cov.forEach(s0, e0, (pos, blk, o) => { const t = blk[o] + blk[o + 1] + blk[o + 2] + blk[o + 3] + blk[o + 4]; if (t > mx) mx = t; if (fpos < 0 && t && flaggedPos(cov, d, pos)) fpos = pos; });
      if (fpos >= 0) drawStack(px, 1, fpos); else if (mx) { ctx.fillStyle = COL.css.cov; ctx.fillRect(px, base - mx * scale, 1, mx * scale); }
    }
  }
  ctx.fillStyle = COL.css.muted; ctx.font = '10.5px ' + MONO; ctx.textAlign = 'left';
  ctx.fillText('0–' + fmtInt(maxD) + '×', 6, top + 10);
  ctx.strokeStyle = COL.css.line; ctx.beginPath(); ctx.moveTo(0, base + 3.5); ctx.lineTo(W, base + 3.5); ctx.stroke();
}

// Track header (only with several tracks): file name, samples, close button; a line above it.
function drawTrackHeader(T, i, g) {
  if (!g.hdr) return;
  ctx.save();
  if (i > 0) { ctx.strokeStyle = COL.css.fg; ctx.globalAlpha = 0.35; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(0, Math.round(g.y0) - 1.5); ctx.lineTo(W, Math.round(g.y0) - 1.5); ctx.stroke(); ctx.globalAlpha = 1; }
  ctx.font = '600 11px ' + SANS; ctx.textAlign = 'left'; ctx.textBaseline = 'middle'; ctx.fillStyle = COL.css.fg;
  const y = g.y0 + g.hdr / 2 + 1;
  let label = T.name + (T.bam.samples.length ? '  ·  ' + T.bam.samples.join(', ') : '');
  const maxW = W - 48;
  while (label.length > 4 && ctx.measureText(label).width > maxW) label = label.slice(0, -2);
  ctx.fillText(label, 8, y);
  if (CUR === S.panels[S.panels.length - 1]) { ctx.fillStyle = COL.css.muted; ctx.font = '13px ' + SANS; ctx.textAlign = 'center'; ctx.fillText('×', W - 14, y); }
  ctx.restore();
}
function drawReadsOverlay(T, g) {
  const sh = T.shown, rt = g.readsTop;
  if (T.msg) { ctx.save(); ctx.font = '12px ' + SANS; ctx.textAlign = 'left'; ctx.textBaseline = 'top'; ctx.fillStyle = T.msg.error ? COL.css.accent : COL.css.muted; ctx.fillText(T.msg.text, 8, rt + 6); ctx.restore(); }
  if (!sh) return;
  ctx.save(); ctx.beginPath(); ctx.rect(0, rt, W, g.bottom - rt); ctx.clip();
  const rh = P.readHeight, gg = groupGap(), bpp = S.view.bpPerPx, ppb = 1 / bpp, d = T.data, scroll = T.scroll, Hh = g.bottom;
  // group labels
  if (gg) sh.groups.forEach((g, gi) => {
    const y = rt + g.start * rh + gi * gg - scroll + gg;
    ctx.strokeStyle = COL.css.line; ctx.beginPath(); ctx.moveTo(0, Math.round(y - gg + 3) + 0.5); ctx.lineTo(W, Math.round(y - gg + 3) + 0.5); ctx.stroke();
    ctx.font = '600 11px ' + SANS; ctx.textAlign = 'left'; ctx.fillStyle = COL.css.fg;
    ctx.fillText(g.key + '  ·  ' + fmtInt(g.count), 8, y - 6);
  });
  // mismatch letters and insertion lengths when zoomed in
  if (ppb >= 7 && rh >= 9) {
    ctx.font = '600 ' + Math.min(12, rh - 1) + 'px ' + MONO; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
    const vs = S.view.start, ve = viewEnd();
    sh.groups.forEach((g, gi) => {
      const firstRow = Math.max(0, Math.floor((scroll - gi * gg) / rh) - g.start - 1);
      const lastRow = Math.min(g.rows.length - 1, Math.ceil((scroll + Hh - rt - gi * gg) / rh) - g.start + 1);
      for (let ri = firstRow; ri <= lastRow; ri++) {
        const y = rt + gg + (g.start + ri) * rh + gi * gg - scroll + (rh - (rh >= 5 ? 1 : 0)) / 2;
        for (const rec of g.rows[ri].flatMap(readsOf)) {
          if (rec.end < vs || rec.pos > ve) continue;
          const a = mapqAlpha(rec.mapq);
          if (P.showSoftClips) for (const m of clipMismatches(rec)) { ctx.globalAlpha = a * qualAlpha(m.q); ctx.fillStyle = '#fff'; ctx.fillText(String.fromCharCode(m.b), (m.pos - vs + 0.5) * ppb, y + 0.5); }
          for (const f of feats(rec)) {
            if (f.t === 0) { ctx.globalAlpha = a * qualAlpha(f.q); ctx.fillStyle = '#fff'; ctx.fillText(String.fromCharCode(f.b), (f.pos - vs + 0.5) * ppb, y + 0.5); }
            else if (f.t === 1 && ppb >= 10) { ctx.globalAlpha = 1; ctx.fillStyle = COL.css.fg; ctx.font = '600 9px ' + MONO; ctx.fillText(String(f.len), (f.pos - vs) * ppb, y - rh / 2 - 1); ctx.font = '600 ' + Math.min(12, rh - 1) + 'px ' + MONO; }
          }
        }
      }
    });
    ctx.globalAlpha = 1; ctx.textBaseline = 'alphabetic';
  }
  // where reads were left out by downsampling (IGV draws the same marks)
  if (d.sampled && d.sampled.length) {
    ctx.fillStyle = COL.css.fg; ctx.globalAlpha = 0.55;
    for (const [a, b] of d.sampled) {
      if (b < S.view.start || a > viewEnd()) continue;
      const x0 = (a - S.view.start) * ppb, x1 = (b - S.view.start) * ppb;
      ctx.fillRect(x0 + (x1 - x0 > 3 ? 0.5 : 0), rt + 0, Math.max(1, x1 - x0 - (x1 - x0 > 3 ? 1 : 0)), 3);
    }
    ctx.globalAlpha = 1;
  }
  // the selected read and its mate (and other alignments of the same read), in every panel
  if (S.hl && S.hl.id === T.id) {
    const vs = S.view.start, ve = viewEnd();
    ctx.strokeStyle = COL.css.accent; ctx.lineWidth = 2;
    sh.groups.forEach((g, gi) => {
      const firstRow = Math.max(0, Math.floor((scroll - gi * gg) / rh) - g.start - 1);
      const lastRow = Math.min(g.rows.length - 1, Math.ceil((scroll + Hh - rt - gi * gg) / rh) - g.start + 1);
      for (let ri = firstRow; ri <= lastRow; ri++) {
        const y = rt + gg + (g.start + ri) * rh + gi * gg - scroll;
        for (const rec of g.rows[ri].flatMap(readsOf)) {
          if (rec.name !== S.hl.name || rec.end < vs || rec.pos > ve) continue;
          const [e0, e1] = readExtent(rec, P.showSoftClips), x0 = (e0 - vs) * ppb, x1 = (e1 - vs) * ppb;
          ctx.strokeRect(x0 - 1, y - 1, x1 - x0 + 2, Math.max(2, rh - (rh >= 5 ? 1 : 0)) + 2);
        }
      }
    });
    ctx.lineWidth = 1;
  }
  ctx.restore();
}

function drawMessage() {
  const el = $('msg');
  // the active panel's message first, then any other panel's (with its locus when there are several)
  const pn = [S.active, ...S.panels].find(p => p.msg) || null;
  if (S.panels.some(p => p.loading)) { el.hidden = false; el.className = 'msg'; el.innerHTML = '<span class="spin"></span>' + esc(t('loading')); return; }
  if (!S.contigs) { el.hidden = true; return; }
  if (pn) {
    el.hidden = false; el.className = 'msg' + (pn.msg.error ? ' err' : '');
    el.textContent = (S.panels.length > 1 ? panelLocus(pn) + ': ' : '') + pn.msg.text;
    if (pn.msg.action) { const b = document.createElement('button'); b.className = 'btn'; b.textContent = pn.msg.action; b.onclick = () => { if (!S.panels.includes(pn)) return; usePanel(pn); S.forceLoad = true; S.msg = null; ensureData(); }; el.append(' ', b); }
    return;
  }
  el.hidden = true;
}

function updateStatus() {
  const st = S.stats;
  let sh = null; st.reads = 0;
  forEachPanel(() => S.tracks.forEach(T => {
    if (T.data) st.reads += T.data.total;
    if (T.shown) { sh = sh || { sampled: false, n: 0 }; sh.sampled = sh.sampled || T.shown.sampled; sh.n += T.shown.reads.length; }
  }));
  let s = '';
  // median of recent frames: one slow frame (first draw, a busy moment) does not stick in the status
  const med = frameTimes.length ? [...frameTimes].sort((a, b) => a - b)[frameTimes.length >> 1] : 0;
  if (sh) s = t('status', { reads: fmtInt(st.reads), sampled: sh.sampled ? t('status.sampled', { n: fmtInt(sh.n) }) : '', load: Math.round(st.loadMs), build: Math.round(st.buildMs), frame: med.toFixed(1) });
  if (sh && med > 40 && st.parts) s += ' (' + Object.entries(st.parts).map(([k, v]) => t('part.' + k) + ' ' + Math.round(v)).join(', ') + ')';
  if (S.refFromMD && sh) s += ' · ' + t('ref.fromMd');
  const el = $('status'); el.textContent = s;
  el.title = 'WebGL: ' + (GL_RENDERER || '—') + (GL_SOFTWARE ? ' · ' + t('gl.software') : '');
}

// ======================= Locus =======================
// Locus of the current panel, as typed in the field.
function locusText() {
  const s = Math.max(0, Math.round(S.view.start)), e = Math.min(refLen(), Math.round(viewEnd()));
  return `${S.contigs[S.view.refId].name}:${fmtInt(s + 1)}-${fmtInt(e)}`;
}
function panelLocus(pn) { const prev = CUR; usePanel(pn); const s = locusText(); usePanel(prev); return s; }
// The field lists every panel's locus, as IGV does; typing several loci opens one panel for each.
function syncLocus() {
  if (!S.contigs || document.activeElement === $('locus')) return;
  $('locus').value = S.panels.map(panelLocus).join('  ');
  const prev = CUR; usePanel(S.active);
  $('span').textContent = fmtBp(Math.min(refLen(), viewEnd()) - Math.max(0, S.view.start));
  usePanel(prev);
}
function parseLocus(str) {
  if (typeof str !== 'string' || str.length > 300) return null;
  str = str.trim().replace(/[  \s]/g, '').replace(/[–—]/g, '-');
  const m = str.match(/^([^:]+?)(?::([\d,.]+)(?:-([\d,.]+))?)?$/);
  if (!m) return null;
  const id = contigId(m[1]);
  if (id < 0) {
    // not a contig: a gene symbol or a MANE transcript id (ENST…, NM_…)
    if (!m[2] && GENES.ready) {
      const tx = findGene(m[1])[0];
      if (tx && !genesActive()) return { error: t('locus.genesGrch38') };
      const gid = tx ? contigId(tx.chr) : -1;
      if (gid >= 0) { const pad = Math.max(50, Math.round((tx.end - tx.start) * 0.05)); return { id: gid, s: Math.max(0, tx.start - pad), e: Math.min(S.contigs[gid].len, tx.end + pad), gene: tx }; }
    }
    return { error: t(GENES.ready ? 'locus.noContigGene' : 'locus.noContig', { name: m[1] }) };
  }
  const num = x => parseInt(x.replace(/[,.]/g, ''), 10);
  const L = S.contigs[id].len;
  if (!m[2]) return { id, s: 0, e: L };
  // positions are 1-based in the text; out-of-range values are clamped to the contig
  const clamp = v => Math.max(0, Math.min(L, v));
  const a = num(m[2]) - 1;
  if (!Number.isSafeInteger(a)) return null;
  if (!m[3]) { const c = clamp(a); return { id, s: c - 50, e: c + 51, center: c }; }
  const b = num(m[3]);
  if (!Number.isSafeInteger(b)) return null;
  // "500-100" means the same region as "100-500"
  const s = clamp(Math.min(a, b - 1)), e = clamp(Math.max(a + 1, b));
  return { id, s, e: Math.max(e, s + 1) };
}

// One locus, or several separated by spaces (one panel each). Returns a list, { error } or null.
function parseLoci(str) {
  const one = parseLocus(str);
  if (one && !one.error) return [one];
  const toks = String(str).trim().split(/[\s;]+/).filter(Boolean);
  if (toks.length < 2) return one;
  const list = [];
  for (const x of toks) { const l = parseLocus(x); if (!l || l.error) return l; list.push(l); }
  return list;
}

// ======================= Input =======================
let drag = null;
// Genome gestures only start on the canvas itself: panels, cards and messages on top of
// the stage keep their own clicks, scrolling and text selection.
const onCanvas = e => e.target === ov;
// The ruler band (coordinates + reference bases): drag selects a region to zoom into, as in IGV;
// a plain click centers the view there.
const inRuler = y => y >= panelHdr() && y < geneTop();
const inGenes = y => geneH() > 0 && y >= geneTop() && y < geneTop() + GENE_H;
const clampX = x => Math.max(0, Math.min(W, x));
function setActive(pn) { if (S.active !== pn) { S.active = pn; syncLocus(); draw(); } }
let ideoDrag = false;
stage.addEventListener('pointerdown', e => {
  if (e.button !== 0 || !S.contigs || !onCanvas(e)) return;
  if (e.offsetY < ideoH()) { stage.setPointerCapture(e.pointerId); ideoDrag = true; ideoGo(e.offsetX); return; } // the chromosome map
  const L = local(e), pn = L.pn;
  setActive(pn);
  if (L.y < panelHdr()) { if (L.x > W - 22) closePanel(pn); return; } // a panel's header: × closes it
  stage.setPointerCapture(e.pointerId);
  if (inRuler(L.y)) { S.rubber = { x0: clampX(L.x), x1: clampX(L.x), moved: false }; draw(); return; }
  const ti = trackAtY(L.y), T = S.tracks[ti];
  drag = { pn, x: e.clientX, y: e.clientY, sx: S.view.start, T, ti, sc: T ? T.scroll : 0, vs: S.vscroll || 0, moved: false, inReads: !!T && L.y > trackGeom(ti).readsTop };
});
stage.addEventListener('pointermove', e => {
  if (ideoDrag) { ideoGo(e.offsetX); return; }
  const rp = S.panels.find(p => p.rubber);
  if (rp) {
    usePanel(rp);
    S.rubber.x1 = clampX(e.offsetX - rp.x0);
    if (Math.abs(S.rubber.x1 - S.rubber.x0) > 3) S.rubber.moved = true;
    draw(); return;
  }
  if (!drag) { if (!S.contigs || !onCanvas(e)) { stage.style.cursor = ''; return; } if (e.offsetY < ideoH()) { stage.style.cursor = 'pointer'; return; } const L = local(e); stage.style.cursor = inRuler(L.y) ? 'crosshair' : L.y < panelHdr() && L.x > W - 22 ? 'pointer' : ''; return; }
  if (!S.panels.includes(drag.pn)) { drag = null; return; }
  usePanel(drag.pn);
  const dx = e.clientX - drag.x, dy = e.clientY - drag.y;
  if (Math.abs(dx) + Math.abs(dy) > 3) { if (!drag.moved) stage.style.cursor = 'grabbing'; drag.moved = true; }
  if (!drag.moved) return;
  S.view.start = drag.sx - dx * S.view.bpPerPx; clampView();
  if (drag.inReads && S.tracks[drag.ti] === drag.T) { drag.T.scroll = drag.sc - dy; clampTrackScroll(drag.T, drag.ti); }
  else if (!drag.inReads) { S.vscroll = drag.vs - dy; clampVScroll(); } // over a coverage track or a header: the whole stack moves
  requestLoad(); draw(); syncLocus();
});
stage.addEventListener('pointercancel', () => { ideoDrag = false; });
stage.addEventListener('pointercancel', () => { for (const p of S.panels) p.rubber = null; drag = null; stage.style.cursor = ''; draw(); });
stage.addEventListener('pointerup', e => {
  if (ideoDrag) { ideoDrag = false; return; }
  const rp = S.panels.find(p => p.rubber);
  if (rp) {
    usePanel(rp);
    const r = S.rubber; S.rubber = null;
    const a = S.view.start + Math.min(r.x0, r.x1) * S.view.bpPerPx, b = S.view.start + Math.max(r.x0, r.x1) * S.view.bpPerPx;
    if (r.moved) goto(S.view.refId, Math.floor(a), Math.ceil(b));
    else { const c = S.view.start + r.x0 * S.view.bpPerPx, span = W * S.view.bpPerPx; goto(S.view.refId, c - span / 2, c + span / 2); }
    return;
  }
  if (!drag) return;
  const d = drag; drag = null; stage.style.cursor = '';
  if (!S.panels.includes(d.pn)) return;
  usePanel(d.pn);
  if (!d.moved) onClick(e.offsetX - d.pn.x0, e.offsetY, e.offsetX);
});
stage.addEventListener('dblclick', e => { if (!S.contigs || !onCanvas(e)) return; const L = local(e); if (L.y >= panelHdr()) zoomAt(L.x, 0.5); });
stage.addEventListener('wheel', e => {
  if (!S.contigs || !onCanvas(e)) return;
  closeMenu();
  e.preventDefault();
  const L = local(e);
  const dy = e.deltaMode === 1 ? e.deltaY * 16 : e.deltaY, dx = e.deltaMode === 1 ? e.deltaX * 16 : e.deltaX;
  if (e.ctrlKey || e.metaKey) { zoomAt(L.x, Math.exp(dy * 0.01)); return; }
  if (e.shiftKey || Math.abs(dx) > Math.abs(dy)) { pan(e.shiftKey ? dy : dx); return; }
  // over the reads: the sample's own rows, and once they reach their end, the whole stack of samples;
  // over coverage or a header: the stack
  const ti = trackAtY(L.y), T = S.tracks[ti];
  if (!T || L.y < trackGeom(ti).readsTop) { vscrollBy(dy); return; }
  const before = T.scroll; T.scroll += dy; clampTrackScroll(T, ti);
  const rest = dy - (T.scroll - before);
  if (Math.abs(rest) >= 1) vscrollBy(rest);
  draw();
}, { passive: false });

function hitTest(x, y) {
  const ti = trackAtY(y), T = S.tracks[ti]; if (!T || !T.shown) return null;
  const sh = T.shown, tg = trackGeom(ti);
  const rt = tg.readsTop, rh = P.readHeight, gg = groupGap();
  if (y < rt || y > tg.bottom) return null;
  for (let gi = 0; gi < sh.groups.length; gi++) {
    const g = sh.groups[gi];
    const gy = rt + gg + g.start * rh + gi * gg - T.scroll;
    const ri = Math.floor((y - gy) / rh);
    if (ri < 0 || ri >= g.rows.length) continue;
    const bp = S.view.start + x * S.view.bpPerPx, slop = 2 * S.view.bpPerPx;
    const row = g.rows[ri];
    let lo = 0, hi = row.length - 1, best = null;
    while (lo <= hi) { const m = (lo + hi) >> 1; if (row[m].pos <= bp + slop) { best = m; lo = m + 1; } else hi = m - 1; }
    for (let k = best; k != null && k >= 0 && k >= best - 3; k--) {
      if (row[k].end + slop <= bp) continue;
      const rs = readsOf(row[k]);
      const ext = r => readExtent(r, P.showSoftClips);
      const hit = rs.find(r => { const [a, b] = ext(r); return a - slop <= bp && b + slop > bp; })
        || rs.reduce((a, r) => { const [ra, rb] = ext(r), [aa, ab] = ext(a); return Math.min(Math.abs(ra - bp), Math.abs(rb - bp)) < Math.min(Math.abs(aa - bp), Math.abs(ab - bp)) ? r : a; });
      return { T, rec: hit, row: g.start + ri, gi };
    }
  }
  return null;
}

// x: inside the current panel; gx: on the whole view (where the card goes).
function onClick(x, y, gx) {
  if (inGenes(y)) { const h = geneHit(x, y); if (h) { S.selected = null; S.hl = null; showGeneCard(h, Math.floor(S.view.start + x * S.view.bpPerPx), gx, y); } else hideCard(); draw(); return; }
  const ti = trackAtY(y);
  if (ti >= 0) { const g = trackGeom(ti); if (g.hdr && y < g.y0 + g.hdr && x > W - 28 && CUR === S.panels[S.panels.length - 1]) { closeTrack(S.tracks[ti]); return; } }
  const h = hitTest(x, y);
  if (!h) { S.selected = null; S.hl = null; hideCard(); draw(); return; }
  S.selected = Object.assign(h, { pn: CUR }); S.hl = { name: h.rec.name, id: h.T.id }; showCard(h.rec, h.T, x, gx, y); draw();
}

const FLAG_BITS = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048];
function esc(s) { return String(s).replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c])); }
// ---- Read details, laid out like IGV's popup: "key = value" lines in blocks ----
function readBaseAt(rec, qi) { const b = rec.seq[qi >> 1]; return '=ACMGRSVTWYHKDBN'.charCodeAt(qi & 1 ? b & 15 : b >> 4); }
function rgSamples(bam) {
  if (bam._rg) return bam._rg;
  const m = new Map();
  for (const line of bam.headerText.split('\n')) if (line.startsWith('@RG')) { const id = line.match(/\tID:([^\t]+)/), sm = line.match(/\tSM:([^\t]+)/); if (id) m.set(id[1], sm ? sm[1] : ''); }
  return (bam._rg = m);
}
// Contig name as shown: the reference's name for contigs matched to it, the BAM's otherwise.
function contigName(bamRef, T) { const v = T.back[bamRef]; return v >= 0 ? S.contigs[v].name : T.bam.refs[bamRef].name; }
// The loaded mate of a read: in the current panel's window, or else in another panel's.
function findMate(rec, T) {
  if (!(rec.flag & 1)) return null;
  const want = rec.flag & 64 ? 128 : 64, look = () => T.data ? T.data.reads.find(r => r !== rec && r.name === rec.name && (r.flag & want) && !(r.flag & 0x900)) || null : null;
  let m = look();
  if (!m) forEachPanel(() => { if (!m) m = look(); });
  return m;
}
// Pair orientation as IGV writes it: strand + read number of the left read, then of the right one.
function pairOrientation2(rec) {
  if (!(rec.flag & 1) || rec.flag & 8 || rec.mateRef !== rec.refID) return '';
  const me = (rec.flag & 16 ? 'R' : 'F') + (rec.flag & 64 ? '1' : '2'), mate = (rec.flag & 32 ? 'R' : 'F') + (rec.flag & 64 ? '2' : '1');
  return rec.pos <= rec.matePos ? me + mate : mate + me;
}
// Lines of one alignment: ['kv', key, value] | ['sep'] | ['txt', text].
function alignmentLines(rec, T) {
  const L = [], kv = (k, v) => L.push(['kv', k, v]), tags = recTags(rec), f = rec.flag;
  const strand = r => (r.flag & 16 ? '(-)' : '(+)');
  kv(t('card.readName'), rec.name);
  if (tags.RG !== undefined) { const sm = rgSamples(T.bam).get(String(tags.RG)); if (sm) kv(t('card.sample'), sm); kv(t('card.rg'), tags.RG); }
  kv(t('card.readLen'), rec.lseq + 'bp');
  kv(t('card.flags'), f);
  L.push(['sep']);
  kv(t('card.mapping'), (f & 256 ? t('card.secondary') : f & 2048 ? t('card.supplementary') : t('card.primary')) + ' @ MAPQ ' + rec.mapq);
  if (f & 1024) kv(t('card.dup'), t('card.yes'));
  if (f & 512) kv(t('card.qcfail'), t('card.yes'));
  kv(t('card.refSpan'), `${contigName(rec.refID, T)}:${fmtInt(rec.pos + 1)}-${fmtInt(rec.end)} ${strand(rec)} = ${rec.end - rec.pos}bp`);
  kv(t('card.cigar'), cigarString(rec));
  const [sl, sr] = clipLens(rec), [hl, hr] = hardClipLens(rec), clip = [];
  if (sl || hl) clip.push(t('card.left') + ' ' + [sl ? sl + ' ' + t('card.soft') : '', hl ? hl + ' ' + t('card.hard') : ''].filter(Boolean).join(', '));
  if (sr || hr) clip.push(t('card.right') + ' ' + [sr ? sr + ' ' + t('card.soft') : '', hr ? hr + ' ' + t('card.hard') : ''].filter(Boolean).join(', '));
  kv(t('card.clipping'), clip.length ? clip.join('; ') : t('card.none'));
  if (f & 1) {
    L.push(['sep']);
    const mapped = !(f & 8) && rec.mateRef >= 0;
    kv(t('card.mateMapped'), mapped ? t('card.yes') : t('card.no'));
    if (mapped) kv(t('card.mateStart'), `${contigName(rec.mateRef, T)}:${fmtInt(rec.matePos + 1)} ${f & 32 ? '(-)' : '(+)'}`);
    kv(t('card.insert'), rec.tlen);
    L.push(['txt', f & 64 ? t('card.first') : f & 128 ? t('card.second') : '']);
    const po = pairOrientation2(rec); if (po) kv(t('card.pairOr'), po);
  }
  L.push(['sep']);
  const hidden = [];
  for (const k of Object.keys(tags)) { if (k === 'RG') { hidden.push(k); continue; } kv(k, Array.isArray(tags[k]) ? tags[k].join(',') : tags[k]); }
  if (hidden.length) L.push(['txt', t('card.hidden') + ': ' + hidden.join(', ')]);
  return L;
}
// The clicked base: Location and "Base = A @ QV 14", as IGV shows them under the clicked read.
function baseLines(rec, pos) {
  const L = [['kv', t('card.location'), `${S.contigs[S.view.refId].name}:${fmtInt(pos + 1)}`]];
  const [cl, cr] = clipLens(rec);
  let qi = null, note = '';
  const at = queryAt(rec, pos);
  if (at && at.del) { L.push(['kv', t('card.base'), t('card.deletion')]); return L; }
  if (at) qi = at.qi;
  else if (P.showSoftClips && pos < rec.pos && pos >= rec.pos - cl) { qi = pos - (rec.pos - cl); note = ' (' + t('card.softclip') + ')'; }
  else if (P.showSoftClips && pos >= rec.end && pos < rec.end + cr) { qi = rec.lseq - cr + (pos - rec.end); note = ' (' + t('card.softclip') + ')'; }
  if (qi != null && rec.lseq) { const q = rec.qual[qi]; L.push(['kv', t('card.base'), String.fromCharCode(readBaseAt(rec, qi)) + (q === 255 ? '' : ' @ QV ' + q) + note]); }
  return L;
}
function linesHtml(lines, cls) {
  return lines.map(l => l[0] === 'sep' ? '<div class="sep"></div>' : l[0] === 'txt' ? (l[1] ? `<div class="ln">${esc(l[1])}</div>` : '')
    : `<div class="ln${cls ? ' ' + cls : ''}">${esc(l[1])} = ${esc(String(l[2]))}</div>`).join('');
}
function linesText(lines) { return lines.map(l => l[0] === 'sep' ? '----------------------' : l[0] === 'txt' ? l[1] : `${l[1]} = ${l[2]}`).filter(x => x !== '').join('\n'); }

let cardState = null;
function showCard(rec, T, x, gx, y) { cardState = { rec, T, x: gx, y, pos: Math.floor(S.view.start + x * S.view.bpPerPx), pn: CUR }; renderCard(true); }
function renderCard(place) {
  if (!S.panels.includes(cardState.pn)) { hideCard(); return; }
  usePanel(cardState.pn);
  const { rec, T, x, y, pos } = cardState, c = $('card');
  const paired = !!(rec.flag & 1) && !(rec.flag & 0x900);
  const mate = paired && P.cardPair ? findMate(rec, T) : null;
  const cols = [];
  if (P.cardPair && paired) {
    const mateRec = mate;
    const left = !mateRec || rec.pos <= mateRec.pos ? rec : mateRec, right = left === rec ? mateRec : rec;
    for (const [r, title] of [[left, t('card.leftAln')], [right, t('card.rightAln')]]) {
      if (!r) { cols.push({ title, lines: [['txt', (rec.flag & 8) ? t('card.mateUnmapped2') : t('card.mateNotLoaded')]], base: null }); continue; }
      cols.push({ title, lines: alignmentLines(r, T), base: r === rec ? baseLines(rec, pos) : null });
    }
  } else cols.push({ title: '', lines: alignmentLines(rec, T), base: baseLines(rec, pos) });
  const seg = paired ? `<div class="seg" role="group" aria-label="${esc(t('card.show'))}"><button class="btn${P.cardPair ? '' : ' on'}" data-pair="0" aria-pressed="${!P.cardPair}">${esc(t('card.showRead'))}</button><button class="btn${P.cardPair ? ' on' : ''}" data-pair="1" aria-pressed="${P.cardPair}">${esc(t('card.showPair'))}</button></div>` : '';
  const goMate = paired && rec.mateRef >= 0 && !(rec.flag & 8);
  c.classList.toggle('wide', cols.length > 1);
  c.innerHTML = seg + `<div class="cols" style="--n:${cols.length}">` + cols.map(col =>
    `<div class="col">${col.title ? `<h3>${esc(col.title)}</h3>` : ''}${linesHtml(col.lines)}${col.base ? '<div class="loc">' + linesHtml(col.base, 'b') + '</div>' : ''}</div>`).join('') + '</div>' +
    `<div class="card-actions">${goMate ? `<button class="btn" id="toMate">${esc(t('card.toMate'))}</button><button class="btn" id="mateSplit">${esc(t('card.mateSplit'))}</button>` : ''}<button class="btn" id="copyRead">${esc(t('card.copy'))}</button><button class="btn ghost" id="closeCard">${esc(t('btn.close'))}</button></div>`;
  c.hidden = false;
  if (place) {
    const r = stage.getBoundingClientRect(), cw = c.offsetWidth, ch = c.offsetHeight;
    c.style.left = Math.max(8, Math.min(r.width - cw - 8, x + 14)) + 'px';
    c.style.top = Math.max(8, Math.min(r.height - ch - 8, y + 14)) + 'px';
  } else { // keep it on screen after a width change
    const r = stage.getBoundingClientRect(), cw = c.offsetWidth;
    c.style.left = Math.max(8, Math.min(r.width - cw - 8, parseFloat(c.style.left) || 8)) + 'px';
  }
  c.querySelectorAll('[data-pair]').forEach(b => b.onclick = () => { P.cardPair = b.dataset.pair === '1'; saveProfile(); renderCard(false); });
  $('closeCard').onclick = () => { S.selected = null; hideCard(); draw(); };
  $('copyRead').onclick = () => {
    const txt = cols.map(col => (col.title ? col.title + '\n' : '') + linesText(col.lines) + (col.base ? '\n' + linesText(col.base) : '')).join('\n\n');
    navigator.clipboard.writeText(txt).then(() => toast(t('copied')), () => toast(t('copyFail')));
  };
  if (goMate) $('mateSplit').onclick = () => { if (!S.panels.includes(cardState.pn)) return; usePanel(cardState.pn); hideCard(); openMatePanel(rec, T); draw(); };
  if (goMate) $('toMate').onclick = () => {
    if (!S.panels.includes(cardState.pn)) return;
    usePanel(cardState.pn); hideCard(); S.selected = null;
    const v = T.back[rec.mateRef];
    if (v < 0) { toast(t('sample.mateNoContig', { name: T.bam.refs[rec.mateRef].name })); draw(); return; }
    goto(v, rec.matePos - 150, rec.matePos + 300);
  };
}
function hideCard() { $('card').hidden = true; cardState = null; }

document.addEventListener('keydown', e => {
  if (e.target.matches('input, textarea, select')) { if (e.key === 'Escape') e.target.blur(); return; }
  if (!S.contigs) return;
  usePanel(S.active); // keys act on the active panel
  if (e.key === 'ArrowLeft') pan(-W * 0.15);
  else if (e.key === 'ArrowRight') pan(W * 0.15);
  else if (e.key === 'PageDown' || e.key === 'PageUp') vscrollBy((e.key === 'PageDown' ? 1 : -1) * (H - covTop()) * 0.9);
  else if (e.key === 'ArrowUp') { for (const T of S.tracks) T.scroll -= 60; clampScroll(); draw(); }
  else if (e.key === 'ArrowDown') { for (const T of S.tracks) T.scroll += 60; clampScroll(); draw(); }
  else if (e.key === '+' || e.key === '=') zoomAt(W / 2, 0.5);
  else if (e.key === '-' || e.key === '_') zoomAt(W / 2, 2);
  else if (e.key === 's' || e.key === 'ы') sortHere(); // 'ы' / 'р': same keys on a Russian layout
  else if (e.key === 'h' || e.key === 'р') { P.hideToolbar = !P.hideToolbar; applyChrome(); saveProfile(); }
  else if (e.key === 'Escape') { hideCard(); closePanels(); closeMenu(); $('toast').hidden = true; S.selected = null; S.hl = null; for (const p of S.panels) p.rubber = null; draw(); }
  else return;
  e.preventDefault();
});

function sortHere() {
  S.sortPos = Math.floor(viewCenter());
  if (P.sortBy === 'start') { P.sortBy = 'base'; syncPanel(); saveProfile(); }
  rebuild(); draw();
}

// ======================= Scrollbar of the sample stack =======================
const vbar = $('vbar'), vthumb = vbar.firstElementChild;
function vbarGeom() { const avail = H - covTop(), total = tracksHeight(), th = Math.max(24, avail * avail / total); return { avail, total, th }; }
function updateVBar() {
  const { avail, total, th } = vbarGeom(), on = S.tracks.length > 1 && total > avail + 1;
  vbar.hidden = !on; if (!on) return;
  vbar.style.top = covTop() + 'px'; vbar.style.height = avail + 'px';
  vthumb.style.height = th + 'px'; vthumb.style.top = (avail - th) * (S.vscroll || 0) / (total - avail) + 'px';
}
let vdrag = null;
vthumb.addEventListener('pointerdown', e => { e.stopPropagation(); vthumb.setPointerCapture(e.pointerId); vdrag = { y: e.clientY, v: S.vscroll || 0 }; });
vthumb.addEventListener('pointermove', e => {
  if (!vdrag) return;
  const { avail, total, th } = vbarGeom();
  S.vscroll = vdrag.v + (e.clientY - vdrag.y) * (total - avail) / Math.max(1, avail - th); clampVScroll(); draw();
});
vthumb.addEventListener('pointerup', () => { vdrag = null; });
// a click on the bar beside the thumb moves a page that way
vbar.addEventListener('pointerdown', e => { if (e.target !== vbar) return; const r = vthumb.getBoundingClientRect(); vscrollBy((e.clientY < r.top ? -1 : 1) * (H - covTop()) * 0.9); });

// ======================= Size =======================
function resize() {
  const r = stage.getBoundingClientRect();
  DPR = Math.min(window.devicePixelRatio || 1, 3);
  FULL_W = Math.max(1, r.width); H = Math.max(1, r.height);
  for (const c of [glCanvas, ov]) { c.width = Math.round(FULL_W * DPR); c.height = Math.round(H * DPR); }
  layoutPanels(false); usePanel(S.active); // each panel keeps its center and scale
  clampScroll(); draw(); syncLocus();
}
new ResizeObserver(resize).observe(stage);

// ======================= Opening data =======================
// Adds a sample as a new track. Its contigs are matched to the current coordinate system — the
// reference, or the first sample's header when no reference is loaded (name variants, then unique
// length); the very first sample without a reference gives the coordinates itself.
// The same file (or link) opened again is not added twice: every track loads its window, so a
// duplicate doubles the work for nothing.
function sourceKey(src) { return src.blob ? `${src.name}|${src.blob.size}|${src.blob.lastModified || 0}` : src.url || src.name; }
async function openBam(bamSrc, baiSrc, opts) {
  opts = opts || {};
  const key = sourceKey(bamSrc);
  if (S.tracks.some(x => x.key === key)) { toast(t('sample.already', { name: bamSrc.name })); return true; }
  S.loading = true; S.msg = null; draw();
  try {
    const bam = /\.cram$/i.test(bamSrc instanceof UrlSource ? urlPath(bamSrc.url) : bamSrc.name)
      ? await new CramFile(bamSrc, baiSrc).open() : await new BamFile(bamSrc, baiSrc).open();
    const back = new Int32Array(bam.refs.length).fill(-1);
    let map, note = '', noteSticky = false;
    const own = !S.fasta && !S.tracks.length; // first sample with no reference: its header is the genome
    if (!own) {
      const link = matchContigs(bam.refs, S.contigs);
      // a large contig with the same name but another length: the sample is from another genome build
      const big = link.conflicts.filter(c => Math.max(c[1], c[3]) > 1e6);
      if (big.length) {
        S.loading = false; draw();
        const eg = big.slice(0, 2).map(c => t('sample.lenPair', { bam: c[0], bl: fmtInt(c[1]), ref: c[2], rl: fmtInt(c[3]) })).join('; ');
        toast(t('sample.build', { name: bamSrc.name, eg }), true);
        return false;
      }
      if (!link.matched.length) {
        S.loading = false; draw();
        const list = (arr, n) => arr.slice(0, n).map(x => `${x.name} (${fmtInt(x.len)} ${t('unit.bp')})`).join(', ') + (arr.length > n ? ` … +${fmtInt(arr.length - n)}` : '');
        toast(t('sample.mismatch', { name: bamSrc.name, bam: list(bam.refs, 3), ref: list(S.contigs, 3) }), true);
        return false;
      }
      const viewIdx = new Map(S.contigs.map((c, i) => [c.name, i])), bamIdx = new Map(bam.refs.map((r, i) => [r.name, i]));
      map = new Int32Array(S.contigs.length).fill(-1);
      for (const [b, f] of link.matched) { map[viewIdx.get(f)] = bamIdx.get(b); back[bamIdx.get(b)] = viewIdx.get(f); }
      const byLen = link.matched.filter(m => m[2] === 'length');
      if (byLen.length) note += t('sample.byLength', { pairs: byLen.slice(0, 3).map(m => m[0] + ' = ' + m[1]).join(', ') + (byLen.length > 3 ? ' …' : '') }) + ' ';
      if (link.conflicts.length) noteSticky = true;
      if (link.conflicts.length) note += t('sample.lenDiff', { names: link.conflicts.slice(0, 3).map(c => t('sample.lenPair', { bam: c[0], bl: fmtInt(c[1]), ref: c[2], rl: fmtInt(c[3]) })).join('; ') }) + ' ';
      if (link.unmatched.length) note += t('sample.partial', { n: fmtInt(link.unmatched.length), names: link.unmatched.slice(0, 3).join(', ') + (link.unmatched.length > 3 ? ' …' : '') });
    } else {
      S.contigs = bam.refs.map(r => ({ name: r.name, len: r.len })); S._ck = null;
      map = Int32Array.from(bam.refs, (r, i) => i); back.forEach((v, i) => { back[i] = i; });
    }
    usePanel(S.active);
    const keepView = !own && S.viewSet && S.panels.every(pn => map[pn.view.refId] >= 0);
    if (!opts.demo) for (const x of S.tracks.filter(x => x.demo)) closeTrack(x); // own data replaces the demo sample
    const T = makeTrack({ id: ++S.trackSeq, key, bam, name: bamSrc.name, demo: !!opts.demo, map, back });
    // CRAM keeps reads as differences from the reference: its bases come from the loaded reference
    if (bam.isCram) bam.refProvider = async (id, a, b) => { const v = T.back[id]; return v >= 0 && S.fasta ? S.fasta.fetch(S.contigs[v].name, a, b) : null; };
    await attachWorker(T);
    usePanel(S.active);
    S.tracks.push(T);
    S.sortPos = S.tracks.length > 1 ? S.sortPos : null; S.selected = null; hideCard();
    updateFileLabel();
    $('examples').hidden = !S.tracks.some(x => x.demo);
    S.loading = false;
    if (note) toast(note.trim(), noteSticky);
    const loc = opts.locus ? parseLocus(opts.locus) : null;
    // every track of a panel loads the same window: each panel loads again with the new track
    forEachPanel(() => { S.loadToken++; S.data = null; });
    if (loc && !loc.error) { resetPanels(); goto(loc.id, loc.s, loc.e); }
    else if (keepView) { requestLoad(); draw(); }
    else { resetPanels(); const i = Math.max(0, map.findIndex(x => x >= 0)); goto(i, 0, Math.min(S.contigs[i].len, 10000)); }
    S.viewSet = true;
    return true;
  } catch (e) {
    usePanel(S.active);
    S.loading = false;
    if (S.tracks.length) toast(t('sample.failed', { name: bamSrc.name, msg: errText(e) }), true);
    else { S.msg = { text: errText(e), error: true }; }
    draw();
    return false;
  }
}
function closeTrack(T) {
  const i = S.tracks.indexOf(T); if (i < 0) return;
  S.tracks.splice(i, 1); forEachPanel(pn => { freeVao(T); pn.ts.delete(T.id); });
  if (T.worker) T.worker.postMessage({ op: 'close', key: T.id });
  if (S.selected && S.selected.T === T) { S.selected = null; hideCard(); }
  $('examples').hidden = !S.tracks.some(x => x.demo);
  updateFileLabel(); clampScroll(); rebuild(); draw();
}
function closeAllTracks() { for (const T of S.tracks) { forEachPanel(pn => { freeVao(T); pn.ts.delete(T.id); }); if (T.worker) T.worker.postMessage({ op: 'close', key: T.id }); } S.tracks = []; S.selected = null; S.hl = null; hideCard(); }
function updateFileLabel() {
  const parts = [];
  if (S.fasta) parts.push(S.fasta.name);
  if (S.tracks.length === 1) { const T = S.tracks[0]; parts.push(T.name + (T.bam.samples.length ? ' · ' + T.bam.samples.join(', ') : '')); }
  else if (S.tracks.length > 1) parts.push(t('tracks.n', { n: S.tracks.length }));
  $('fileName').textContent = parts.join('  ·  ');
  $('fileName').title = S.tracks.map(T => T.name).join('\n');
}

function b64ToBlob(b64) { const bin = atob(b64); const u = new Uint8Array(bin.length); for (let i = 0; i < bin.length; i++) u[i] = bin.charCodeAt(i); return new Blob([u]); }
async function openDemo() {
  if (typeof DEMO === 'undefined') return;
  const ok = await openReference({ fa: new FileSource(b64ToBlob(DEMO.fa), 'demo.fa'), fai: new FileSource(b64ToBlob(DEMO.fai), 'demo.fa.fai') }, { quiet: true, noGoto: true });
  if (ok) openBam(new FileSource(b64ToBlob(DEMO.bam), 'demo.bam'), new FileSource(b64ToBlob(DEMO.bai), 'demo.bam.bai'), { demo: true, locus: 'demo:11,941-12,060' });
}

// Pairs each alignment file with its index by name: sample.bam ↔ sample.bam.bai or sample.bai,
// sample.cram ↔ sample.cram.crai or sample.crai (case ignored). A single file with a single index of
// the right kind and any name is a pair too.
const IDX_EXT = { '.bam': '.bai', '.cram': '.crai' };
const alnExt = name => { const n = name.toLowerCase(); return n.endsWith('.bam') ? '.bam' : n.endsWith('.cram') ? '.cram' : null; };
function pairBams(files) {
  const L = f => f.name.toLowerCase();
  const bams = files.filter(f => alnExt(f.name));
  const idx = new Map(files.filter(f => /\.(bai|crai)$/.test(L(f))).map(f => [L(f), f]));
  const pairs = [], missing = [];
  for (const bam of bams) {
    const n = L(bam), ext = alnExt(n), ie = IDX_EXT[ext];
    const bai = idx.get(n + ie) || idx.get(n.slice(0, -ext.length) + ie);
    if (bai) pairs.push([bam, bai]); else missing.push(bam);
  }
  if (bams.length === 1 && missing.length === 1) {
    const ie = IDX_EXT[alnExt(bams[0].name)], only = [...idx.values()].filter(f => L(f).endsWith(ie));
    if (only.length === 1) { pairs.push([bams[0], only[0]]); missing.length = 0; }
  }
  return { bams, pairs, missing };
}
// Files chosen or dropped: any number of samples (each with its index) and optionally a reference.
// A folder with several samples first lists them, so a whole run directory is not opened at once.
// ---- The index next to a chosen file (File System Access: Chrome, Edge) ----
// A web page reads only the files the user picks. Where the browser can grant a folder, choosing just
// the .bam or .cram is enough: CGV asks once to read its folder and takes the index from there.
// Granted folders are read-only and kept for this session only (never stored), so later files from
// the same folder open without asking again.
const FSA = typeof showOpenFilePicker === 'function' && typeof showDirectoryPicker === 'function';
const grantedDirs = [];
// Is this the folder the file came from? Same name, size and modification time.
async function dirHolds(dir, f) {
  try { const g = await (await dir.getFileHandle(f.name)).getFile(); return g.size === f.size && g.lastModified === f.lastModified; } catch (e) { return false; }
}
// The index of an alignment file in a granted folder, by the same names as pairBams.
async function indexInDirs(f) {
  const n = f.name.toLowerCase(), ext = alnExt(n); if (!ext) return null;
  const want = [n + IDX_EXT[ext], n.slice(0, -ext.length) + IDX_EXT[ext]];
  for (const dir of grantedDirs) {
    if (!(await dirHolds(dir, f))) continue;
    for await (const [name, h] of dir.entries()) if (h.kind === 'file' && want.includes(name.toLowerCase())) return h.getFile();
  }
  return null;
}
// Asks for the folder of the files left without an index; the picker opens in that folder.
function askFolder(missing, files, fromFolder) {
  const el = $('toast'), withHandle = missing.find(f => f._h);
  el.innerHTML = `<span>${esc(t('files.askFolder', { names: missing.map(f => f.name).join(', ') }))}</span> <button class="btn primary" id="allowDir">${esc(t('files.allowFolder'))}</button><button class="toast-x" aria-label="${esc(t('btn.close'))}">×</button>`;
  el.hidden = false; clearTimeout(toastTimer);
  $('allowDir').onclick = async ev => {
    ev.stopPropagation(); el.hidden = true;
    let dir;
    try { dir = await showDirectoryPicker({ id: 'cgv-data', mode: 'read', startIn: withHandle ? withHandle._h : 'documents' }); }
    catch (e) { return; } // closed without choosing
    if (!(await dirHolds(dir, missing[0]))) { toast(t('files.dirWrong', { name: missing[0].name }), true); return; }
    grantedDirs.push(dir);
    openLocalFiles(files, fromFolder);
  };
}
async function openLocalFiles(fileList, fromFolder) {
  const files = [...fileList], L = f => f.name.toLowerCase();
  $('toast').hidden = true; // a message about an earlier choice no longer applies
  if (FSA) { // indexes from folders this session can already read
    for (const f of pairBams(files).missing) { const idx = await indexInDirs(f); if (idx) files.push(idx); }
    const left = pairBams(files).missing;
    if (left.length && left.some(f => f._h)) { askFolder(left, files, fromFolder); return; }
  }
  const { bams, pairs, missing } = pairBams(files);
  const ref = pickRefFiles(files);
  if (!bams.length) {
    if (ref) { openReference(ref); return; }
    toast(t(files.some(f => /\.(bai|crai)$/.test(L(f))) ? 'files.needBam' : 'files.needPair'));
    return;
  }
  if (fromFolder && bams.length > 1) { showFolderList(files, pairs, missing); return; }
  if (!pairs.length) {
    const csi = files.some(f => L(f).endsWith('.csi'));
    toast(csi ? t('files.csi') : t('files.noIndex', { names: missing.map(f => f.name).join(', ') }), true);
    return;
  }
  closePanels();
  const notes = [];
  if (missing.length) notes.push(t('files.noIndexSome', { names: missing.map(f => f.name).join(', ') }));
  // a reference chosen together with the samples is loaded first; the samples are then matched to it
  (ref ? openReference(ref, { noGoto: true }) : Promise.resolve(true)).then(async ok => {
    if (!ok) return;
    for (const [bam, bai] of pairs) await openBam(new FileSource(bam), new FileSource(bai));
    if (notes.length) toast(notes.join(' '), true); // after the reference's own message
  });
}
// A folder with several samples: a checklist in the Open panel.
function showFolderList(files, pairs, missing) {
  const box = $('folderList');
  const rows = pairs.map(([bam], i) => `<label class="chk"><input type="checkbox" data-i="${i}"${i < 4 ? ' checked' : ''}> ${esc(bam.name)}</label>`).join('') +
    missing.map(f => `<div class="hint">${esc(f.name)} — ${esc(t('folder.noIndex'))}</div>`).join('');
  box.innerHTML = `<p class="hint">${esc(t('folder.found', { n: pairs.length + missing.length }))}</p>${rows}` +
    (pairs.length ? `<div class="row-actions"><button class="btn primary" id="folderOpen">${esc(t('folder.open'))}</button></div>` : '');
  box.hidden = false;
  $('openPanel').hidden = false; $('openBtn').setAttribute('aria-expanded', 'true');
  if (!pairs.length) return;
  $('folderOpen').onclick = () => {
    const chosen = [...box.querySelectorAll('input[data-i]:checked')].map(x => pairs[+x.dataset.i]);
    if (!chosen.length) { toast(t('folder.none')); return; }
    box.hidden = true; box.innerHTML = '';
    const ref = pickRefFiles(files);
    openLocalFiles([...chosen.flat(), ...(ref ? files.filter(f => /\.(fa|fasta|fna)(\.gz)?$|\.fai$|\.gzi$/i.test(f.name)) : [])]);
  };
}
// Files of a dropped folder (and its subfolders, a few levels deep).
async function droppedFiles(dt) {
  // handles must be asked for during the drop event itself; they let CGV find an index next to a file
  const hp = [...dt.items].map(it => it.kind === 'file' && it.getAsFileSystemHandle ? it.getAsFileSystemHandle().catch(() => null) : null);
  const entries = [...dt.items].map(it => it.webkitGetAsEntry && it.webkitGetAsEntry()).filter(Boolean);
  if (!entries.some(e => e.isDirectory)) {
    const files = [...dt.files], hs = await Promise.all(hp);
    hs.forEach(h => { if (!h) return; if (h.kind === 'directory') grantedDirs.push(h); else { const f = files.find(x => x.name === h.name); if (f) f._h = h; } });
    return { files, folder: false };
  }
  const out = [];
  const walk = async (e, depth) => {
    if (e.isFile) { out.push(await new Promise((res, rej) => e.file(res, rej))); return; }
    if (depth > 3) return;
    const reader = e.createReader();
    for (;;) { const batch = await new Promise((res, rej) => reader.readEntries(res, rej)); if (!batch.length) break; for (const c of batch) await walk(c, depth + 1); }
  };
  for (const e of entries) await walk(e, 0);
  return { files: out, folder: true };
}

// Links: one per line (or separated by spaces). The index is looked for next to each BAM:
// sample.bam.bai, then sample.bai. The index field applies when there is one link.
function openUrls() {
  const links = $('urlBam').value.split(/\s+/).map(x => x.trim()).filter(Boolean);
  $('toast').hidden = true;
  if (!links.length) { toast(t('url.need')); return; }
  const auth = $('urlAuth').value.trim();
  const headers = auth ? { Authorization: /^bearer\s/i.test(auth) ? auth : 'Bearer ' + auth } : {};
  const bams = links;
  const given = bams.length === 1 ? $('urlBai').value.trim() : '';
  closePanels();
  (async () => {
    for (const u of bams) {
      const cram = /\.cram$/i.test(urlPath(u)), ext = cram ? '.cram' : '.bam', ie = IDX_EXT[ext];
      const cands = given ? [given] : [addSuffix(u, ie), swapExt(u, ext, ie)].filter((x, i, a) => x && a.indexOf(x) === i);
      await openBam(new UrlSource(u, headers), cands.map(x => new UrlSource(x, headers)));
    }
  })();
}
function urlPath(u) { try { return new URL(u).pathname; } catch (e) { return u.split(/[?#]/)[0]; } }
function swapExt(url, from, to) {
  try { const x = new URL(url); if (!x.pathname.toLowerCase().endsWith(from)) return null; x.pathname = x.pathname.slice(0, -from.length) + to; return x.toString(); }
  catch (e) { return url.toLowerCase().endsWith(from) ? url.slice(0, -from.length) + to : null; }
}
function addSuffix(url, sfx) { try { const x = new URL(url); x.pathname += sfx; return x.toString(); } catch (e) { return url + sfx; } }


// ======================= Reference =======================
// Public assemblies hosted by the IGV team on S3 (CORS enabled); same list igv.js uses.
const REF_PRESETS = {
  hg38: { label: 'GRCh38 / hg38', fa: 'https://igv-genepattern-org.s3.amazonaws.com/genomes/seq/hg38/hg38.fa' },
  hg19: { label: 'GRCh37 / hg19', fa: 'https://igv-genepattern-org.s3.amazonaws.com/genomes/seq/hg19/hg19.fasta' },
  chm13: { label: 'T2T-CHM13 v2.0', fa: 'https://s3.amazonaws.com/igv.org.genomes/chm13v2.0/chm13v2.0.fa' },
};
function pickRefFiles(files) {
  const L = f => f.name.toLowerCase();
  const fa = [...files].find(f => /\.(fa|fasta|fna)(\.gz)?$/.test(L(f)));
  if (!fa) return null;
  const fai = [...files].find(f => L(f).endsWith('.fai')) || null;
  const gzi = [...files].find(f => L(f).endsWith('.gzi')) || null;
  return { fa: new FileSource(fa), fai: fai ? new FileSource(fai) : null, gzi: gzi ? new FileSource(gzi) : null, size: fa.size };
}
// ref = { fa, fai?, gzi?, size? } — sources. Without a .fai (or a .fa.gz without .gzi)
// the FASTA is indexed here in one pass, with progress.
let refBusy = false;
// Loading a reference switches the genome: samples opened against the previous one are closed.
async function openReference(ref, opts) {
  opts = opts || {};
  if (refBusy) return false;
  refBusy = true;
  const name = ref.fa.name, gz = /\.gz$/i.test(name);
  try {
    let fasta = null;
    if (ref.fai && (!gz || ref.gzi)) {
      try { fasta = await new IndexedFasta(ref.fa, ref.fai, ref.gzi).open(); }
      catch (e) { if (!(ref.fa instanceof UrlSource) || !['notfound', 'http', 'auth', 'network'].includes(e.code) || e.vars.name === name) throw e; fasta = null; }
    }
    if (!fasta) {
      const fmt = n => n >= 1 << 30 ? (n / (1 << 30)).toFixed(1) + ' GB' : (n / (1 << 20)).toFixed(0) + ' MB';
      toast(t('ref.indexing', { name, done: '0%' }), true);
      fasta = await new IndexedFasta(ref.fa).openUnindexed(ref.size, (done, total) =>
        toast(t('ref.indexing', { name, done: total ? Math.min(100, Math.round(done / total * 100)) + '%' : fmt(done) }), true));
    }
    if (opts.ifIdle && (S.tracks.length || S.fasta)) return 'busy'; // opened at start, but the user was quicker
    const closed = S.tracks.length > 0;
    resetPanels();
    S.loadToken++;
    S.fasta = fasta; S.contigs = [...fasta.idx].map(([n, e]) => ({ name: n, len: e.len })); S._ck = null;
    closeAllTracks();
    S.data = null; S.sortPos = null; S.msg = null; S.loading = false;
    hideCard(); updateRefInfo(); updateFileLabel(); $('examples').hidden = true; closePanels();
    if (!opts.quiet) {
      let note = fasta.built ? t('ref.indexed', { name, n: fmtInt(fasta.idx.size) }) : t('ref.loaded', { name });
      if (closed) note += ' ' + t('ref.closedSamples');
      toast(note, note.length > 90);
    }
    // a new genome means a new coordinate system: the old view no longer applies
    S.view.refId = 0; S.view.start = 0; S.viewSet = false;
    if (!opts.noGoto) { goto(0, 0, Math.min(S.contigs[0].len, 1000)); S.viewSet = true; }
    return true;
  } catch (e) { toast(t('ref.failed', { msg: errText(e) })); return false; }
  finally { refBusy = false; }
}
function updateRefInfo() {
  const el = $('refInfo'); if (!el) return;
  el.textContent = S.fasta ? t('ref.current', { name: S.fasta.name }) : t('ref.none');
}
$('pickRef').onchange = e => { const r = pickRefFiles(e.target.files); e.target.value = ''; if (!r) { toast(t('ref.failed', { msg: t('err.faEmpty') })); return; } openReference(r); };
$('refPreset').onchange = e => { const p = REF_PRESETS[e.target.value]; if (p) $('refUrl').value = p.fa; };
$('loadRef').onclick = () => {
  const u = $('refUrl').value.trim(); if (!u) return;
  const gz = /\.gz(\?|$)/i.test(u);
  openReference({ fa: new UrlSource(u), fai: new UrlSource(addSuffix(u, '.fai')), gzi: gz ? new UrlSource(addSuffix(u, '.gzi')) : null });
};

// ======================= Language =======================
function buildLanguageSelect() {
  const sel = $('language');
  sel.innerHTML = LANGS.map(([code, name]) => `<option value="${code}">${esc(name)}</option>`).join('');
}
function applyLanguage() {
  LANG = STRINGS[P.language] ? P.language : 'en';
  applyI18n(); updateRefInfo(); syncPanel(); syncLocus();
  if (S.data) rebuild(); draw();
  if (S.selected) { hideCard(); S.selected = null; }
}

// ======================= Context menu (right click on the reads) =======================
// Pair view, coloring and sorting live here, as in IGV; sorting uses the clicked position.
// x: inside the current panel; gx, y: on the whole view. hit: the read under the pointer, if any.
let menuState = null;
function openMenu(x, gx, y, hit) {
  const m = $('ctx'), pos = Math.max(0, Math.floor(S.view.start + x * S.view.bpPerPx));
  const item = (act, val, on, label, kind) => `<button class="mi" role="${kind}" aria-checked="${on}" data-act="${act}" data-val="${val}"><span class="mk" aria-hidden="true">${on ? (kind === 'menuitemcheckbox' ? '✓' : '●') : ''}</span>${esc(label)}</button>`;
  const r0 = hit && hit.rec, mateOk = r0 && (r0.flag & 1) && !(r0.flag & 8) && r0.mateRef >= 0;
  menuState = { pn: CUR, hit: mateOk ? hit : null };
  m.innerHTML =
    (mateOk ? `<button class="mi" role="menuitem" data-act="mateSplit" data-val=""><span class="mk" aria-hidden="true"></span>${esc(t('menu.mateSplit'))}</button><div class="msep"></div>` : '') +
    item('pairs', '', !!P.viewPairs, t('menu.pairs'), 'menuitemcheckbox') +
    item('clips', '', !!P.showSoftClips, t('menu.clips'), 'menuitemcheckbox') +
    `<div class="ms">${esc(t('menu.colorBy'))}</div>` +
    ['none', 'strand', 'pair', 'tag'].map(v => item('color', v, P.colorBy === v, t('color.' + v), 'menuitemradio')).join('') +
    `<div class="mt"><label for="ctxTag">${esc(t('view.colorTag'))}</label><input id="ctxTag" type="text" maxlength="2" spellcheck="false" autocomplete="off" value="${esc(P.colorTag)}"></div>` +
    `<div class="ms">${esc(t('menu.sortBy', { pos: fmtInt(pos + 1) }))}</div>` +
    ['start', 'base', 'mapq', 'strand', 'insert'].map(v => item('sort', v, P.sortBy === v && (v === 'start' || S.sortPos === pos), t('sort.' + v), 'menuitemradio')).join('');
  m.dataset.pos = pos; m.hidden = false;
  const r = stage.getBoundingClientRect(), w = m.offsetWidth, h = m.offsetHeight;
  m.style.left = Math.max(8, Math.min(r.width - w - 8, gx + 2)) + 'px';
  m.style.top = Math.max(8, Math.min(r.height - h - 8, y + 2)) + 'px';
  const tag = $('ctxTag');
  tag.oninput = () => { const v = tag.value.trim(); if (v.length === 2) { P.colorTag = v; saveProfile(); if (P.colorBy === 'tag') { rebuild(); draw(); } } };
  tag.onkeydown = e => { if (e.key === 'Enter' || e.key === 'Escape') { e.preventDefault(); closeMenu(); } };
  (m.querySelector('[aria-checked="true"]') || m.querySelector('.mi')).focus({ preventScroll: true });
}
function closeMenu() { $('ctx').hidden = true; }
$('ctx').addEventListener('click', e => {
  const b = e.target.closest('[data-act]'); if (!b) return;
  const act = b.dataset.act, val = b.dataset.val, pos = +$('ctx').dataset.pos; let reloadNeeded = false;
  if (!menuState || !S.panels.includes(menuState.pn)) { closeMenu(); return; }
  usePanel(menuState.pn);
  if (act === 'mateSplit') { const h = menuState.hit; closeMenu(); hideCard(); if (h) openMatePanel(h.rec, h.T); draw(); return; }
  if (act === 'pairs') P.viewPairs = !P.viewPairs;
  else if (act === 'clips') { P.showSoftClips = !P.showSoftClips; reloadNeeded = true; }
  else if (act === 'color') P.colorBy = val;
  else if (act === 'sort') { P.sortBy = val; S.sortPos = val === 'start' ? null : pos; }
  saveProfile(); syncPanel(); closeMenu(); hideCard(); S.selected = null;
  if (reloadNeeded) reloadData(); else { rebuild(); draw(); }
});
$('ctx').addEventListener('keydown', e => {
  const items = [...$('ctx').querySelectorAll('.mi')], i = items.indexOf(document.activeElement);
  if (e.key === 'ArrowDown' || e.key === 'ArrowUp') { e.preventDefault(); const n = items.length; items[((i < 0 ? -1 : i) + (e.key === 'ArrowDown' ? 1 : -1) + n) % n].focus(); }
  else if (e.key === 'Escape') { e.preventDefault(); closeMenu(); }
});
stage.addEventListener('contextmenu', e => {
  if (!onCanvas(e)) return;
  e.preventDefault();
  if (!S.tracks.length || !S.contigs) return;
  const L = local(e); if (L.y < panelHdr()) return;
  setActive(L.pn);
  hideCard(); openMenu(L.x, e.offsetX, L.y, hitTest(L.x, L.y));
});
document.addEventListener('pointerdown', e => { if (!$('ctx').hidden && !$('ctx').contains(e.target)) closeMenu(); }, true);

// ======================= Panels and settings =======================
function closePanels() { $('panel').hidden = true; $('openPanel').hidden = true; $('viewBtn').setAttribute('aria-expanded', 'false'); $('openBtn').setAttribute('aria-expanded', 'false'); }
function togglePanel(id, btn) { const el = $(id), was = el.hidden; closePanels(); el.hidden = !was; btn.setAttribute('aria-expanded', String(was)); if (id === 'panel' && was) syncPanel(); }
$('viewBtn').onclick = () => togglePanel('panel', $('viewBtn'));
$('openBtn').onclick = () => togglePanel('openPanel', $('openBtn'));
document.querySelectorAll('[data-close]').forEach(b => b.onclick = closePanels);
$('pickFiles').onchange = e => { if (e.target.files.length) openLocalFiles(e.target.files); e.target.value = ''; };
// Chrome, Edge: the browser's own file picker with handles, so the folder can be asked for when needed
document.querySelectorAll('label[for="pickFiles"]').forEach(lb => lb.addEventListener('click', async e => {
  if (!FSA) return;
  e.preventDefault();
  let hs;
  try { hs = await showOpenFilePicker({ multiple: true, id: 'cgv-data' }); }
  catch (err) { if (err.name !== 'AbortError') $('pickFiles').click(); return; } // not allowed here (an embedded page): the plain picker
  const files = await Promise.all(hs.map(async h => { const f = await h.getFile(); f._h = h; return f; }));
  openLocalFiles(files);
}));
$('pickFolder').onchange = e => { if (e.target.files.length) openLocalFiles(e.target.files, true); e.target.value = ''; };
$('openUrl').onclick = openUrls;
$('openDemo').onclick = () => { openDemo(); closePanels(); };
stage.addEventListener('dragover', e => { e.preventDefault(); stage.classList.add('drop'); });
stage.addEventListener('dragleave', () => stage.classList.remove('drop'));
stage.addEventListener('drop', async e => {
  e.preventDefault(); stage.classList.remove('drop');
  try { const { files, folder } = await droppedFiles(e.dataTransfer); if (files.length) openLocalFiles(files, folder); }
  catch (err) { if (e.dataTransfer.files.length) openLocalFiles(e.dataTransfer.files); }
});

// Go to the typed locus: Enter in the field or the Go button.
function goLocus() {
  if (!S.contigs) return;
  const el = $('locus'), l = parseLoci(el.value);
  if (!l) { toast(t('locus.bad')); return; }
  if (l.error) { toast(l.error); return; }
  if (l.length > maxPanels()) toast(t('split.full', { n: maxPanels() }));
  el.blur(); setPanels(l);
}
$('locus').addEventListener('keydown', e => { if (e.key === 'Enter') goLocus(); });
$('goBtn').onclick = goLocus;
$('examples').onchange = e => { const v = e.target.value; e.target.value = ''; if (!v || !S.contigs) return; const l = parseLocus(v); if (l && !l.error) setPanels([l]); };

// profile controls: id → [profile path, type, what to recompute]
const CONTROLS = [
  ['readHeight', 'readHeight', 'num', 'draw'],
  ['colorBy', 'colorBy', 'str', 'rebuild'], ['colorTag', 'colorTag', 'str', 'rebuild'],
  ['groupBy', 'groupBy', 'str', 'rebuild'], ['groupTag', 'groupTag', 'str', 'rebuild'],
  ['sortBy', 'sortBy', 'str', 'rebuild'],
  ['mqLow', 'mapqAlpha.low', 'num', 'rebuild'], ['mqHigh', 'mapqAlpha.high', 'num', 'rebuild'], ['mqMinAlpha', 'mapqAlpha.minAlpha', 'num', 'rebuild'],
  ['bqShade', 'baseQualityShading', 'bool', 'rebuild'], ['softClips', 'showSoftClips', 'bool', 'reload'],
  ['minMapq', 'filters.minMapq', 'num', 'reload'], ['hideDup', 'filters.hideDuplicates', 'bool', 'reload'],
  ['hideSec', 'filters.hideSecondary', 'bool', 'reload'], ['hideSupp', 'filters.hideSupplementary', 'bool', 'reload'], ['hideQc', 'filters.hideQcFail', 'bool', 'reload'],
  ['dsOn', 'downsample.enabled', 'bool', 'reload'], ['dsMax', 'downsample.maxReads', 'num', 'reload'], ['dsWin', 'downsample.window', 'num', 'reload'],
  ['covHeight', 'coverage.height', 'num', 'draw'], ['trackMinH', 'trackMinHeight', 'num', 'draw'], ['alleleThr', 'coverage.alleleThreshold', 'num', 'rebuild'],
  ['hideToolbar', 'hideToolbar', 'bool', 'chrome'], ['showGenes', 'showGenes', 'bool', 'draw'], ['startRef', 'startReference', 'str', 'none'], ['language', 'language', 'str', 'lang'],
];
const getPath = (o, p) => p.split('.').reduce((a, k) => a[k], o);
const setPath = (o, p, v) => { const ks = p.split('.'), last = ks.pop(); ks.reduce((a, k) => a[k], o)[last] = v; };
function syncPanel() {
  for (const [id, path, type] of CONTROLS) { const el = $(id); if (!el) continue; const v = getPath(P, path); if (type === 'bool') el.checked = v; else el.value = v; }
  $('readHeightOut').textContent = (+P.readHeight).toFixed(1) + ' px';
  $('mqMinAlphaOut').textContent = Math.round(P.mapqAlpha.minAlpha * 100) + '%';
  $('alleleThrOut').textContent = Math.round(P.coverage.alleleThreshold * 100) + '%';
  $('covHeightOut').textContent = P.coverage.height + ' px';
  $('trackMinHOut').textContent = P.trackMinHeight + ' px';
  $('groupTagRow').hidden = P.groupBy !== 'tag';
  $('profileJson').value = JSON.stringify(P, null, 2);
  updateMqLegend();
}
function updateMqLegend() {
  const el = $('mqLegend'); if (!el) return;
  el.style.background = `linear-gradient(90deg, rgb(${COL.read} / ${P.mapqAlpha.minAlpha}), rgb(${COL.read} / 1))`;
  $('mqLegendLo').textContent = 'MQ ≤ ' + P.mapqAlpha.low; $('mqLegendHi').textContent = 'MQ ≥ ' + P.mapqAlpha.high;
}
let rebuildTimer = 0;
for (const [id, path, type, effect] of CONTROLS) {
  const el = $(id); if (!el) continue;
  el.addEventListener(type === 'bool' || el.tagName === 'SELECT' ? 'change' : 'input', () => {
    let v = type === 'bool' ? el.checked : type === 'num' ? parseFloat(el.value) : el.value.trim();
    if (type === 'num' && !isFinite(v)) return;
    setPath(P, path, v); saveProfile(); syncPanel();
    if (effect === 'draw') { clampScroll(); draw(); }
    else if (effect === 'chrome') applyChrome();
    else if (effect === 'lang') applyLanguage();
    else if (effect === 'reload') { clearTimeout(rebuildTimer); rebuildTimer = setTimeout(reloadData, 250); }
    else { clearTimeout(rebuildTimer); rebuildTimer = setTimeout(() => { usePanel(S.active); if (id === 'sortBy' && S.sortPos == null) S.sortPos = Math.floor(viewCenter()); rebuild(); draw(); }, 30); }
  });
}
$('copyProfile').onclick = () => navigator.clipboard.writeText(JSON.stringify(P, null, 2)).then(() => toast(t('profile.copied')), () => { $('profileJson').select(); toast(t('copyFail')); });
$('applyProfile').onclick = () => {
  try { P = mergeProfile(DEFAULT_PROFILE, JSON.parse($('profileJson').value)); saveProfile(); applyLanguage(); syncPanel(); applyChrome(); rebuild(); draw(); toast(t('profile.applied')); }
  catch (e) { toast(t('profile.bad', { msg: e.message })); }
};
$('resetProfile').onclick = () => { P = clone(DEFAULT_PROFILE); saveProfile(); applyLanguage(); syncPanel(); applyChrome(); rebuild(); draw(); toast(t('profile.reset')); };
$('showToolbar').onclick = () => { P.hideToolbar = false; saveProfile(); applyChrome(); };
function applyChrome() { document.body.classList.toggle('bare', P.hideToolbar); $('showToolbar').hidden = !P.hideToolbar; requestAnimationFrame(resize); }

// Messages at the bottom: closed by a click (× or anywhere on it) or Escape. Short ones go by
// themselves; important ones (sticky) stay until closed.
let toastTimer = 0;
$('toast').onclick = () => { $('toast').hidden = true; };
function toast(msg, sticky) {
  const el = $('toast');
  el.innerHTML = `<span>${esc(msg)}</span><button class="toast-x" aria-label="${esc(t('btn.close'))}">×</button>`;
  el.hidden = false; clearTimeout(toastTimer);
  if (!sticky) toastTimer = setTimeout(() => el.hidden = true, msg.length > 90 ? 9000 : 3500);
}

// ======================= Theme =======================
function onTheme() { readColors(); if (S.data) rebuild(); draw(); updateMqLegend(); }
matchMedia('(prefers-color-scheme: dark)').addEventListener('change', onTheme);
new MutationObserver(onTheme).observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] });

// ======================= Start =======================
const SANS = getComputedStyle(document.documentElement).getPropertyValue('--sans') || 'system-ui, sans-serif';
const MONO = getComputedStyle(document.documentElement).getPropertyValue('--mono') || 'ui-monospace, monospace';
LANG = STRINGS[P.language] ? P.language : 'en';
buildLanguageSelect(); applyI18n(); readColors(); syncPanel(); applyChrome(); updateRefInfo();
if (!GL_OK) { $('msg').hidden = false; $('msg').className = 'msg err'; $('msg').textContent = t('nowebgl'); }
else { warmUpGL(); if (GL_SOFTWARE) setTimeout(() => toast(t('gl.software')), 1500); }

// ======================= Genes =======================
// Transcripts of the current panel drawn in two rows: line = intron (arrows give the strand),
// thin boxes = UTR, thick boxes = coding exons, the symbol under the transcript where it fits.
function drawGenes(pn) {
  pn.geneHits = [];
  if (!geneH()) return;
  const y0 = geneTop(), vs = S.view.start, ve = viewEnd(), bpp = S.view.bpPerPx;
  const txs = genesIn(S.contigs[S.view.refId].name, Math.floor(vs), Math.ceil(ve));
  ctx.save(); ctx.beginPath(); ctx.rect(0, y0, W, GENE_H); ctx.clip();
  ctx.font = '600 10.5px ' + SANS; ctx.textAlign = 'center'; ctx.textBaseline = 'alphabetic';
  const ends = [-Infinity, -Infinity], px = p => (p - vs) / bpp;
  for (const tx of txs) {
    const x0 = px(tx.start), x1 = px(tx.end), a = Math.max(x0, 0), b = Math.min(x1, W);
    const label = tx.sym + (tx.clinical ? ' (Plus Clinical)' : ''), lw = ctx.measureText(label).width;
    const cx = Math.max(lw / 2 + 2, Math.min(W - lw / 2 - 2, (a + b) / 2));
    const e0 = Math.min(a, cx - lw / 2), e1 = Math.max(b, cx + lw / 2);
    let row = ends.findIndex(e => e + 6 <= e0), withLabel = true;
    if (row < 0) { row = ends[0] <= ends[1] ? 0 : 1; withLabel = false; }
    ends[row] = Math.max(ends[row], withLabel ? e1 : b);
    const top = y0 + row * GENE_ROW, mid = top + 6.5, col = tx.clinical ? COL.css.geneClin : COL.css.gene;
    ctx.fillStyle = col; ctx.strokeStyle = col; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(a, mid + 0.5); ctx.lineTo(b, mid + 0.5); ctx.stroke();
    if (b - a > 30) { // strand arrows on the line
      ctx.beginPath();
      for (let x = a + 10; x < b - 6; x += 22) { const d = tx.strand > 0 ? 2.5 : -2.5; ctx.moveTo(x - d, mid - 2); ctx.lineTo(x + d, mid + 0.5); ctx.lineTo(x - d, mid + 3); }
      ctx.stroke();
    }
    for (let i = 0; i < tx.exons.length; i += 2) {
      const s = tx.exons[i], e = tx.exons[i + 1]; if (e < vs || s > ve) continue;
      const ex0 = px(s), ex1 = px(e);
      ctx.fillRect(ex0, mid - 3, Math.max(1, ex1 - ex0), 7); // UTR height
      const c0 = Math.max(s, tx.cdsStart), c1 = Math.min(e, tx.cdsEnd);
      if (c1 > c0) ctx.fillRect(px(c0), mid - 5.5, Math.max(1, px(c1) - px(c0)), 12);
      // exon number (in the transcript's direction) inside the exon when there is room
      if (ex1 - ex0 >= 26) {
        const n = tx.strand > 0 ? i / 2 + 1 : tx.exons.length / 2 - i / 2, xa = Math.max(ex0, 0), xb = Math.min(ex1, W);
        if (xb - xa >= 26) { ctx.save(); ctx.font = '600 9px ' + SANS; ctx.fillStyle = COL.css.panel; ctx.textBaseline = 'middle'; ctx.fillText(t('gene.exShort', { n }), (xa + xb) / 2, mid + 0.5); ctx.restore(); }
      }
    }
    if (withLabel) ctx.fillText(label, cx, top + GENE_ROW - 1);
    pn.geneHits.push({ tx, row, x0: Math.min(a, withLabel ? e0 : a), x1: Math.max(b, withLabel ? e1 : b) });
  }
  ctx.strokeStyle = COL.css.line; ctx.beginPath(); ctx.moveTo(0, y0 + GENE_H - 1.5); ctx.lineTo(W, y0 + GENE_H - 1.5); ctx.stroke();
  ctx.restore();
}
function geneHit(x, y) {
  const row = Math.floor((y - geneTop()) / GENE_ROW), hits = CUR.geneHits || [];
  return (hits.find(h => h.row === row && x >= h.x0 - 2 && x <= h.x1 + 2) || hits.find(h => x >= h.x0 - 2 && x <= h.x1 + 2) || {}).tx || null;
}
// Details of a transcript, and where the clicked position falls in it (exon / intron number).
function showGeneCard(tx, pos, gx, y) {
  const c = $('card'), name = S.contigs[S.view.refId].name, part = txPart(tx, pos);
  const L = [['kv', t('gene.symbol'), tx.sym + (tx.hgnc ? ' (' + tx.hgnc + ')' : '')], ['txt', tx.name], ['sep'],
    ['kv', t('gene.set'), tx.clinical ? 'MANE Plus Clinical' : 'MANE Select'], ['kv', 'Ensembl', tx.enst], ['kv', 'RefSeq', tx.nm],
    ['kv', t('gene.loc'), `${name}:${fmtInt(tx.start + 1)}-${fmtInt(tx.end)} (${tx.strand > 0 ? '+' : '−'})`],
    ['kv', t('gene.exons'), tx.exons.length / 2]];
  if (tx.cdsEnd > tx.cdsStart) L.push(['kv', t('gene.cds'), `${fmtInt(tx.cdsStart + 1)}-${fmtInt(tx.cdsEnd)}`]);
  const B = [['kv', t('card.location'), `${name}:${fmtInt(pos + 1)}`]];
  if (part) B.push(['kv', t('gene.where'), part.exon ? t('gene.exonN', { n: part.exon, of: part.of }) + (part.coding ? '' : ' (UTR)') : t('gene.intronN', { n: part.intron, of: part.of })]);
  c.classList.remove('wide');
  c.innerHTML = `<div class="cols" style="--n:1"><div class="col">${linesHtml(L)}<div class="loc">${linesHtml(B, 'b')}</div></div></div>` +
    `<div class="card-actions"><button class="btn" id="geneZoom">${esc(t('gene.zoom'))}</button><button class="btn" id="copyGene">${esc(t('card.copy'))}</button><button class="btn ghost" id="closeCard">${esc(t('btn.close'))}</button></div>`;
  c.hidden = false;
  const r = stage.getBoundingClientRect();
  c.style.left = Math.max(8, Math.min(r.width - c.offsetWidth - 8, gx + 14)) + 'px';
  c.style.top = Math.max(8, Math.min(r.height - c.offsetHeight - 8, y + 14)) + 'px';
  const pn = CUR;
  $('closeCard').onclick = () => { hideCard(); };
  $('copyGene').onclick = () => navigator.clipboard.writeText(linesText(L) + '\n' + linesText(B)).then(() => toast(t('copied')), () => toast(t('copyFail')));
  $('geneZoom').onclick = () => { if (!S.panels.includes(pn)) return; usePanel(pn); hideCard(); const pad = Math.max(50, Math.round((tx.end - tx.start) * 0.05)); goto(contigId(tx.chr), tx.start - pad, tx.end + pad); };
}
// Gene names offered while typing in the locus field (for the word being typed).
function updateGeneSuggest() {
  const dl = $('geneList'); if (!dl) return;
  const v = $('locus').value, word = v.slice(v.lastIndexOf(' ') + 1);
  const list = /[:\-]/.test(word) ? [] : geneSuggest(word, 12);
  dl.replaceChildren(...list.map(sym => { const o = document.createElement('option'); o.value = v.slice(0, v.length - word.length) + sym; return o; }));
}
$('locus').addEventListener('input', updateGeneSuggest);

// ======================= Chromosome map =======================
// The chromosome of every panel, side by side across the top when panels show different ones, with
// a red box where each panel is (numbered like the panel headers; the active panel's box filled).
// GRCh38 has its bands and centromere; other genomes get a plain bar. A click or drag moves the
// active panel there.
const IDEO = { segs: [] };
const STAIN = { gneg: 0.02, gpos25: 0.28, gpos50: 0.48, gpos75: 0.68, gpos100: 0.86, gvar: 0.36, stalk: 0.3 };
function drawIdeogram() {
  IDEO.segs = [];
  if (!ideoH()) return;
  const chroms = []; for (const pn of S.panels) if (!chroms.includes(pn.view.refId)) chroms.push(pn.view.refId);
  const k = chroms.length, gap = 14, segW = (FULL_W - gap * (k - 1)) / k;
  const bg = hexToRgb(COL.css.bg), fg = hexToRgb(COL.css.fg), mix = d => `rgb(${bg.map((c, i) => Math.round(c + (fg[i] - c) * d)).join(' ')})`;
  ctx.save();
  ctx.fillStyle = COL.css.bg; ctx.fillRect(0, 0, FULL_W, IDEO_H);
  chroms.forEach((rid, j) => {
    const x = j * (segW + gap), name = S.contigs[rid].name, len = S.contigs[rid].len;
    ctx.font = '600 11px ' + SANS; ctx.textAlign = 'left'; ctx.textBaseline = 'middle'; ctx.fillStyle = COL.css.fg;
    const lw = Math.min(ctx.measureText(name).width + 10, segW * 0.35);
    ctx.save(); ctx.beginPath(); ctx.rect(x, 0, lw - 4, IDEO_H); ctx.clip(); ctx.fillText(name, x + 4, 19); ctx.restore();
    const bx = x + lw, bw = Math.max(8, segW - lw - 6), by = 13, bh = 12;
    const bands = genesActive() ? cytoIn(name) : null;
    ctx.save();
    ctx.beginPath(); ctx.roundRect(bx, by, bw, bh, bh / 2); ctx.clip();
    ctx.fillStyle = mix(STAIN.gneg); ctx.fillRect(bx, by, bw, bh);
    if (bands) for (const b of bands) {
      const x0 = bx + b.s / len * bw, x1 = bx + b.e / len * bw;
      if (b.stain === 'acen') { ctx.fillStyle = '#c4544e'; ctx.fillRect(x0, by + 3, Math.max(0.5, x1 - x0), bh - 6); ctx.fillStyle = COL.css.bg; ctx.fillRect(x0, by, x1 - x0, 3); ctx.fillRect(x0, by + bh - 3, x1 - x0, 3); continue; }
      if (b.stain === 'gneg') continue;
      ctx.fillStyle = mix(STAIN[b.stain] ?? 0.4); ctx.fillRect(x0, by, Math.max(0.5, x1 - x0), bh);
    }
    ctx.restore();
    ctx.strokeStyle = mix(0.55); ctx.lineWidth = 1; ctx.beginPath(); ctx.roundRect(bx + 0.5, by + 0.5, bw - 1, bh - 1, bh / 2); ctx.stroke();
    // where the panels are
    S.panels.forEach((pn, i) => {
      if (pn.view.refId !== rid) return;
      const s = Math.max(0, pn.view.start), e = Math.min(len, pn.view.start + pn.w * pn.view.bpPerPx);
      let m0 = bx + s / len * bw, m1 = bx + e / len * bw;
      if (m1 - m0 < 3) { const c = (m0 + m1) / 2; m0 = c - 1.5; m1 = c + 1.5; }
      ctx.strokeStyle = '#e0312b'; ctx.lineWidth = 2; ctx.strokeRect(m0, by - 2, m1 - m0, bh + 4);
      if (pn === S.active) { ctx.fillStyle = 'rgb(224 49 43 / .28)'; ctx.fillRect(m0, by - 2, m1 - m0, bh + 4); }
      if (S.panels.length > 1) { ctx.font = '700 9px ' + SANS; ctx.textAlign = 'center'; ctx.fillStyle = '#e0312b'; ctx.fillText(String(i + 1), (m0 + m1) / 2, 6); }
    });
    IDEO.segs.push({ rid, x, w: segW, bx, bw, len });
  });
  ctx.strokeStyle = COL.css.line; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(0, IDEO_H - 0.5); ctx.lineTo(FULL_W, IDEO_H - 0.5); ctx.stroke();
  ctx.restore();
}
// A click or drag on the map: the active panel moves there (to that chromosome if it is another one).
function ideoGo(x) {
  const seg = IDEO.segs.find(g => x >= g.x && x < g.x + g.w + 14) || IDEO.segs[IDEO.segs.length - 1]; if (!seg) return;
  usePanel(S.active);
  const pos = Math.max(0, Math.min(seg.len, (x - seg.bx) / seg.bw * seg.len)), span = W * S.view.bpPerPx;
  goto(seg.rid, pos - span / 2, pos + span / 2);
}

// ======================= Startup =======================
// The bundled genes load in the background; the chosen reference (GRCh38 by default) opens.
async function startup() {
  const cyto = $('cgv-cyto'); if (cyto) parseCyto(cyto.textContent);
  const mane = $('cgv-mane');
  if (mane && mane.textContent.trim()) loadGenes(mane.textContent.trim()).then(() => { S._g38 = null; draw(); }).catch(e => console.error(e));
  const what = P.startReference;
  if (what === 'demo') { openDemo(); return; }
  const pr = REF_PRESETS[what]; if (!pr) return;
  $('refUrl').value = pr.fa;
  const ok = await openReference({ fa: new UrlSource(pr.fa), fai: new UrlSource(addSuffix(pr.fa, '.fai')) }, { quiet: true, noGoto: true, ifIdle: true });
  if (ok === 'busy') return; // the user opened something first
  if (!ok) { toast(t('start.refFailed', { name: pr.label }), true); openDemo(); return; }
  toast(t('start.ref', { name: pr.label }));
  const i = Math.max(0, contigId('chr1'));
  // waits for the genes (a moment) so the first view can open on a gene
  for (let k = 0; k < 20 && !GENES.ready; k++) await new Promise(r => setTimeout(r, 100));
  const l = what === 'hg38' ? parseLocus(START_GENE) : null;
  if (S.tracks.length) return;
  if (l && !l.error) goto(l.id, l.s, l.e); else goto(i, 0, Math.min(S.contigs[i].len, 100000));
  S.viewSet = true;
}
const START_GENE = 'BRCA1';
startup();
