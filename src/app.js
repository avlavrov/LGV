// CGV app: state, loading, WebGL read rendering, Canvas2D for ruler/coverage, UI.
'use strict';

// ======================= Display profile =======================
const DEFAULT_PROFILE = {
  version: 1,
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
  maxReads: 300000,
  maxLoadBp: 3000000,
  hideToolbar: false,
  language: 'en',                // en | ru
};
const PROFILE_KEY = 'cgv.profile.v1';
function clone(o) { return JSON.parse(JSON.stringify(o)); }
function mergeProfile(base, over) {
  const out = clone(base);
  for (const k of Object.keys(over || {})) {
    if (!(k in base)) continue;
    if (base[k] && typeof base[k] === 'object' && !Array.isArray(base[k])) out[k] = Object.assign({}, base[k], over[k] && typeof over[k] === 'object' ? over[k] : {});
    else if (typeof over[k] === typeof base[k]) out[k] = over[k];
  }
  return out;
}
let P = DEFAULT_PROFILE;
try { const s = localStorage.getItem(PROFILE_KEY); if (s) P = mergeProfile(DEFAULT_PROFILE, JSON.parse(s)); } catch (e) { P = clone(DEFAULT_PROFILE); }
if (P === DEFAULT_PROFILE) P = clone(DEFAULT_PROFILE);
function saveProfile() { try { localStorage.setItem(PROFILE_KEY, JSON.stringify(P)); } catch (e) { /* storage unavailable: the profile lives until reload */ } }

// ======================= State =======================
const S = {
  bam: null, fasta: null, bamName: '', isDemo: false,
  view: { refId: 0, start: 0, bpPerPx: 1 },
  data: null,          // { refId, start, end, reads, refArr, insert }
  shown: null,         // { reads, cov, groups, inst }
  scroll: 0, sortPos: null, selected: null,
  loading: false, loadToken: 0, forceLoad: false,
  stats: { reads: 0, loadMs: 0, buildMs: 0, drawMs: 0, bytes: 0 },
  msg: null,
};

// ======================= DOM =======================
const $ = id => document.getElementById(id);
const stage = $('stage'), glCanvas = $('gl'), ov = $('ov'), ctx = ov.getContext('2d');
let W = 0, H = 0, DPR = 1;

// track geometry (CSS px)
const RULER_H = 26, REF_H = 14;
function covTop() { return RULER_H + REF_H + 4; }
function readsTop() { return covTop() + P.coverage.height + 6; }
function groupGap() { return P.groupBy === 'none' ? 0 : 20; }

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
    css: { bg: g('--bg'), fg: g('--fg'), muted: g('--muted'), line: g('--line'), cov: g('--cov'), accent: g('--accent'), panel: g('--panel') },
    read: hexToRgb(g('--read')), fwd: hexToRgb(g('--read-fwd')), rev: hexToRgb(g('--read-rev')),
    del: hexToRgb(g('--del')), skip: hexToRgb(g('--muted')), ins: hexToRgb(g('--ins')), clip: hexToRgb(g('--clip')),
    big: hexToRgb(g('--pair-large')), small: hexToRgb(g('--pair-small')), rf: hexToRgb(g('--pair-rf')), ff: hexToRgb(g('--pair-ff')), mu: hexToRgb(g('--pair-unmapped')),
    base: { 65: hexToRgb(g('--base-a')), 67: hexToRgb(g('--base-c')), 71: hexToRgb(g('--base-g')), 84: hexToRgb(g('--base-t')), 78: hexToRgb(g('--muted')) },
    cat: [1, 2, 3, 4, 5, 6, 7, 8].map(i => hexToRgb(g('--cat-' + i))),
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
  uniform float u_viewStart, u_bpPerPx, u_rowH, u_gap, u_grpGap, u_scroll, u_top, u_mmFade;
  uniform vec2 u_size;
  out vec4 v_color;
  void main() {
    // minpx < 0: an isolated mismatch that fades when zoomed out
    float minpx = abs(a_y.z);
    float fade = a_y.z < 0.0 ? clamp(u_mmFade / u_bpPerPx, 0.0, 1.0) : 1.0;
    float px0 = (a_x.x - u_viewStart) / u_bpPerPx;
    float px1 = (a_x.y - u_viewStart) / u_bpPerPx;
    if (px1 - px0 < minpx) { float c = (px0 + px1) * 0.5; px0 = c - minpx * 0.5; px1 = c + minpx * 0.5; }
    float h = max(u_rowH - u_gap, 1.0);
    float y0 = u_top + a_x.z * u_rowH + a_x.w * u_grpGap - u_scroll;
    float ya = y0 + h * a_y.x, yb = y0 + h * a_y.y;
    if (yb - ya < 1.0) yb = ya + 1.0;
    float x = mix(px0, px1, a_corner.x), y = mix(ya, yb, a_corner.y);
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
  GL.u = {}; for (const n of ['u_viewStart', 'u_bpPerPx', 'u_rowH', 'u_gap', 'u_grpGap', 'u_scroll', 'u_top', 'u_size', 'u_mmFade']) GL.u[n] = gl.getUniformLocation(prog, n);
  GL.vao = gl.createVertexArray(); gl.bindVertexArray(GL.vao);
  const quad = gl.createBuffer(); gl.bindBuffer(gl.ARRAY_BUFFER, quad);
  gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([0, 0, 1, 0, 0, 1, 1, 1]), gl.STATIC_DRAW);
  gl.enableVertexAttribArray(0); gl.vertexAttribPointer(0, 2, gl.FLOAT, false, 0, 0);
  GL.fbuf = gl.createBuffer(); GL.cbuf = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, GL.fbuf);
  gl.enableVertexAttribArray(1); gl.vertexAttribPointer(1, 4, gl.FLOAT, false, 28, 0); gl.vertexAttribDivisor(1, 1);
  gl.enableVertexAttribArray(2); gl.vertexAttribPointer(2, 3, gl.FLOAT, false, 28, 16); gl.vertexAttribDivisor(2, 1);
  gl.bindBuffer(gl.ARRAY_BUFFER, GL.cbuf);
  gl.enableVertexAttribArray(3); gl.vertexAttribPointer(3, 4, gl.UNSIGNED_BYTE, true, 0, 0); gl.vertexAttribDivisor(3, 1);
  gl.enable(gl.BLEND); gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
  GL.n = 0;
  return true;
}
const GL_OK = (() => { try { return initGL(); } catch (e) { console.error(e); return false; } })();

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

function readColor(rec) {
  switch (P.colorBy) {
    case 'strand': return rec.flag & F_REV ? COL.rev : COL.fwd;
    case 'tag': { const v = recTags(rec)[P.colorTag]; return v === undefined ? COL.read : COL.cat[hashStr(v) % 8]; }
    case 'pair': {
      const o = pairOrientation(rec), ins = S.data.insert;
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

function refAtFactory() { const d = S.data; return p => { const j = p - d.start; return j >= 0 && j < d.refArr.length ? d.refArr[j] : 0; }; }
function feats(rec) { if (!rec.feats) rec.feats = computeFeatures(rec, refAtFactory()); return rec.feats; }

// Full rebuild: filters → downsampling → coverage → groups → layout → sort → instances.
function rebuild() {
  const d = S.data; if (!d) { S.shown = null; GL.n = 0; return; }
  const t0 = performance.now();
  let reads = d.reads.filter(r => passFilters(r, P.filters));
  const covReads = reads;
  let sampled = false;
  if (reads.length > P.maxReads) { // even downsampling with a fixed step, deterministic
    const step = reads.length / P.maxReads, out = [];
    for (let i = 0; i < P.maxReads; i++) out.push(reads[Math.floor(i * step)]);
    reads = out; sampled = true;
  }
  const cov = computeCoverage(covReads, d.start, d.end - d.start);
  // groups
  const gmap = new Map();
  for (const r of reads) { const k = groupKey(r); if (!gmap.has(k)) gmap.set(k, []); gmap.get(k).push(r); }
  const keys = [...gmap.keys()].sort((a, b) => a.localeCompare(b, LANG));
  const groups = []; let rowCursor = 0;
  for (const k of keys) {
    const rows = packRows(gmap.get(k), 3);
    sortRows(rows);
    groups.push({ key: k, rows, start: rowCursor, count: gmap.get(k).length });
    rowCursor += rows.length;
  }
  // instances: bodies → lines → features on top
  const ib = new InstBuf(Math.max(1024, reads.length * 3));
  const bodies = [], marks = [];
  groups.forEach((g, gi) => g.rows.forEach((row, ri) => { const R = g.start + ri; for (const rec of row) { bodies.push([rec, R, gi]); } }));
  for (const [rec, R, gi] of bodies) {
    const col = readColor(rec), a = mapqAlpha(rec.mapq);
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
  }
  for (const [t, x0, x1, R, gi, a] of marks) ib.add(x0 - d.start, x1 - d.start, R, gi, t === 2 ? 0.44 : 0.48, t === 2 ? 0.56 : 0.52, 0, t === 2 ? COL.del : COL.skip, a);
  for (const [rec, R, gi] of bodies) {
    const a = mapqAlpha(rec.mapq);
    for (const f of feats(rec)) {
      if (f.t === 0) ib.add(f.pos - d.start, f.pos + 1 - d.start, R, gi, 0, 1, flaggedPos(cov, d, f.pos) ? 1 : -1, COL.base[f.b] || COL.base[78], a * qualAlpha(f.q));
      else if (f.t === 1) ib.add(f.pos - d.start, f.pos - d.start, R, gi, -0.08, 1.08, 2, COL.ins, Math.max(a, 0.6));
      else if ((f.t === 4 || f.t === 5) && P.showSoftClips && f.len >= 1) ib.add(f.pos - d.start, f.pos - d.start, R, gi, 0, 1, 1.5, COL.clip, Math.max(a, 0.7));
    }
  }
  S.shown = { reads, cov, groups, totalRows: rowCursor, sampled, total: covReads.length };
  if (GL_OK) {
    gl.bindBuffer(gl.ARRAY_BUFFER, GL.fbuf); gl.bufferData(gl.ARRAY_BUFFER, ib.f.subarray(0, ib.n * 7), gl.STATIC_DRAW);
    gl.bindBuffer(gl.ARRAY_BUFFER, GL.cbuf); gl.bufferData(gl.ARRAY_BUFFER, ib.c.subarray(0, ib.n * 4), gl.STATIC_DRAW);
    GL.n = ib.n;
  }
  S.stats.buildMs = performance.now() - t0;
  S.stats.reads = covReads.length;
  clampScroll();
}

// A position where non-reference alleles pass the threshold (same rule as the coverage highlight).
function flaggedPos(cov, d, pos) {
  const j = (pos - d.start) * 6, refc = d.refArr[pos - d.start];
  if (j < 0 || j >= cov.length || !refc) return false;
  const t = cov[j] + cov[j + 1] + cov[j + 2] + cov[j + 3] + cov[j + 4]; if (t < 2) return false;
  const ri = { 65: 0, 67: 1, 71: 2, 84: 3 }[refc]; if (ri === undefined) return false;
  return (t - cov[j + ri]) / t >= P.coverage.alleleThreshold;
}

function sortRows(rows) {
  if (P.sortBy === 'start') return;
  const pos = S.sortPos != null ? S.sortPos : Math.floor(viewCenter());
  const refAt = refAtFactory(), rb = refAt(pos);
  const keyOf = row => {
    const rec = row.find(r => r.pos <= pos && r.end > pos);
    if (!rec) return [9, 0];
    switch (P.sortBy) {
      case 'base': { const b = baseAtRefPos(rec, pos); if (!b) return [8, 0]; if (b === 45) return [2, 0]; return b === rb ? [3, 0] : [1, b]; }
      case 'mapq': return [1, -rec.mapq];
      case 'strand': return [1, rec.flag & F_REV ? 1 : 0];
      case 'insert': return [1, -Math.abs(rec.tlen)];
    }
    return [1, 0];
  };
  const keyed = rows.map(r => [keyOf(r), r]);
  keyed.sort((a, b) => a[0][0] - b[0][0] || a[0][1] - b[0][1]);
  for (let i = 0; i < rows.length; i++) rows[i] = keyed[i][1];
}

// ======================= View =======================
function refLen() { return S.bam ? S.bam.refs[S.view.refId].len : 0; }
function viewEnd() { return S.view.start + W * S.view.bpPerPx; }
function viewCenter() { return S.view.start + W * S.view.bpPerPx / 2; }
function clampView() {
  const L = refLen(); if (!L || !W) return;
  const minBpp = 1 / 48, maxBpp = Math.max(minBpp, L / W);
  S.view.bpPerPx = Math.min(maxBpp, Math.max(minBpp, S.view.bpPerPx));
  const span = W * S.view.bpPerPx;
  S.view.start = Math.min(L - span * 0.5, Math.max(-span * 0.5, S.view.start));
}
function contentHeight() { if (!S.shown) return 0; return S.shown.totalRows * P.readHeight + S.shown.groups.length * groupGap(); }
function clampScroll() { const maxS = Math.max(0, contentHeight() - (H - readsTop()) + 24); S.scroll = Math.min(maxS, Math.max(0, S.scroll)); }
function goto(refId, start, end) {
  S.view.refId = refId;
  const span = Math.max(20, end - start);
  S.view.bpPerPx = span / Math.max(1, W);
  S.view.start = start;
  clampView(); S.scroll = 0; requestLoad(); draw(); syncLocus();
}
function zoomAt(px, factor) {
  const bp = S.view.start + px * S.view.bpPerPx;
  S.view.bpPerPx *= factor; clampView();
  S.view.start = bp - px * S.view.bpPerPx; clampView();
  requestLoad(); draw(); syncLocus();
}
function pan(dpx) { S.view.start += dpx * S.view.bpPerPx; clampView(); requestLoad(); draw(); syncLocus(); }

// ======================= Loading =======================
let loadTimer = 0;
function requestLoad() { clearTimeout(loadTimer); loadTimer = setTimeout(ensureData, 120); }
async function ensureData() {
  if (!S.bam || !W) return;
  const vs = Math.max(0, Math.floor(S.view.start)), ve = Math.min(refLen(), Math.ceil(viewEnd()));
  const d = S.data;
  if (d && d.refId === S.view.refId && vs >= d.start && ve <= d.end) return;
  const span = ve - vs;
  if (span > P.maxLoadBp && !S.forceLoad) {
    S.data = null; S.shown = null; GL.n = 0;
    S.msg = { text: t('wide', { span: fmtBp(span) }), action: t('wide.btn', { span: fmtBp(span) }) };
    draw(); return;
  }
  S.forceLoad = false;
  const ls = Math.max(0, vs - Math.min(span, 50000)), le = Math.min(refLen(), ve + Math.min(span, 50000));
  const token = ++S.loadToken; S.loading = true; S.msg = null; draw();
  const t0 = performance.now(); const b0 = S.bam.src.bytes;
  try {
    const reads = await S.bam.fetch(S.view.refId, ls, le);
    let refArr = null;
    const refName = S.bam.refs[S.view.refId].name;
    if (S.fasta) { try { refArr = await S.fasta.fetch(refName, ls, le); } catch (e) { refArr = null; } }
    if (token !== S.loadToken) return;
    if (!refArr || refArr.length !== le - ls) { refArr = new Uint8Array(le - ls); for (const r of reads) fillRefFromMD(r, refArr, ls); S.refFromMD = true; } else S.refFromMD = false;
    S.data = { refId: S.view.refId, start: ls, end: le, reads, refArr, insert: insertSizeBounds(reads) };
    S.stats.loadMs = performance.now() - t0; S.stats.bytes = S.bam.src.bytes - b0;
    S.loading = false; S.selected = null; hideCard();
    rebuild(); draw();
  } catch (e) {
    if (token !== S.loadToken) return;
    S.loading = false; S.msg = { text: errText(e), error: true }; draw();
  }
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
function fmtInt(n) { return Math.round(n).toLocaleString(NUM_LOCALE[LANG]); }

function render() {
  const t0 = performance.now();
  ctx.setTransform(DPR, 0, 0, DPR, 0, 0);
  ctx.clearRect(0, 0, W, H);
  const rt = readsTop();
  // ---- reads (WebGL) ----
  if (GL_OK) {
    gl.viewport(0, 0, glCanvas.width, glCanvas.height);
    gl.clearColor(0, 0, 0, 0); gl.clear(gl.COLOR_BUFFER_BIT);
    if (S.shown && GL.n) {
      gl.enable(gl.SCISSOR_TEST);
      gl.scissor(0, 0, glCanvas.width, Math.max(0, Math.floor((H - rt) * DPR)));
      gl.useProgram(GL.prog); gl.bindVertexArray(GL.vao);
      gl.uniform1f(GL.u.u_viewStart, S.view.start - S.data.start);
      gl.uniform1f(GL.u.u_bpPerPx, S.view.bpPerPx);
      gl.uniform1f(GL.u.u_rowH, P.readHeight);
      gl.uniform1f(GL.u.u_gap, P.readHeight >= 5 ? 1 : 0);
      gl.uniform1f(GL.u.u_grpGap, groupGap());
      gl.uniform1f(GL.u.u_scroll, S.scroll);
      gl.uniform1f(GL.u.u_top, rt + groupGap());
      gl.uniform2f(GL.u.u_size, W, H);
      gl.uniform1f(GL.u.u_mmFade, 1.5);
      gl.drawArraysInstanced(gl.TRIANGLE_STRIP, 0, 4, GL.n);
      gl.disable(gl.SCISSOR_TEST);
    }
  }
  drawRuler();
  drawRefTrack();
  drawCoverage();
  drawReadsOverlay(rt);
  drawMessage();
  S.stats.drawMs = performance.now() - t0;
  updateStatus();
}

function niceStep(raw) { const p = Math.pow(10, Math.floor(Math.log10(raw))); for (const m of [1, 2, 5, 10]) if (m * p >= raw) return m * p; return 10 * p; }

function drawRuler() {
  if (!S.bam) return;
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

function drawRefTrack() {
  const d = S.data; if (!d) return;
  const bpp = S.view.bpPerPx, ppb = 1 / bpp, y = RULER_H + 2;
  if (ppb < 1.2) return;
  const b0 = Math.max(d.start, Math.floor(S.view.start)), b1 = Math.min(d.end, Math.ceil(viewEnd()));
  ctx.textAlign = 'center'; ctx.font = '600 11px ' + MONO;
  for (let b = b0; b < b1; b++) {
    const c = d.refArr[b - d.start]; if (!c) continue;
    const x = (b - S.view.start) * ppb;
    if (ppb >= 8) { ctx.fillStyle = COL.baseCss[c] || COL.css.muted; ctx.fillText(String.fromCharCode(c), x + ppb / 2, y + 10); }
    else { ctx.fillStyle = COL.baseCss[c] || COL.css.muted; ctx.fillRect(x, y + 2, Math.max(1, ppb - (ppb > 3 ? 1 : 0)), REF_H - 6); }
  }
}

function drawCoverage() {
  const sh = S.shown, d = S.data; if (!sh || !d) return;
  const top = covTop(), h = P.coverage.height, bpp = S.view.bpPerPx, ppb = 1 / bpp;
  const cov = sh.cov, thr = P.coverage.alleleThreshold;
  const b0 = Math.max(d.start, Math.floor(S.view.start)), b1 = Math.min(d.end, Math.ceil(viewEnd()));
  let maxD = 1;
  for (let b = b0; b < b1; b++) { const j = (b - d.start) * 6; const t = cov[j] + cov[j + 1] + cov[j + 2] + cov[j + 3] + cov[j + 4] + cov[j + 5]; if (t > maxD) maxD = t; }
  const scale = (h - 10) / maxD, base = top + h;
  const flagged = (j) => flaggedPos(cov, d, d.start + j / 6);
  const drawStack = (x, w, j) => {
    let y = base;
    for (const [k, c] of [[0, 65], [1, 67], [2, 71], [3, 84]]) { const hh = cov[j + k] * scale; if (hh > 0) { ctx.fillStyle = COL.baseCss[c]; ctx.fillRect(x, y - hh, w, hh); y -= hh; } }
  };
  if (ppb >= 1) {
    for (let b = b0; b < b1; b++) {
      const j = (b - d.start) * 6, x = (b - S.view.start) * ppb, w = Math.max(1, ppb - (ppb > 4 ? 1 : 0));
      if (flagged(j)) drawStack(x, w, j);
      else { const t = cov[j] + cov[j + 1] + cov[j + 2] + cov[j + 3] + cov[j + 4]; ctx.fillStyle = COL.css.cov; ctx.fillRect(x, base - t * scale, w, t * scale); }
    }
  } else {
    for (let px = 0; px < W; px++) {
      const s = Math.floor(S.view.start + px * bpp), e = Math.floor(S.view.start + (px + 1) * bpp);
      let mx = 0, fj = -1;
      for (let b = Math.max(s, d.start); b < Math.min(Math.max(e, s + 1), d.end); b++) {
        const j = (b - d.start) * 6, t = cov[j] + cov[j + 1] + cov[j + 2] + cov[j + 3] + cov[j + 4];
        if (t > mx) mx = t;
        if (fj < 0 && flagged(j)) fj = j;
      }
      if (fj >= 0) drawStack(px, 1, fj); else if (mx) { ctx.fillStyle = COL.css.cov; ctx.fillRect(px, base - mx * scale, 1, mx * scale); }
    }
  }
  ctx.fillStyle = COL.css.muted; ctx.font = '10.5px ' + MONO; ctx.textAlign = 'left';
  ctx.fillText('0–' + fmtInt(maxD) + '×', 6, top + 10);
  ctx.strokeStyle = COL.css.line; ctx.beginPath(); ctx.moveTo(0, base + 3.5); ctx.lineTo(W, base + 3.5); ctx.stroke();
}

function drawReadsOverlay(rt) {
  const sh = S.shown; if (!sh) return;
  ctx.save(); ctx.beginPath(); ctx.rect(0, rt, W, H - rt); ctx.clip();
  const rh = P.readHeight, gg = groupGap(), bpp = S.view.bpPerPx, ppb = 1 / bpp, d = S.data;
  // group labels
  if (gg) sh.groups.forEach((g, gi) => {
    const y = rt + g.start * rh + gi * gg - S.scroll + gg;
    ctx.strokeStyle = COL.css.line; ctx.beginPath(); ctx.moveTo(0, Math.round(y - gg + 3) + 0.5); ctx.lineTo(W, Math.round(y - gg + 3) + 0.5); ctx.stroke();
    ctx.font = '600 11px ' + SANS; ctx.textAlign = 'left'; ctx.fillStyle = COL.css.fg;
    ctx.fillText(g.key + '  ·  ' + fmtInt(g.count), 8, y - 6);
  });
  // mismatch letters and insertion lengths when zoomed in
  if (ppb >= 7 && rh >= 9) {
    ctx.font = '600 ' + Math.min(12, rh - 1) + 'px ' + MONO; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
    const vs = S.view.start, ve = viewEnd();
    sh.groups.forEach((g, gi) => {
      const firstRow = Math.max(0, Math.floor((S.scroll - gi * gg) / rh) - g.start - 1);
      const lastRow = Math.min(g.rows.length - 1, Math.ceil((S.scroll + H - rt - gi * gg) / rh) - g.start + 1);
      for (let ri = firstRow; ri <= lastRow; ri++) {
        const y = rt + gg + (g.start + ri) * rh + gi * gg - S.scroll + (rh - (rh >= 5 ? 1 : 0)) / 2;
        for (const rec of g.rows[ri]) {
          if (rec.end < vs || rec.pos > ve) continue;
          const a = mapqAlpha(rec.mapq);
          for (const f of feats(rec)) {
            if (f.t === 0) { ctx.globalAlpha = a * qualAlpha(f.q); ctx.fillStyle = '#fff'; ctx.fillText(String.fromCharCode(f.b), (f.pos - vs + 0.5) * ppb, y + 0.5); }
            else if (f.t === 1 && ppb >= 10) { ctx.globalAlpha = 1; ctx.fillStyle = COL.css.fg; ctx.font = '600 9px ' + MONO; ctx.fillText(String(f.len), (f.pos - vs) * ppb, y - rh / 2 - 1); ctx.font = '600 ' + Math.min(12, rh - 1) + 'px ' + MONO; }
          }
        }
      }
    });
    ctx.globalAlpha = 1; ctx.textBaseline = 'alphabetic';
  }
  // selected read
  if (S.selected) {
    const L = S.selected;
    const y = rt + gg + L.row * rh + L.gi * gg - S.scroll;
    const x0 = (L.rec.pos - S.view.start) * ppb, x1 = (L.rec.end - S.view.start) * ppb;
    ctx.strokeStyle = COL.css.accent; ctx.lineWidth = 2; ctx.strokeRect(x0 - 1, y - 1, x1 - x0 + 2, Math.max(2, rh - (rh >= 5 ? 1 : 0)) + 2); ctx.lineWidth = 1;
  }
  ctx.restore();
}

function drawMessage() {
  const el = $('msg');
  if (S.loading) { el.hidden = false; el.className = 'msg'; el.innerHTML = '<span class="spin"></span>' + esc(t('loading')); return; }
  if (!S.bam) { el.hidden = true; return; }
  if (S.msg) {
    el.hidden = false; el.className = 'msg' + (S.msg.error ? ' err' : '');
    el.textContent = S.msg.text;
    if (S.msg.action) { const b = document.createElement('button'); b.className = 'btn'; b.textContent = S.msg.action; b.onclick = () => { S.forceLoad = true; S.msg = null; ensureData(); }; el.append(' ', b); }
    return;
  }
  el.hidden = true;
}

function updateStatus() {
  const st = S.stats, sh = S.shown;
  let s = '';
  if (sh) s = t('status', { reads: fmtInt(st.reads), sampled: sh.sampled ? t('status.sampled', { n: fmtInt(sh.reads.length) }) : '', load: Math.round(st.loadMs), build: Math.round(st.buildMs), frame: st.drawMs.toFixed(1) });
  if (S.refFromMD && sh) s += ' · ' + t('ref.fromMd');
  $('status').textContent = s;
}

// ======================= Locus =======================
function syncLocus() {
  if (!S.bam || document.activeElement === $('locus')) return;
  const s = Math.max(0, Math.round(S.view.start)), e = Math.min(refLen(), Math.round(viewEnd()));
  $('locus').value = `${S.bam.refs[S.view.refId].name}:${fmtInt(s + 1)}-${fmtInt(e)}`;
  $('span').textContent = fmtBp(e - s);
}
function parseLocus(str) {
  str = str.trim().replace(/[  \s]/g, '').replace(/[–—]/g, '-');
  const m = str.match(/^([^:]+?)(?::([\d,.]+)(?:-([\d,.]+))?)?$/);
  if (!m) return null;
  const id = S.bam.refId(m[1]); if (id < 0) return { error: t('locus.noContig', { name: m[1] }) };
  const num = x => parseInt(x.replace(/[,.]/g, ''), 10);
  const L = S.bam.refs[id].len;
  if (!m[2]) return { id, s: 0, e: L };
  const a = num(m[2]) - 1;
  if (!m[3]) return { id, s: a - 50, e: a + 51, center: a };
  return { id, s: a, e: num(m[3]) };
}

// ======================= Input =======================
let drag = null;
// Genome gestures only start on the canvas itself: panels, cards and messages on top of
// the stage keep their own clicks, scrolling and text selection.
const onCanvas = e => e.target === ov;
stage.addEventListener('pointerdown', e => {
  if (e.button !== 0 || !S.bam || !onCanvas(e)) return;
  stage.setPointerCapture(e.pointerId);
  drag = { x: e.clientX, y: e.clientY, sx: S.view.start, sc: S.scroll, moved: false, inReads: e.offsetY > readsTop() };
});
stage.addEventListener('pointermove', e => {
  if (!drag) return;
  const dx = e.clientX - drag.x, dy = e.clientY - drag.y;
  if (Math.abs(dx) + Math.abs(dy) > 3) drag.moved = true;
  if (!drag.moved) return;
  S.view.start = drag.sx - dx * S.view.bpPerPx; clampView();
  if (drag.inReads) { S.scroll = drag.sc - dy; clampScroll(); }
  requestLoad(); draw(); syncLocus();
});
stage.addEventListener('pointerup', e => {
  if (!drag) return;
  const d = drag; drag = null;
  if (!d.moved) onClick(e.offsetX, e.offsetY);
});
stage.addEventListener('dblclick', e => { if (S.bam && onCanvas(e)) zoomAt(e.offsetX, 0.5); });
stage.addEventListener('wheel', e => {
  if (!S.bam || !onCanvas(e)) return;
  e.preventDefault();
  const dy = e.deltaMode === 1 ? e.deltaY * 16 : e.deltaY, dx = e.deltaMode === 1 ? e.deltaX * 16 : e.deltaX;
  if (e.ctrlKey || e.metaKey) { zoomAt(e.offsetX, Math.exp(dy * 0.01)); return; }
  if (e.shiftKey || Math.abs(dx) > Math.abs(dy)) { pan(e.shiftKey ? dy : dx); return; }
  S.scroll += dy; clampScroll(); draw();
}, { passive: false });

function hitTest(x, y) {
  const sh = S.shown; if (!sh) return null;
  const rt = readsTop(), rh = P.readHeight, gg = groupGap();
  if (y < rt) return null;
  for (let gi = 0; gi < sh.groups.length; gi++) {
    const g = sh.groups[gi];
    const gy = rt + gg + g.start * rh + gi * gg - S.scroll;
    const ri = Math.floor((y - gy) / rh);
    if (ri < 0 || ri >= g.rows.length) continue;
    const bp = S.view.start + x * S.view.bpPerPx, slop = 2 * S.view.bpPerPx;
    const row = g.rows[ri];
    let lo = 0, hi = row.length - 1, best = null;
    while (lo <= hi) { const m = (lo + hi) >> 1; if (row[m].pos <= bp + slop) { best = m; lo = m + 1; } else hi = m - 1; }
    for (let k = best; k != null && k >= 0 && k >= best - 3; k--) if (row[k].end + slop > bp) return { rec: row[k], row: g.start + ri, gi };
  }
  return null;
}

function onClick(x, y) {
  const h = hitTest(x, y);
  if (!h) { S.selected = null; hideCard(); draw(); return; }
  S.selected = h; showCard(h.rec, x, y); draw();
}

const FLAG_BITS = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048];
function esc(s) { return String(s).replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c])); }
function cardText(rec) {
  const refs = S.bam.refs, tags = recTags(rec);
  const lines = [
    [t('card.name'), rec.name], [t('card.pos'), `${refs[rec.refID].name}:${fmtInt(rec.pos + 1)}-${fmtInt(rec.end)}`], ['CIGAR', cigarString(rec)],
    ['MAPQ', rec.mapq], [t('card.flag'), `${rec.flag} (${FLAG_BITS.filter(b => rec.flag & b).map(b => t('flag.' + b)).join(', ') || '—'})`],
    [t('card.len'), rec.lseq + ' ' + t('unit.bp')],
  ];
  if (rec.flag & 1) lines.push([t('card.mate'), rec.mateRef >= 0 ? `${refs[rec.mateRef].name}:${fmtInt(rec.matePos + 1)}` : t('card.mateUnmapped')], [t('card.tlen'), rec.tlen]);
  for (const k of Object.keys(tags)) lines.push([k, Array.isArray(tags[k]) ? tags[k].join(',') : tags[k]]);
  return lines;
}
function showCard(rec, x, y) {
  const c = $('card');
  const rows = cardText(rec).map(([k, v]) => `<dt>${esc(k)}</dt><dd>${esc(v)}</dd>`).join('');
  const mate = (rec.flag & 1) && rec.mateRef >= 0;
  c.innerHTML = `<dl>${rows}</dl><div class="card-actions">${mate ? `<button class="btn" id="toMate">${esc(t('card.toMate'))}</button>` : ''}<button class="btn" id="copyRead">${esc(t('card.copy'))}</button><button class="btn ghost" id="closeCard">${esc(t('btn.close'))}</button></div>`;
  c.hidden = false;
  const r = stage.getBoundingClientRect(), cw = c.offsetWidth, ch = c.offsetHeight;
  c.style.left = Math.max(8, Math.min(r.width - cw - 8, x + 14)) + 'px';
  c.style.top = Math.max(8, Math.min(r.height - ch - 8, y + 14)) + 'px';
  $('closeCard').onclick = () => { S.selected = null; hideCard(); draw(); };
  $('copyRead').onclick = () => {
    const txt = cardText(rec).map(([k, v]) => k + '\t' + v).join('\n');
    navigator.clipboard.writeText(txt).then(() => toast(t('copied')), () => toast(t('copyFail')));
  };
  if (mate) $('toMate').onclick = () => { hideCard(); S.selected = null; goto(rec.mateRef, rec.matePos - 150, rec.matePos + 300); };
}
function hideCard() { $('card').hidden = true; }

document.addEventListener('keydown', e => {
  if (e.target.matches('input, textarea, select')) { if (e.key === 'Escape') e.target.blur(); return; }
  if (!S.bam) return;
  if (e.key === 'ArrowLeft') pan(-W * 0.15);
  else if (e.key === 'ArrowRight') pan(W * 0.15);
  else if (e.key === 'ArrowUp') { S.scroll -= 60; clampScroll(); draw(); }
  else if (e.key === 'ArrowDown') { S.scroll += 60; clampScroll(); draw(); }
  else if (e.key === '+' || e.key === '=') zoomAt(W / 2, 0.5);
  else if (e.key === '-' || e.key === '_') zoomAt(W / 2, 2);
  else if (e.key === 's' || e.key === 'ы') sortHere(); // 'ы' / 'р': same keys on a Russian layout
  else if (e.key === 'h' || e.key === 'р') { P.hideToolbar = !P.hideToolbar; applyChrome(); saveProfile(); }
  else if (e.key === 'Escape') { hideCard(); closePanels(); S.selected = null; draw(); }
  else return;
  e.preventDefault();
});

function sortHere() {
  S.sortPos = Math.floor(viewCenter());
  if (P.sortBy === 'start') { P.sortBy = 'base'; syncPanel(); saveProfile(); }
  rebuild(); draw();
}

// ======================= Size =======================
function resize() {
  const r = stage.getBoundingClientRect();
  DPR = Math.min(window.devicePixelRatio || 1, 3);
  const oldW = W; W = Math.max(1, r.width); H = Math.max(1, r.height);
  for (const c of [glCanvas, ov]) { c.width = Math.round(W * DPR); c.height = Math.round(H * DPR); }
  if (oldW && S.bam) { const c = S.view.start + oldW * S.view.bpPerPx / 2; S.view.start = c - W * S.view.bpPerPx / 2; clampView(); }
  clampScroll(); draw(); syncLocus();
}
new ResizeObserver(resize).observe(stage);

// ======================= Opening data =======================
async function openBam(bamSrc, baiSrc, faSrc, faiSrc, opts) {
  opts = opts || {};
  S.loadToken++; S.loading = true; S.msg = null; draw();
  try {
    const bam = await new BamFile(bamSrc, baiSrc).open();
    let fasta = null;
    if (faSrc && faiSrc) { try { fasta = await new IndexedFasta(faSrc, faiSrc, opts.gzi).open(); } catch (e) { toast(t('ref.failed', { msg: errText(e) })); } }
    const keep = opts.keepRef && S.fasta && bam.refs.some(r => S.fasta.has(r.name));
    S.bam = bam; if (fasta || !keep) S.fasta = fasta; updateRefInfo(); S.data = null; S.shown = null; GL.n = 0; S.isDemo = !!opts.demo; S.sortPos = null;
    S.bamName = bamSrc.name; $('fileName').textContent = bamSrc.name + (bam.samples.length ? ' · ' + bam.samples.join(', ') : '');
    $('examples').hidden = !S.isDemo;
    S.loading = false;
    const loc = opts.locus ? parseLocus(opts.locus) : null;
    if (loc && !loc.error) goto(loc.id, loc.s, loc.e);
    else { const L = bam.refs[0].len; goto(0, 0, Math.min(L, 20000)); }
  } catch (e) {
    S.loading = false; S.msg = { text: errText(e), error: true }; draw();
  }
}

function b64ToBlob(b64) { const bin = atob(b64); const u = new Uint8Array(bin.length); for (let i = 0; i < bin.length; i++) u[i] = bin.charCodeAt(i); return new Blob([u]); }
function openDemo() {
  if (typeof DEMO === 'undefined') return;
  openBam(new FileSource(b64ToBlob(DEMO.bam), 'demo.bam'), new FileSource(b64ToBlob(DEMO.bai), 'demo.bam.bai'),
    new FileSource(b64ToBlob(DEMO.fa), 'demo.fa'), new FileSource(b64ToBlob(DEMO.fai), 'demo.fa.fai'),
    { demo: true, locus: 'demo:11,941-12,060' });
}

function openLocalFiles(files) {
  const by = ext => [...files].find(f => f.name.toLowerCase().endsWith(ext));
  const cram = by('.cram');
  if (cram) { toast(t('files.cram')); return; }
  const bam = by('.bam'), bai = by('.bai');
  const ref = pickRefFiles(files);
  if (!bam && ref) { openReference(ref); return; }
  if (by('.csi') && !bai) { toast(t('files.csi')); return; }
  if (!bam || !bai) { toast(t('files.needPair')); return; }
  closePanels();
  openBam(new FileSource(bam), new FileSource(bai), null, null, { keepRef: !ref }).then(() => { if (ref) openReference(ref); });
}

function openUrls() {
  const u = $('urlBam').value.trim(); if (!u) { toast(t('url.need')); return; }
  if (/\.cram(\?|$)/i.test(u)) { toast(t('files.cram')); return; }
  const idx = $('urlBai').value.trim() || addSuffix(u, '.bai');
  const auth = $('urlAuth').value.trim();
  const headers = auth ? { Authorization: /^bearer\s/i.test(auth) ? auth : 'Bearer ' + auth } : {};
  openBam(new UrlSource(u, headers), new UrlSource(idx, headers), null, null, { keepRef: true });
  closePanels();
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
async function openReference(ref) {
  if (!S.bam) { toast(t('ref.needBam')); return; }
  if (refBusy) return;
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
    if (!S.bam.refs.some(r => fasta.has(r.name))) { toast(t('ref.mismatch')); return; }
    S.fasta = fasta; S.data = null; S.shown = null; GL.n = 0; updateRefInfo(); closePanels();
    toast(fasta.built ? t('ref.indexed', { name, n: fmtInt(fasta.idx.size) }) : t('ref.loaded', { name }));
    ensureData();
  } catch (e) { toast(t('ref.failed', { msg: errText(e) })); }
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

// ======================= Panels and settings =======================
function closePanels() { $('panel').hidden = true; $('openPanel').hidden = true; $('viewBtn').setAttribute('aria-expanded', 'false'); $('openBtn').setAttribute('aria-expanded', 'false'); }
function togglePanel(id, btn) { const el = $(id), was = el.hidden; closePanels(); el.hidden = !was; btn.setAttribute('aria-expanded', String(was)); if (id === 'panel' && was) syncPanel(); }
$('viewBtn').onclick = () => togglePanel('panel', $('viewBtn'));
$('openBtn').onclick = () => togglePanel('openPanel', $('openBtn'));
document.querySelectorAll('[data-close]').forEach(b => b.onclick = closePanels);
$('pickFiles').onchange = e => { if (e.target.files.length) openLocalFiles(e.target.files); e.target.value = ''; };
$('openUrl').onclick = openUrls;
$('openDemo').onclick = () => { openDemo(); closePanels(); };
stage.addEventListener('dragover', e => { e.preventDefault(); stage.classList.add('drop'); });
stage.addEventListener('dragleave', () => stage.classList.remove('drop'));
stage.addEventListener('drop', e => { e.preventDefault(); stage.classList.remove('drop'); if (e.dataTransfer.files.length) openLocalFiles(e.dataTransfer.files); });

$('locus').addEventListener('keydown', e => {
  if (e.key !== 'Enter' || !S.bam) return;
  const l = parseLocus(e.target.value);
  if (!l) { toast(t('locus.bad')); return; }
  if (l.error) { toast(l.error); return; }
  e.target.blur(); goto(l.id, l.s, l.e);
});
$('examples').onchange = e => { const v = e.target.value; e.target.value = ''; if (!v || !S.bam) return; const l = parseLocus(v); if (l && !l.error) goto(l.id, l.s, l.e); };

// profile controls: id → [profile path, type, what to recompute]
const CONTROLS = [
  ['readHeight', 'readHeight', 'num', 'draw'],
  ['colorBy', 'colorBy', 'str', 'rebuild'], ['colorTag', 'colorTag', 'str', 'rebuild'],
  ['groupBy', 'groupBy', 'str', 'rebuild'], ['groupTag', 'groupTag', 'str', 'rebuild'],
  ['sortBy', 'sortBy', 'str', 'rebuild'],
  ['mqLow', 'mapqAlpha.low', 'num', 'rebuild'], ['mqHigh', 'mapqAlpha.high', 'num', 'rebuild'], ['mqMinAlpha', 'mapqAlpha.minAlpha', 'num', 'rebuild'],
  ['bqShade', 'baseQualityShading', 'bool', 'rebuild'], ['softClips', 'showSoftClips', 'bool', 'rebuild'],
  ['minMapq', 'filters.minMapq', 'num', 'rebuild'], ['hideDup', 'filters.hideDuplicates', 'bool', 'rebuild'],
  ['hideSec', 'filters.hideSecondary', 'bool', 'rebuild'], ['hideSupp', 'filters.hideSupplementary', 'bool', 'rebuild'], ['hideQc', 'filters.hideQcFail', 'bool', 'rebuild'],
  ['covHeight', 'coverage.height', 'num', 'draw'], ['alleleThr', 'coverage.alleleThreshold', 'num', 'rebuild'],
  ['hideToolbar', 'hideToolbar', 'bool', 'chrome'], ['language', 'language', 'str', 'lang'],
];
const getPath = (o, p) => p.split('.').reduce((a, k) => a[k], o);
const setPath = (o, p, v) => { const ks = p.split('.'), last = ks.pop(); ks.reduce((a, k) => a[k], o)[last] = v; };
function syncPanel() {
  for (const [id, path, type] of CONTROLS) { const el = $(id); if (!el) continue; const v = getPath(P, path); if (type === 'bool') el.checked = v; else el.value = v; }
  $('readHeightOut').textContent = (+P.readHeight).toFixed(1) + ' px';
  $('mqMinAlphaOut').textContent = Math.round(P.mapqAlpha.minAlpha * 100) + '%';
  $('alleleThrOut').textContent = Math.round(P.coverage.alleleThreshold * 100) + '%';
  $('covHeightOut').textContent = P.coverage.height + ' px';
  $('colorTagRow').hidden = P.colorBy !== 'tag'; $('groupTagRow').hidden = P.groupBy !== 'tag';
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
    else { clearTimeout(rebuildTimer); rebuildTimer = setTimeout(() => { if (id === 'sortBy' && S.sortPos == null) S.sortPos = Math.floor(viewCenter()); rebuild(); draw(); }, 30); }
  });
}
$('sortHere').onclick = sortHere;
$('copyProfile').onclick = () => navigator.clipboard.writeText(JSON.stringify(P, null, 2)).then(() => toast(t('profile.copied')), () => { $('profileJson').select(); toast(t('copyFail')); });
$('applyProfile').onclick = () => {
  try { P = mergeProfile(DEFAULT_PROFILE, JSON.parse($('profileJson').value)); saveProfile(); applyLanguage(); syncPanel(); applyChrome(); rebuild(); draw(); toast(t('profile.applied')); }
  catch (e) { toast(t('profile.bad', { msg: e.message })); }
};
$('resetProfile').onclick = () => { P = clone(DEFAULT_PROFILE); saveProfile(); applyLanguage(); syncPanel(); applyChrome(); rebuild(); draw(); toast(t('profile.reset')); };
$('showToolbar').onclick = () => { P.hideToolbar = false; saveProfile(); applyChrome(); };
function applyChrome() { document.body.classList.toggle('bare', P.hideToolbar); $('showToolbar').hidden = !P.hideToolbar; requestAnimationFrame(resize); }

let toastTimer = 0;
function toast(msg, sticky) { const el = $('toast'); el.textContent = msg; el.hidden = false; clearTimeout(toastTimer); if (!sticky) toastTimer = setTimeout(() => el.hidden = true, msg.length > 90 ? 7000 : 3500); }

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
openDemo();
