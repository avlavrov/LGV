// LGV: drawing signal tracks (TDF, bigWig) — part of the page script, after app.js's helpers.
// The drawing follows UCSC's wiggle rules (kent/src/hg/hgTracks/wigTrack.c): values are gathered per
// pixel column (bases with data, sum, sum of squares, minimum, maximum); the windowing function picks
// the column's value (maximum: the one of min and max farther from zero; minimum: the nearer one;
// mean); "mean+whiskers" draws the mean dark, mean + 1 SD medium and the maximum light; negative
// values take altColor; the log transform is ln(1 + x) (−ln(1 − x) below zero);
// auto-scale spans the values in view, alwaysZero keeps 0 in the range.
// Two departures: a point is one screen pixel (UCSC: 3 px); points drawn from full-resolution data
// show every value, also when several fall into one pixel column, so a distribution (allele
// frequencies) stays visible.
'use strict';

const SIG_BIG_DOTS = 0.8; // points drawn 2 × 2 px while fewer than this × the track's width in pixels
const SIG_ROWS_MAX = 40, SIG_PREF_KEY = 'lgv.signal.v1';
const WF_UCSC = { maximum: 'max', minimum: 'min', mean: 'mean', 'mean+whiskers': 'whiskers' };
function sigPrefs() { try { const o = JSON.parse(localStorage.getItem(SIG_PREF_KEY) || '{}'); return o && typeof o === 'object' ? o : {}; } catch (e) { return {}; } }
function sigPrefKey(a) { return a.name + '|' + (a.src.blob ? a.src.blob.size : a.src.url || ''); }
function saveSigPref(a) {
  const all = sigPrefs(), u = a.sig.user; const k = sigPrefKey(a);
  if (Object.keys(u).length) all[k] = u; else delete all[k];
  const keys = Object.keys(all); for (const x of keys.slice(0, Math.max(0, keys.length - 100))) delete all[x];
  try { localStorage.setItem(SIG_PREF_KEY, JSON.stringify(all)); } catch (e) { /* kept until reload */ }
}

// Settings: the file's track line first (UCSC defaults for what it leaves out), then LGV's rule for a
// file that declares nothing, then what the user chose in the track's menu (remembered per file).
function sigSettings(a) {
  const p = a.ix.props || {}, d = a.def, u = a.sig ? a.sig.user : (sigPrefs()[sigPrefKey(a)] || {});
  const s = {
    graph: p.visibility === 'dense' && !p.graphType ? 'heatmap' : d.graph,
    wf: WF_UCSC[p.windowing] || 'max',
    auto: p.autoScale !== false, lo: p.viewLimits ? p.viewLimits[0] : 0, hi: p.viewLimits ? p.viewLimits[1] : 1,
    zero: !!p.alwaysZero, log: !!p.log, smooth: p.smoothing || 0,
    yLine: p.yLineOnOff ? (p.yLineMark ?? 0) : null,
    color: p.color ? `rgb(${p.color.join(' ')})` : COL.css.gene, alt: p.altColor ? `rgb(${p.altColor.join(' ')})` : COL.css.muted,
    height: p.height || 0, why: d.why, user: Object.assign({}, u),
  };
  if (d.fixed01 && p.autoScale === undefined && !p.viewLimits) { s.auto = false; s.lo = 0; s.hi = 1; }
  if (['bar', 'points', 'line', 'heatmap'].includes(u.graph)) s.graph = u.graph;
  if (['max', 'mean', 'min', 'whiskers'].includes(u.wf)) s.wf = u.wf;
  if (u.scale === 'auto') s.auto = true; else if (u.scale === '01') { s.auto = false; s.lo = 0; s.hi = 1; } else if (u.scale === 'file' && p.viewLimits) { s.auto = false; s.lo = p.viewLimits[0]; s.hi = p.viewLimits[1]; }
  if (typeof u.zero === 'boolean') s.zero = u.zero;
  if (typeof u.log === 'boolean') s.log = u.log;
  if (u.h >= 12 && u.h <= 200) s.height = u.h;
  return s;
}
async function openSignalTrack(src, key) {
  const f = await openSignal(src);
  const a = { id: ++ANN.seq, key, name: src.name, kind: 'signal', ix: f, src, def: await signalDefaults(f), cache: new WeakMap() };
  a.sig = sigSettings(a);
  f.onData = e => { if (e && !a.err) { a.err = e; toast(t('track.failed', { name: src.name, msg: errText(e) }), true); } draw(); };
  return a;
}
function sigRowH(a) {
  const s = a.sig;
  if (s.height) return s.height;
  if (s.graph === 'heatmap') return 12;
  return a.ix.tracks.length > 3 ? 24 : 40;
}
function sigRows(a) { return Math.min(SIG_ROWS_MAX, a.ix.tracks.length); }
function sigHeight(a) { return ANN_HDR + sigRows(a) * (sigRowH(a) + 1); }

// ---- per pixel column statistics for one panel's view ----
// raw: values one by one; summary: values precomputed per window (zoom levels).
function sigView(a) {
  const f = a.ix, name = S.contigs[S.view.refId].name, vs = S.view.start, bpp = S.view.bpPerPx, ve = vs + W * bpp;
  // values one by one for points, and for a function the file has no summaries of (TDF); the reader
  // falls back to summaries when that would be too much to read for this view
  const wf = a.sig.wf === 'whiskers' ? 'max' : a.sig.wf;
  const needRaw = a.sig.graph === 'points' || (f.format === 'tdf' && wf !== 'mean' && !f.windowFunctions.includes(wf));
  const plan = f.failed ? { failed: true } : f.plan(name, Math.max(0, vs), ve, bpp, refLen(), needRaw ? 'raw' : 'mean');
  if (plan.failed || plan.none || plan.tooMuch || plan.wait) return { plan, cols: null };
  const got = f.blocks(plan, Math.max(0, vs), ve);
  const key = [vs, bpp, W, a.sig.graph, a.sig.wf, a.sig.log, a.sig.smooth, got.blocks.length, got.missing, plan.mode, plan.z ?? plan.zoom].join('|');
  const memo = CUR.sigMemo && CUR.sigMemo.get(a.id);
  if (memo && memo.key === key) return memo.v;
  const nt = sigRows(a), cols = [], pts = [];
  for (let t = 0; t < nt; t++) {
    const C = { n: new Float64Array(W), sum: new Float64Array(W), sq: new Float64Array(W), mn: new Float64Array(W).fill(Infinity), mx: new Float64Array(W).fill(-Infinity) };
    const P = plan.mode === 'raw' && a.sig.graph === 'points' ? [] : null;
    for (const bl of got.blocks) {
      const S0 = bl.start, E0 = bl.end, m = S0.length; if (!m || S0[0] >= ve) continue;
      const raw = !!bl.v, v = raw ? bl.v[t] : null, sm = raw ? null : bl.s;
      const mean = sm && sm.mean ? sm.mean[t] : null, mnA = sm && sm.min ? sm.min[t] : null, mxA = sm && sm.max ? sm.max[t] : null, sdA = sm && sm.sd ? sm.sd[t] : null;
      if (!raw && !mean && !mnA && !mxA) continue;
      // first item that may reach the view (positions are in order within a block)
      let lo = 0, hi = m; while (lo < hi) { const md = (lo + hi) >> 1; if (E0[md] <= vs) lo = md + 1; else hi = md; }
      for (let k = lo; k < m; k++) {
        const s = S0[k], e = E0[k]; if (s >= ve) break; if (e <= vs) continue;
        let x; const val = raw ? v[k] : (mean ? mean[k] : (mxA || mnA)[k]);
        if (Number.isNaN(val)) continue;
        const vmin = raw ? val : (mnA ? mnA[k] : val), vmax = raw ? val : (mxA ? mxA[k] : val);
        if (P) P.push((s + e) / 2, val);
        const x0 = Math.max(0, Math.floor((s - vs) / bpp)), x1 = Math.min(W - 1, Math.floor((e - 1 - vs) / bpp));
        const sd = sdA ? sdA[k] : 0;
        for (x = x0; x <= x1; x++) {
          const cs = vs + x * bpp, ce = cs + bpp, w = Math.max(1e-9, Math.min(e, ce) - Math.max(s, cs));
          C.n[x] += w; C.sum[x] += val * w; C.sq[x] += (val * val + sd * sd) * w;
          if (vmin < C.mn[x]) C.mn[x] = vmin; if (vmax > C.mx[x]) C.mx[x] = vmax;
        }
      }
    }
    cols.push(C); pts.push(P);
  }
  const v = { plan, cols, pts, missing: got.missing };
  (CUR.sigMemo || (CUR.sigMemo = new Map())).set(a.id, { key, v });
  return v;
}
// one image buffer and canvas for painting rows, reused
const SIG_OFF = document.createElement('canvas'); let SIG_IMG = null;
function sigImage(w, h) {
  if (SIG_OFF.width < w || SIG_OFF.height < h) { SIG_OFF.width = Math.max(SIG_OFF.width, w); SIG_OFF.height = Math.max(SIG_OFF.height, h); }
  if (!SIG_IMG || SIG_IMG.width !== w || SIG_IMG.height !== h) SIG_IMG = new ImageData(w, h); else SIG_IMG.data.fill(0);
  return SIG_IMG;
}
const RGB_CACHE = new Map();
function cssRGB(c) {
  let v = RGB_CACHE.get(c); if (v) return v;
  const g = document.createElement('canvas').getContext('2d'); g.fillStyle = '#000'; g.fillStyle = c; g.fillRect(0, 0, 1, 1);
  const p = g.getImageData(0, 0, 1, 1).data; v = [p[0], p[1], p[2]]; RGB_CACHE.set(c, v); return v;
}
const logish = x => x >= 0 ? Math.log1p(x) : -Math.log1p(-x);
// the value a column shows, by the windowing function (UCSC's rules)
function colValue(C, x, wf) {
  if (!C.n[x]) return NaN;
  if (wf === 'max') return Math.abs(C.mn[x]) > Math.abs(C.mx[x]) ? C.mn[x] : C.mx[x];
  if (wf === 'min') return Math.abs(C.mn[x]) < Math.abs(C.mx[x]) ? C.mn[x] : C.mx[x];
  return C.sum[x] / C.n[x];
}

// Draws the rows of a signal track in [y, y + h); returns the text for the name line.
function drawSignal(a, y, h) {
  const s = a.sig, f = a.ix, rh = sigRowH(a), nt = sigRows(a);
  const V = sigView(a), tf = s.log ? logish : (x => x);
  const note = V.plan.failed ? t('sig.failed') : V.plan.none ? t('sig.noContig') : V.plan.tooMuch ? t('sig.tooMuch') : V.plan.wait || (V.missing && !V.cols) ? t('sig.loading') : '';
  let scaleTxt = '';
  for (let r = 0; r < nt; r++) {
    const top = y + r * (rh + 1), C = V.cols && V.cols[r], P = V.pts && V.pts[r];
    if (r > 0) { ctx.strokeStyle = COL.css.line; ctx.beginPath(); ctx.moveTo(0, top - 0.5); ctx.lineTo(W, top - 0.5); ctx.stroke(); }
    if (!C) continue;
    // values per column, smoothed (UCSC: a moving mean over the window, weighted by bases)
    const val = new Float64Array(W).fill(NaN);
    for (let x = 0; x < W; x++) val[x] = tf(colValue(C, x, s.wf === 'whiskers' ? 'mean' : s.wf));
    if (s.smooth > 1) {
      const out = new Float64Array(W).fill(NaN), half = s.smooth >> 1;
      for (let x = 0; x < W; x++) { if (Number.isNaN(val[x])) continue; let sw = 0, sv = 0; for (let k = Math.max(0, x - half); k <= Math.min(W - 1, x + half); k++) if (!Number.isNaN(val[k])) { sw += C.n[k]; sv += val[k] * C.n[k]; } out[x] = sv / sw; }
      val.set(out);
    }
    // the scale
    let lo = Infinity, hi = -Infinity;
    if (s.auto) {
      if (P && P.length) for (let k = 1; k < P.length; k += 2) { const v = tf(P[k]); if (v < lo) lo = v; if (v > hi) hi = v; }
      else for (let x = 0; x < W; x++) { if (Number.isNaN(val[x])) continue; const a1 = s.wf === 'whiskers' ? tf(C.mx[x]) : val[x], a0 = s.wf === 'whiskers' ? tf(C.mn[x]) : val[x]; if (a1 > hi) hi = a1; if (a0 < lo) lo = a0; }
      if (s.zero) { if (hi < 0) hi = 0; else if (lo > 0) lo = 0; }
      if (!(hi > lo)) { if (hi > 0) lo = 0; else if (hi < 0) { lo = hi; hi = 0; } else { lo = -1; hi = 1; } }
    } else { lo = tf(s.lo); hi = tf(s.hi); }
    if (!isFinite(lo) || !isFinite(hi) || !(hi > lo)) { lo = 0; hi = 1; }
    const sc = (+lo.toPrecision(3)) + '–' + (+hi.toPrecision(3)) + (s.log ? ' ln' : '');
    if (nt === 1) scaleTxt = sc;
    // the row is painted into an image, one column per screen pixel and the rows at device pixels,
    // then drawn in one go: the cost does not grow with the number of values or columns drawn
    const R = Math.max(2, Math.round((rh - 2) * DPR)), img = sigImage(W, R), d = img.data, cPos = cssRGB(s.color), cNeg = cssRGB(s.alt);
    const ry = v => Math.round((R - 1) - (Math.max(lo, Math.min(hi, v)) - lo) / (hi - lo) * (R - 1));
    const span = (x, ya, yb, c, al) => {
      const a0 = Math.max(0, Math.min(ya, yb)), a1 = Math.min(R - 1, Math.max(ya, yb)), A = Math.round(al * 255);
      for (let yy = a0, i = (a0 * W + x) * 4; yy <= a1; yy++, i += W * 4) { d[i] = c[0]; d[i + 1] = c[1]; d[i + 2] = c[2]; d[i + 3] = A; }
    };
    const zr = ry(0), px = Math.max(1, Math.round(DPR));
    // a point is one screen pixel; when the points in view are few (fewer than SIG_BIG_DOTS × the
    // track's width), 2 × 2 px, easier on the eyes
    let nIn = 0;
    if (P) { const vs = S.view.start, ve = vs + W * S.view.bpPerPx; for (let k = 0; k < P.length; k += 2) if (P[k] >= vs && P[k] < ve) nIn++; }
    else if (s.graph === 'points') for (let x = 0; x < W; x++) if (!Number.isNaN(val[x])) nIn++;
    const big = nIn < SIG_BIG_DOTS * W, dot = (big ? 2 : 1) * px - 1, dw = big ? 2 : 1;
    const point = (x, y0, c) => { span(x, y0 - (big ? px : 0), y0 - (big ? px : 0) + dot, c, 1); if (dw === 2 && x + 1 < W) span(x + 1, y0 - px, y0 - px + dot, c, 1); };
    if (s.graph === 'heatmap') {
      for (let x = 0; x < W; x++) { const v = val[x]; if (Number.isNaN(v)) continue; span(x, 0, R - 1, v < 0 ? cNeg : cPos, 0.12 + 0.88 * (Math.max(lo, Math.min(hi, v)) - lo) / (hi - lo)); }
    } else if (P) {
      // every value its own point; a spot already marked is skipped (most are, zoomed out)
      const vs = S.view.start, bpp = S.view.bpPerPx, seen = new Uint8Array(W * R), k1 = (R - 1) / (hi - lo);
      for (let k = 0; k < P.length; k += 2) {
        const x = Math.floor((P[k] - vs) / bpp); if (x < 0 || x >= W) continue;
        const v = s.log ? tf(P[k + 1]) : P[k + 1], y0 = Math.round((R - 1) - (Math.max(lo, Math.min(hi, v)) - lo) * k1), cell = y0 * W + x;
        if (seen[cell]) continue; seen[cell] = 1;
        point(x, y0, v < 0 ? cNeg : cPos);
      }
    } else if (s.graph === 'line') {
      let prev = null;
      for (let x = 0; x < W; x++) { const v = val[x]; if (Number.isNaN(v)) { prev = null; continue; } const y0 = ry(v); span(x, prev === null ? y0 : prev, y0, v < 0 ? cNeg : cPos, 1); prev = y0; }
    } else {
      const wh = s.wf === 'whiskers';
      for (let x = 0; x < W; x++) {
        const v = val[x]; if (Number.isNaN(v)) continue;
        const c = v < 0 ? cNeg : cPos;
        if (wh) {
          const mean = C.sum[x] / C.n[x], sd = Math.sqrt(Math.max(0, C.sq[x] / C.n[x] - mean * mean));
          if (s.graph === 'points') { // light: minimum to maximum; medium: mean ± 1 SD; dark: the mean
            span(x, ry(tf(C.mx[x])), ry(tf(C.mn[x])), c, 0.3); span(x, ry(tf(mean + sd)), ry(tf(mean - sd)), c, 0.55); span(x, ry(v), ry(v), c, 1);
          } else { // light: zero to the maximum (minimum below zero); medium: to mean ± 1 SD; dark: to the mean
            span(x, zr, ry(v < 0 ? tf(C.mn[x]) : tf(C.mx[x])), c, 0.3); span(x, zr, ry(tf(v < 0 ? mean - sd : mean + sd)), c, 0.55); span(x, zr, ry(v), c, 1);
          }
        } else if (s.graph === 'points') {
          point(x, ry(v), c);
        } else span(x, zr, ry(v), c, 1);
      }
    }
    if (s.yLine !== null && s.graph !== 'heatmap') { const v = tf(s.yLine); if (v >= lo && v <= hi) { const yl = ry(v), cm = cssRGB(COL.css.muted); for (let x = 0; x < W; x++) if (x % 6 < 3 && !d[(yl * W + x) * 4 + 3]) span(x, yl, yl, cm, 1); } }
    SIG_OFF.getContext('2d').putImageData(img, 0, 0);
    ctx.save(); ctx.imageSmoothingEnabled = false; ctx.drawImage(SIG_OFF, 0, 0, W, R, 0, top + 1, W, rh - 2); ctx.restore();
    ctx.save();
    // the row's name (several tracks) and scale, on a light backing at the top left
    if (nt > 1 || s.graph !== 'heatmap') {
      ctx.font = '8px ' + SANS; ctx.textBaseline = 'top'; ctx.textAlign = 'left';
      const lbl = (nt > 1 ? f.tracks[r] + '  ' : '') + (nt > 1 && s.graph !== 'heatmap' ? sc : '');
      if (lbl) { const tw = ctx.measureText(lbl).width; ctx.globalAlpha = 0.8; ctx.fillStyle = COL.css.panel; ctx.fillRect(2, top + 1, tw + 6, 9); ctx.globalAlpha = 1; ctx.fillStyle = COL.css.muted; ctx.fillText(lbl, 5, top + 2); }
    }
    ctx.restore();
  }
  const extra = [];
  if (scaleTxt) extra.push(scaleTxt);
  if (f.tracks.length > nt) extra.push(t('track.nSamples', { n: fmtInt(f.tracks.length), shown: nt }));
  if (V.plan.mode === 'summary') {
    const win = f.format === 'tdf' ? refLen() / (2 ** V.plan.z) / 700 : V.plan.zoom;
    const fields = f.format === 'tdf' ? Object.keys(V.plan.fields) : ['mean', 'min', 'max'];
    // what is shown when the file has no summary for the chosen function
    const used = s.wf === 'max' && fields.includes('max') ? 'max' : s.wf === 'min' && fields.includes('min') ? 'min' : 'mean';
    extra.push(t('sig.summary', { what: t('sig.wf.' + used), bp: fmtBp(win) }));
  }
  if (note) extra.push(note);
  return extra.join('  ·  ');
}

// row under y in a signal track
function sigRowAt(a, y) { const r = Math.floor((y - a.y0 - ANN_HDR) / (sigRowH(a) + 1)); return r >= 0 && r < sigRows(a) ? r : -1; }
function sigClick(a, x, y, gx) {
  const r = sigRowAt(a, y); S.selected = null; S.hl = null;
  const V = CUR.sigMemo && CUR.sigMemo.get(a.id); const C = V && V.v.cols && V.v.cols[r];
  if (r < 0 || !C) { hideCard(); draw(); return; }
  const xi = Math.max(0, Math.min(W - 1, Math.floor(x))), bpp = S.view.bpPerPx;
  const pos = Math.floor(S.view.start + x * bpp), attrs = {};
  // the value at the position itself when the data are one value per position (zoomed in)
  let near = null;
  const P = V.v.pts && V.v.pts[r];
  if (P) { let best = 6 * bpp; for (let k = 0; k < P.length; k += 2) { const d = Math.abs(P[k] - pos); if (d < best) { best = d; near = { at: Math.floor(P[k]), v: P[k + 1] }; } } }
  if (!C.n[xi] && !near) { hideCard(); draw(); return; }
  const fmt = v => String(+v.toPrecision(6));
  const f = { start: Math.floor(S.view.start + xi * bpp), end: Math.max(Math.floor(S.view.start + (xi + 1) * bpp), Math.floor(S.view.start + xi * bpp) + 1), name: a.ix.tracks[r], attrs };
  // one value: the card is about its position; several in the pixel column: their statistics too
  const several = C.n[xi] && (V.v.plan.mode === 'summary' || C.mx[xi] > C.mn[xi] || C.n[xi] > (near ? 1 : 0) && !near);
  if (near) { f.start = near.at; f.end = near.at + 1; f.value = fmt(near.v); }
  if (several) {
    if (near) attrs[t('sig.column')] = `${fmtInt(Math.floor(S.view.start + xi * bpp) + 1)}-${fmtInt(Math.floor(S.view.start + (xi + 1) * bpp))}`;
    attrs[t('sig.mean')] = fmt(C.sum[xi] / C.n[xi]); attrs[t('sig.min')] = fmt(C.mn[xi]); attrs[t('sig.max')] = fmt(C.mx[xi]);
    if (V.v.plan.mode === 'raw') attrs[t('sig.bases')] = fmtInt(Math.round(C.n[xi]));
  } else if (!near) f.value = fmt(C.sum[xi] / C.n[xi]);
  showAnnCard(a, f, gx, y); draw();
}

// ---- the track's menu (right click): how to draw it ----
function openSignalMenu(a, gx, y) {
  const m = $('ctx'), s = a.sig, u = s.user, p = a.ix.props || {};
  const item = (act, val, on, label) => `<button class="mi" role="menuitemradio" aria-checked="${on}" data-sig="${act}" data-val="${val}"><span class="mk" aria-hidden="true">${on ? '●' : ''}</span>${esc(label)}</button>`;
  const chk = (act, on, label) => `<button class="mi" role="menuitemcheckbox" aria-checked="${on}" data-sig="${act}" data-val=""><span class="mk" aria-hidden="true">${on ? '✓' : ''}</span>${esc(label)}</button>`;
  const scale = s.auto ? 'auto' : s.lo === 0 && s.hi === 1 ? '01' : 'file';
  menuState = { pn: CUR, sig: a };
  m.innerHTML = `<div class="ms">${esc(t('sig.graph'))}</div>` + ['bar', 'points', 'line', 'heatmap'].map(g => item('graph', g, s.graph === g, t('sig.g.' + g))).join('') +
    `<div class="ms">${esc(t('sig.windowing'))}</div>` + ['max', 'mean', 'min', 'whiskers'].map(w => item('wf', w, s.wf === w, t('sig.wf.' + w))).join('') +
    `<div class="ms">${esc(t('sig.scale'))}</div>` + item('scale', 'auto', scale === 'auto', t('sig.s.auto')) + item('scale', '01', scale === '01', t('sig.s.01')) +
    (p.viewLimits ? item('scale', 'file', scale === 'file', t('sig.s.file', { lo: p.viewLimits[0], hi: p.viewLimits[1] })) : '') +
    chk('zero', s.zero, t('sig.zero')) + chk('log', s.log, t('sig.log')) +
    `<div class="ms">${esc(t('sig.height'))}</div>` + item('h', 'less', false, t('sig.h.less')) + item('h', 'more', false, t('sig.h.more')) +
    (Object.keys(u).length ? '<div class="msep"></div>' + `<button class="mi" role="menuitem" data-sig="reset" data-val=""><span class="mk" aria-hidden="true"></span>${esc(t('sig.reset'))}</button>` : '');
  m.hidden = false;
  const r = stage.getBoundingClientRect(), w = m.offsetWidth, hh = m.offsetHeight;
  m.style.left = Math.max(8, Math.min(r.width - w - 8, gx + 2)) + 'px';
  m.style.top = Math.max(8, Math.min(r.height - hh - 8, y + 2)) + 'px';
  (m.querySelector('[aria-checked="true"]') || m.querySelector('.mi')).focus({ preventScroll: true });
}
function signalMenuAct(a, act, val) {
  const u = a.sig.user;
  if (act === 'graph') u.graph = val;
  else if (act === 'wf') u.wf = val;
  else if (act === 'scale') u.scale = val;
  else if (act === 'zero') u.zero = !a.sig.zero;
  else if (act === 'log') u.log = !a.sig.log;
  else if (act === 'h') u.h = Math.max(12, Math.min(200, sigRowH(a) + (val === 'more' ? 12 : -12)));
  else if (act === 'reset') for (const k of Object.keys(u)) delete u[k];
  a.sig = sigSettings(a); a.sig.user = u;
  saveSigPref(a);
  for (const pn of S.panels) if (pn.sigMemo) pn.sigMemo.delete(a.id);
  clampScroll(); draw();
}
