// Checks of the signal readers (TDF, bigWig) on real files made by igvtools, UCSC tools and bigtools
// (test/signal, from the IGV and igv.js test sets, MIT license), against what the files say about
// themselves and against the text files they were made from; then fuzzing of both readers.
// Run: node tools/test-signal.js [rounds]
'use strict';
const fs = require('fs'), path = require('path');
const C = require('../src/core.js'); Object.assign(globalThis, C);
const S = require('../src/signal.js');
const D = path.join(__dirname, '..', 'test', 'signal');
const ROUNDS = +process.argv[2] || 300;
let fails = 0;
const ok = (cond, msg) => { console.log((cond ? 'OK   ' : 'FAIL ') + msg); if (!cond) fails++; };
const src = (u8, name) => new C.FileSource(new Blob([u8]), name);
const file = name => src(fs.readFileSync(path.join(D, name)), name);

// Everything a view needs, waiting for the background reads; returns { plan, blocks }.
async function fetchView(f, contig, a, b, bpp, chrLen, want) {
  for (let i = 0; i < 1000; i++) {
    if (f.failed) throw f.failed;
    const plan = f.plan(contig, a, b, bpp, chrLen, want);
    if (plan.wait) { await idle(f); continue; }
    if (plan.none || plan.tooMuch) return { plan, blocks: [] };
    const r = f.blocks(plan, a, b);
    if (!r.missing) return { plan, blocks: r.blocks };
    await idle(f);
  }
  throw new Error('never finished');
}
// all raw values of a contig: in windows small enough to be read as raw, without repeats
async function rawAll(f, contig, len, t) {
  let w = len, plan;
  for (;;) { ({ plan } = await fetchView(f, contig, 0, w, 1, len, 'raw')); if (plan.mode === 'raw' || plan.none || w < 1000) break; w = Math.ceil(w / 4); }
  if (plan.mode !== 'raw') return { mode: plan.mode, rows: [] };
  const seen = new Map();
  for (let a = 0; a < len; a += w) {
    const { blocks } = await fetchView(f, contig, a, Math.min(len, a + w), 1, len, 'raw');
    for (const bl of blocks) for (let k = 0; k < bl.start.length; k++) for (let tt = 0; tt < bl.v.length; tt++) { if (t !== undefined && tt !== t) continue; const v = bl.v[tt][k]; if (!Number.isNaN(v)) seen.set(bl.start[k] + '|' + bl.end[k] + '|' + tt, [bl.start[k], bl.end[k], v]); }
  }
  return { mode: 'raw', rows: [...seen.values()].sort((x, y) => x[0] - y[0]) };
}
function idle(f) { return new Promise(res => { const tick = () => { if (f.failed || !f.loader.busy) res(); else setTimeout(tick, 1); }; tick(); }); }
// raw values of track t as [start, end, value] in order
function rawOf(blocks, t, a, b) {
  const out = [];
  for (const bl of blocks) for (let k = 0; k < bl.start.length; k++) { const v = bl.v[t][k]; if (!Number.isNaN(v) && bl.end[k] > a && bl.start[k] < b) out.push([bl.start[k], bl.end[k], v]); }
  return out.sort((x, y) => x[0] - y[0]);
}

(async () => {
  // ---------- bigWig ----------
  // 1. fixedStep.bw against the wiggle file it was made from
  {
    const f = await S.openSignal(file('fixedStep.bw'));
    const wig = fs.readFileSync(path.join(D, 'fixedStep.wig'), 'utf8').split('\n');
    const exp = []; let pos = 0, step = 1, span = 1;
    for (const l of wig) {
      const m = l.match(/^fixedStep chrom=\S+ start=(\d+) step=(\d+)(?: span=(\d+))?/);
      if (m) { pos = +m[1] - 1; step = +m[2]; span = +(m[3] || 1); continue; }
      if (l.trim() && isFinite(+l)) { exp.push([pos, pos + span, +l]); pos += step; }
    }
    const { plan, blocks } = await fetchView(f, 'chr1', 0, 248956422, 1, 248956422, 'raw');
    const got = rawOf(blocks, 0, 0, 248956422);
    ok(plan.mode === 'raw' && got.length === exp.length && got.every((g, i) => g[0] === exp[i][0] && g[1] === exp[i][1] && Math.abs(g[2] - exp[i][2]) < 1e-6),
      `fixedStep.bw: ${got.length} values match fixedStep.wig (${exp.length})`);
  }
  // 2. every bigWig: the values read agree with the file's own total summary (count, min, max, mean),
  //    and its coarsest zoom level agrees with the values (min, max, mean)
  for (const name of ['fixedStep.bw', 'bigEndian.bw', 'uncompressed.bw', 'chromTreeByteOrder.bw', 'bigtools.bigWig', 'nozooms.bw']) {
    const f = await S.openSignal(file(name));
    let n = 0, mn = Infinity, mx = -Infinity, sum = 0;
    const zs = [];
    for (const c of f.chrs.values()) {
      const { blocks } = await fetchView(f, c.name, 0, c.size, 1, c.size, 'raw');
      for (const [s, e, v] of rawOf(blocks, 0, 0, c.size)) { const w = e - s; n += w; sum += v * w; if (v < mn) mn = v; if (v > mx) mx = v; }
      if (f.zooms.length) { const r = await fetchView(f, c.name, 0, c.size, f.zooms[f.zooms.length - 1].reduction * 2, c.size); zs.push(r); }
    }
    const T = f.total, close = (x, y) => Math.abs(x - y) <= 1e-4 * Math.max(1, Math.abs(y));
    ok(T && T.n === n && close(T.min, mn) && close(T.max, mx) && close(T.mean, sum / n), `${name} (${f.le ? 'little' : 'big'}-endian, ${f.uncompressBufSize ? 'zlib' : 'uncompressed'}): ${n} bases, min ${mn}, max ${+mx.toPrecision(6)}, mean ${+(sum / n).toPrecision(6)} agree with the total summary`);
    if (zs.length) {
      let zn = 0, zmn = Infinity, zmx = -Infinity, zsum = 0;
      for (const r of zs) for (const bl of r.blocks) for (let k = 0; k < bl.start.length; k++) { const c = bl.s.n[0][k]; zn += c; zsum += bl.s.mean[0][k] * c; zmn = Math.min(zmn, bl.s.min[0][k]); zmx = Math.max(zmx, bl.s.max[0][k]); }
      ok(zs.every(r => r.plan.mode === 'summary') && zn === n && close(zmn, mn) && close(zmx, mx) && close(zsum / zn, sum / n), `${name}: coarsest zoom level (${f.zooms[f.zooms.length - 1].reduction} bp) agrees with the values`);
    }
  }
  // ---------- TDF ----------
  // 3. wholeGenomeTest.tdf against the bedGraph it was made from (raw datasets)
  {
    const f = await S.openSignal(file('wholeGenomeTest.tdf'));
    const exp = new Map();
    for (const l of fs.readFileSync(path.join(D, 'wholeGenomeTest.bedgraph'), 'utf8').split('\n')) {
      const p = l.split('\t'); if (p.length < 4 || l[0] === '#') continue;
      if (!exp.has(p[0])) exp.set(p[0], []); exp.get(p[0]).push([+p[1], +p[2], +p[3]]);
    }
    // igvtools stored the long bedGraph intervals cut into windows: every window must lie inside an
    // interval of the same value, and the windows must cover the intervals
    let all = 0, good = 0;
    for (const [chr, rows] of exp) {
      const end = Math.max(...rows.map(r => r[1]));
      const { mode, rows: got } = await rawAll(f, chr, end, 0); const plan = { mode };
      let covered = 0, inside = got.length > 0;
      for (const [s, e, v] of got) { const r = rows.find(r => r[0] <= s && e <= r[1]); if (!r || r[2] !== v) inside = false; covered += e - s; }
      const want = rows.reduce((t, r) => t + r[1] - r[0], 0);
      all++; if (plan.mode === 'raw' && inside && covered === want) good++;
    }
    ok(good === all, `wholeGenomeTest.tdf: raw values of ${good}/${all} contigs lie in the bedGraph intervals of the same value and cover them`);
  }
  // 4. every TDF: max and min over the raw values of all tracks agree with the file's root statistics
  //    ("Maximum", "Minimum"), and the zoom-0 means lie within them; all tracks are read
  for (const name of fs.readdirSync(D).filter(n => n.endsWith('.tdf'))) {
    const f = await S.openSignal(file(name));
    let mn = Infinity, mx = -Infinity, n = 0, zOut = 0, zN = 0;
    for (const c of f.chrs.values()) {
      if (/^(chr)?(m|mt)$/i.test(c)) continue; // igvtools leaves the mitochondrial genome out of its statistics
      const { mode, rows } = await rawAll(f, c, 3e8);
      if (mode !== 'raw') continue;
      for (const [, , x] of rows) { n++; if (x < mn) mn = x; if (x > mx) mx = x; }
    }
    for (const c of f.chrs.values()) {
      if (/^(chr)?(m|mt)$/i.test(c)) continue;
      const r = await fetchView(f, c, 0, 3e8, 3e8 / 700, 3e8);
      for (const bl of r.blocks) for (const v of bl.s.mean) for (const x of v) if (!Number.isNaN(x)) { zN++; if (x < mn - 1e-3 || x > mx + 1e-3) zOut++; }
    }
    const R = f.root, close = (x, y) => Math.abs(x - y) <= 1e-3 * Math.max(1, Math.abs(y));
    const stat = R.Maximum !== undefined ? close(mx, +R.Maximum) && close(mn, +R.Minimum) : true;
    ok(n > 0 && stat && zN > 0 && !zOut, `${name}: ${f.tracks.length} track(s), ${n} raw values, min ${+mn.toPrecision(6)} max ${+mx.toPrecision(6)}${R.Maximum !== undefined ? ` = root ${R.Minimum}…${R.Maximum}` : ''}; ${zN} zoom-0 means inside`);
  }
  // 5. the track line
  {
    const p = S.parseTrackLine('track type=wiggle_0 name="My track" graphType=points viewLimits=0:1 autoScale=off alwaysZero=on maxHeightPixels=64:40:11 yLineMark=0.5 yLineOnOff=on windowingFunction=mean+whiskers smoothingWindow=4 transformFunc=LOG color=255,0,0 altColor=0,0,255 junk=1 viewLimits2=x');
    ok(p.name === 'My track' && p.graphType === 'points' && p.viewLimits.join() === '0,1' && p.autoScale === false && p.alwaysZero === true && p.height === 40 && p.yLineMark === 0.5 && p.yLineOnOff && p.windowing === 'mean+whiskers' && p.smoothing === 4 && p.log === true && p.color.join() === '255,0,0' && p.altColor.join() === '0,0,255' && !('junk' in p),
      'track line: every UCSC setting read, unknown keys ignored');
    const bad = S.parseTrackLine('track viewLimits=5:1 color=300,0,0 graphType=pie maxHeightPixels=1:1:1 smoothingWindow=40 __proto__=x');
    ok(!bad.viewLimits && !bad.color && !bad.graphType && !bad.height && bad.smoothing === undefined && Object.getPrototypeOf(bad) === Object.prototype, 'track line: values out of range ignored');
  }
  // 6. the default way to draw: allele frequencies (0–1 at separate positions) as points, coverage as bars
  {
    const cov = await S.signalDefaults(await S.openSignal(file('BRAF_C.bam.tdf')));
    const af = await S.signalDefaults(await S.openSignal(file('af_sample.tdf')));
    const dense = await S.signalDefaults(await S.openSignal(file('dense01_sample.tdf')));
    ok(cov.graph === 'bar' && af.graph === 'points' && af.fixed01 && dense.graph === 'bar', `default drawing: coverage ${cov.graph}, allele frequencies ${af.graph} (0–1 scale), dense 0–1 values ${dense.graph}`);
  }

  // ---------- fuzzing: damaged files end in a clear error or in data, never in anything else ----------
  let seed = 4242;
  const rnd = n => { seed ^= seed << 13; seed >>>= 0; seed ^= seed >>> 17; seed ^= seed << 5; seed >>>= 0; return seed % n; };
  for (const name of ['BRAF_C.bam.tdf', 'recombRate.igv.txt.tdf', 'fixedStep.bw', 'bigEndian.bw', 'uncompressed.bw', 'af_sample.tdf']) {
    const base = fs.readFileSync(path.join(D, name));
    const codes = {}, bad = []; let fine = 0, slow = 0; const t0 = Date.now();
    for (let i = 0; i < ROUNDS; i++) {
      const u = Uint8Array.from(base), kind = rnd(4);
      if (kind === 0) for (let k = 0; k < 1 + rnd(8); k++) u[rnd(u.length)] = rnd(256);
      else if (kind === 1) { const p = rnd(Math.min(u.length, 4096)); for (let k = 0; k < 4; k++) u[p + k] = 255; }
      else if (kind === 2) { const p = rnd(u.length); u.fill(0, p, Math.min(u.length, p + 1 + rnd(64))); }
      const data = kind === 3 ? u.subarray(0, rnd(u.length)) : u, t = Date.now();
      try {
        const f = await S.openSignal(src(data, name));
        for (const c of [...f.chrs.values()].slice(0, 3)) {
          const nm = c.name || c, len = c.size || 3e8;
          for (const bpp of [1, 100, len / 700]) for (const want of ['raw', 'mean']) await fetchView(f, nm, 0, len, bpp, len, want);
        }
        await S.signalDefaults(f);
        fine++;
      } catch (e) {
        if (e instanceof C.CgvError) codes[e.code] = (codes[e.code] || 0) + 1;
        else bad.push(`#${i}: ${e && e.constructor && e.constructor.name}: ${String(e && e.message).slice(0, 120)}`);
      }
      if (Date.now() - t > 5000) slow++;
    }
    ok(!bad.length && !slow, `fuzz ${name}: ${ROUNDS} cases, ${fine} read, clear errors ${JSON.stringify(codes)}, unexpected ${bad.length}, slow ${slow} (${Date.now() - t0} ms)`);
    for (const b of bad.slice(0, 5)) console.log('       ' + b);
  }
  console.log(fails ? `\n${fails} check(s) failed` : '\nAll signal checks passed');
  process.exit(fails ? 1 : 0);
})();
