// LGV: signal tracks — numeric values along the genome, read by region with range reads:
//  - TDF, the tiled format written by igvtools. There is no published specification; the layout here
//    is the one IGV itself reads and writes (header, master index, datasets, groups, tiles).
//  - bigWig, the UCSC format (Kent et al., Bioinformatics 2010; kent/src/inc/bbiFile.h): header, zoom
//    levels of summaries, chromosome B+ tree, R tree index, data sections; either byte order.
// How a track is drawn follows its UCSC track line (wiggle_0 settings) where the file has one: a TDF
// keeps the line of the wiggle or bedGraph it was made from; a bigWig has none of its own.
//
// Safety (see SECURITY.md): every count, offset and size is checked against the bytes actually read
// before use; blocks are capped before and after decompression; tree walks are capped in depth and
// nodes; one view never reads more than a fixed number of bytes of full-resolution data, less over
// a network (coarser data is used instead), so a large file cannot keep the page busy.
'use strict';
if (typeof contigKey === 'undefined' && typeof require !== 'undefined') Object.assign(globalThis, require('./core.js'));
if (typeof bedFeature === 'undefined' && typeof require !== 'undefined') Object.assign(globalThis, require('./tracks.js'));

const SIG_LIMITS = {
  header: 16 << 20, index: 64 << 20, block: 64 << 20, string: 100000, items: 1e8, depth: 24, nodes: 100000,
  chroms: 1000000, viewBytes: 16 << 20, viewBytesUrl: 4 << 20, viewBlocks: 3000, rawTiles: 3000, cacheBytes: 192 << 20, queue: 64, parallel: 6, tracks: 10000,
};
Object.assign(CGV_ERR_EN, {
  tdfBad: 'The TDF file is damaged or not TDF: {what}',
  bwBad: 'The bigWig file is damaged or not bigWig: {what}',
  bwBigBed: '{name} is a bigBed file (features), not bigWig (values)',
  bbBad: 'the bigBed file is damaged: {what}',
});

// ---- binary reading with bounds checks ----
class SigBin {
  constructor(u8, le, code) { this.u8 = u8; this.dv = new DataView(u8.buffer, u8.byteOffset, u8.byteLength); this.p = 0; this.le = le !== false; this.code = code || 'tdfBad'; }
  bad(what) { return new CgvError(this.code, { what }); }
  need(n) { if (!(n >= 0) || this.p + n > this.u8.length) throw this.bad('cut short'); }
  get left() { return this.u8.length - this.p; }
  u8v() { this.need(1); return this.u8[this.p++]; }
  u16() { this.need(2); const v = this.dv.getUint16(this.p, this.le); this.p += 2; return v; }
  i32() { this.need(4); const v = this.dv.getInt32(this.p, this.le); this.p += 4; return v; }
  u32() { this.need(4); const v = this.dv.getUint32(this.p, this.le); this.p += 4; return v; }
  f32() { this.need(4); const v = this.dv.getFloat32(this.p, this.le); this.p += 4; return v; }
  f64() { this.need(8); const v = this.dv.getFloat64(this.p, this.le); this.p += 8; return v; }
  // 64-bit offset or count as a Number; -1 for a negative (TDF: an empty tile); larger than 2^53 is damage
  i64() {
    this.need(8);
    const a = this.dv.getUint32(this.p, this.le), b = this.dv.getUint32(this.p + 4, this.le); this.p += 8;
    const lo = this.le ? a : b, hi = this.le ? b : a;
    if (hi & 0x80000000) return -1;
    if (hi > 0x1fffff) throw this.bad('an offset is too large');
    return hi * 4294967296 + lo;
  }
  str() {
    let e = this.p; const L = this.u8.length, max = Math.min(L, this.p + SIG_LIMITS.string);
    while (e < max && this.u8[e] !== 0) e++;
    if (e >= max) throw this.bad('a text field does not end');
    const s = new TextDecoder().decode(this.u8.subarray(this.p, e)); this.p = e + 1; return s;
  }
  fixedStr(n) { this.need(n); let e = this.p; while (e < this.p + n && this.u8[e] !== 0) e++; const s = new TextDecoder().decode(this.u8.subarray(this.p, e)); this.p += n; return s; }
  // a count of items of `each` bytes that must fit in what is left
  count(each, max, what, v) { const n = v === undefined ? this.i32() : v; if (!(n >= 0) || n > max || (each && n * each > this.left)) throw this.bad(what || 'a count is out of range'); return n; }
}

// ---- UCSC track line (wiggle_0 settings; https://genome.ucsc.edu/goldenpath/help/wiggle.html and
// bigWig.html). Only known keys are kept, each checked; anything else is ignored. ----
const onOff = v => /^on$/i.test(v) ? true : /^off$/i.test(v) ? false : undefined;
function parseTrackLine(line) {
  const out = {}, raw = {};
  let s = String(line || '').slice(0, 20000).trim();
  if (/^track\b/i.test(s)) s = s.slice(5);
  const re = /([A-Za-z_][\w.]*)\s*=\s*(?:"([^"]*)"|'([^']*)'|(\S+))/g;
  let m, n = 0;
  while ((m = re.exec(s)) && n++ < 200) raw[m[1]] = (m[2] ?? m[3] ?? m[4] ?? '').slice(0, 1000);
  const num = v => (v !== undefined && v !== '' && isFinite(+v) ? +v : undefined);
  const range = v => { const p = String(v || '').split(':'); const a = num(p[0]), b = num(p[1]); return p.length === 2 && a !== undefined && b !== undefined && b > a ? [a, b] : undefined; };
  const rgb = v => { if (!/^\d{1,3},\d{1,3},\d{1,3}$/.test(v || '')) return undefined; const c = v.split(',').map(Number); return c.every(x => x <= 255) ? c : undefined; };
  if (raw.name) out.name = raw.name;
  if (raw.description) out.description = raw.description;
  const g = String(raw.graphType || '').toLowerCase();
  if (['bar', 'points', 'line', 'heatmap'].includes(g)) out.graphType = g; // UCSC: bar, points; line and heatmap as IGV adds them
  if (/^(full|dense|hide|pack|squish)$/i.test(raw.visibility || '')) out.visibility = raw.visibility.toLowerCase();
  const vl = range(raw.viewLimits); if (vl) out.viewLimits = vl;
  for (const k of ['autoScale', 'alwaysZero', 'gridDefault', 'yLineOnOff']) { const v = onOff(raw[k]); if (v !== undefined) out[k] = v; }
  if (num(raw.yLineMark) !== undefined) out.yLineMark = num(raw.yLineMark);
  if (raw.maxHeightPixels) { const p = raw.maxHeightPixels.split(':').map(Number); const d = p.length === 3 ? p[1] : p[0]; if (d >= 4 && d <= 1000) out.height = Math.round(d); }
  const wf = String(raw.windowingFunction || '').toLowerCase();
  if (['mean+whiskers', 'maximum', 'mean', 'minimum'].includes(wf)) out.windowing = wf;
  const sw = String(raw.smoothingWindow || '').toLowerCase();
  if (sw === 'off') out.smoothing = 0; else if (num(sw) >= 2 && num(sw) <= 16) out.smoothing = Math.round(num(sw));
  const tf = String(raw.transformFunc || '').toUpperCase(); if (tf === 'LOG' || tf === 'NONE') out.log = tf === 'LOG';
  const c = rgb(raw.color); if (c) out.color = c;
  const a = rgb(raw.altColor); if (a) out.altColor = a;
  return out;
}

// ---- a cache of decoded blocks shared by all signal files, limited by bytes ----
const SIG_CACHE = { map: new Map(), bytes: 0 };
function sigCacheGet(k) { const v = SIG_CACHE.map.get(k); if (v) { SIG_CACHE.map.delete(k); SIG_CACHE.map.set(k, v); } return v; }
function sigCachePut(k, v, bytes) {
  SIG_CACHE.map.set(k, v); v._bytes = bytes; SIG_CACHE.bytes += bytes;
  for (const [key, old] of SIG_CACHE.map) { if (SIG_CACHE.bytes <= SIG_LIMITS.cacheBytes) break; SIG_CACHE.map.delete(key); SIG_CACHE.bytes -= old._bytes; }
}

// A decoded block: positions and per-track values, all typed arrays.
//   raw:     { start, end, v: [Float32Array per track], names? }
//   summary: { start, end, s: { mean?, min?, max?, sd?, n? } each [Float32Array per track] }
const blockBytes = b => b.start.byteLength * 2 + (b.v ? b.v.length * b.start.length * 4 : Object.keys(b.s).length * b.start.length * 4 * (b.ntr || 1));

// Background loading shared by both formats: at most a few reads at a time; when the view moves on,
// requests still waiting are dropped (they are asked for again if still needed).
class SigLoader {
  constructor(file) { this.file = file; this.queue = []; this.active = 0; this.pending = new Set(); }
  // keys: what the job fills (one key, or the blocks of one read)
  want(keys, job) {
    if (typeof keys === 'string') keys = [keys];
    if (this.file.failed || keys.some(k => this.pending.has(k))) return;
    for (const k of keys) this.pending.add(k);
    this.queue.push({ keys, job });
    while (this.queue.length > SIG_LIMITS.queue) for (const k of this.queue.shift().keys) this.pending.delete(k);
    this.pump();
  }
  pump() {
    while (this.active < SIG_LIMITS.parallel && this.queue.length) {
      const { keys, job } = this.queue.pop(); // newest first: what the view shows now
      this.active++;
      job().then(() => { this.active--; for (const k of keys) this.pending.delete(k); this.pump(); if (this.file.onData) this.file.onData(); },
        e => { this.active--; this.file.failed = e; this.queue = []; this.pending.clear(); if (this.file.onData) this.file.onData(e); });
    }
  }
  get busy() { return this.active > 0 || this.queue.length > 0; }
}
// Reads blocks that are not in the cache yet: blocks next to each other in the file come in one read
// (up to 2 MB), so a view needs few reads, also over a network. items: { key, pos, size, decode(u8) }
function sigFetch(file, items) {
  items = items.filter(it => !file.loader.pending.has(it.key)).sort((x, y) => x.pos - y.pos);
  let run = [];
  const flush = () => {
    if (!run.length) return;
    const r = run; run = [];
    const p0 = r[0].pos, p1 = r[r.length - 1].pos + r[r.length - 1].size;
    file.loader.want(r.map(x => x.key), async () => {
      const buf = await sigRead(file.src, p0, p1 - p0, Math.max(SIG_LIMITS.block, 2 << 20), file.code);
      for (const it of r) { const blk = await it.decode(buf.subarray(it.pos - p0, it.pos - p0 + it.size)); sigCachePut(it.key, blk, blockBytes(blk)); }
    });
  };
  for (const it of items) {
    if (!(it.size > 0) || it.size > SIG_LIMITS.block) { file.failed = new CgvError(file.code, { what: 'a block is too large' }); return; }
    const last = run[run.length - 1];
    if (last && (it.pos < last.pos + last.size || it.pos - (last.pos + last.size) > 65536 || it.pos + it.size - run[0].pos > (2 << 20))) flush();
    run.push(it);
  }
  flush();
}

const sigBudget = f => (f.src.url ? SIG_LIMITS.viewBytesUrl : SIG_LIMITS.viewBytes); // bytes of full-resolution data per view
async function sigRead(src, pos, size, max, code, what) {
  if (!(pos >= 0) || !(size >= 0) || size > max) throw new CgvError(code, { what: what || 'a block is too large' });
  if (size === 0) return new Uint8Array(0);
  const b = await src.read(pos, size);
  if (b.length < size) throw new CgvError(code, { what: 'cut short' });
  return b.length > size ? b.subarray(0, size) : b;
}

// ======================= TDF =======================
// Header (little-endian): magic "TDF" (or "IBF"), version, index position (8), index size, header
// size; then (version 2+) window functions; track type; track line; track names; (version 3+) genome
// and flags (bit 0: tiles compressed with zlib). The master index lists datasets
// ("/chr1/raw", "/chr1/z3/mean"; version 1 without the function) and groups ("/" holds chromosomes,
// maxZoom and more). A dataset lists tiles; a tile is fixedStep, variableStep, bed or bedWithName
// with one row of values per track. Zoom z of a contig has 2^z tiles of 700 bins.
class TdfFile {
  constructor(src) { this.src = src; this.name = src.name; this.format = 'tdf'; this.code = 'tdfBad'; this.ds = new Map(); this.loader = new SigLoader(this); this.id = ++TdfFile.seq; }
  bad(what) { return new CgvError('tdfBad', { what }); }
  async open() {
    const h = await this.src.read(0, 24);
    if (h.length < 24) throw this.bad('cut short');
    const magic = String.fromCharCode(h[0], h[1], h[2]);
    if (magic !== 'TDF' && magic !== 'IBF') throw this.bad('no TDF signature');
    const r = new SigBin(h); r.p = 4;
    this.version = r.i32();
    const indexPos = r.i64(), indexSize = r.i32(), headerSize = r.i32();
    if (this.version < 1 || this.version > 100) throw this.bad('version ' + this.version);
    const hb = new SigBin(await sigRead(this.src, 24, headerSize, SIG_LIMITS.header, 'tdfBad', 'the header is too large'));
    if (this.version >= 2) { const n = hb.count(1, 1000, 'window functions'); this.windowFunctions = []; for (let i = 0; i < n; i++) this.windowFunctions.push(hb.str()); }
    else this.windowFunctions = ['mean'];
    this.trackType = hb.str(); this.trackLine = hb.str().trim();
    const nt = hb.count(1, SIG_LIMITS.tracks, 'track names'); this.trackNames = [];
    for (let i = 0; i < nt; i++) this.trackNames.push(hb.str());
    if (!nt) throw this.bad('no tracks');
    if (this.version > 2) { this.genomeId = hb.str(); this.compressed = (hb.i32() & 1) !== 0; } else this.compressed = false;
    const ix = new SigBin(await sigRead(this.src, indexPos, indexSize, SIG_LIMITS.index, 'tdfBad', 'the index is too large'));
    this.dsIndex = new Map(); this.grIndex = new Map();
    for (const map of [this.dsIndex, this.grIndex]) {
      const n = ix.count(13, SIG_LIMITS.items, 'index');
      for (let i = 0; i < n; i++) { const k = ix.str(); map.set(k, { pos: ix.i64(), size: ix.i32() }); }
    }
    this.root = (await this.group('/')) || {};
    this.maxZoom = /^-?\d+$/.test(this.root.maxZoom || '') ? +this.root.maxZoom : -1;
    // contigs: from the root group and from dataset names, matched as LGV matches names (chr7 = 7)
    this.chrs = new Map();
    const addChr = c => { if (c && c.toLowerCase() !== 'all' && !this.chrs.has(contigKey(c))) this.chrs.set(contigKey(c), c); };
    for (const c of String(this.root.chromosomes || '').split(',')) addChr(c);
    for (const k of this.dsIndex.keys()) addChr(k.split('/')[1]);
    this.props = parseTrackLine(this.trackLine);
    this.tracks = this.trackNames.slice();
    return this;
  }
  async group(name) {
    const e = this.grIndex.get(name); if (!e) return null;
    const r = new SigBin(await sigRead(this.src, e.pos, e.size, SIG_LIMITS.header, 'tdfBad'));
    const g = Object.create(null), n = r.count(2, 100000, 'group');
    for (let i = 0; i < n; i++) { const k = r.str(); g[k] = r.str(); }
    return g;
  }
  // dataset: attributes, data type, tile width, tile positions and sizes; null when there is none
  dataset(name) {
    if (this.ds.has(name)) return this.ds.get(name);
    const e = this.dsIndex.get(name);
    if (!e) { this.ds.set(name, null); return null; }
    this.loader.want('ds:' + name, async () => {
      const r = new SigBin(await sigRead(this.src, e.pos, e.size, SIG_LIMITS.index, 'tdfBad'));
      const na = r.count(2, 100000, 'dataset'); const attrs = Object.create(null);
      for (let i = 0; i < na; i++) { const k = r.str(); attrs[k] = r.str(); }
      r.str(); const tileWidth = r.f32(), n = r.count(12, SIG_LIMITS.items, 'tiles');
      if (!(tileWidth >= 1 && tileWidth < 1e10)) throw this.bad('tile width');
      const pos = new Float64Array(n), size = new Int32Array(n);
      for (let i = 0; i < n; i++) { pos[i] = r.i64(); size[i] = r.i32(); }
      this.ds.set(name, { name, attrs, tileWidth, pos, size });
    });
    return undefined; // being read
  }
  async readTile(d, i) { return this.decodeTile(await sigRead(this.src, d.pos[i], d.size[i], SIG_LIMITS.block, 'tdfBad')); }
  async decodeTile(u8) {
    if (this.compressed) u8 = await inflateStream(u8, 'deflate', SIG_LIMITS.block, 'tdfBad', { what: 'a tile does not decompress' });
    return decodeTdfTile(u8, this.tracks.length);
  }
  // Which dataset a view uses. zoom: what igvtools precomputed for a 700 px wide view of the contig.
  // For values one by one (points, or a function the file did not precompute) the raw dataset is used
  // if the view needs no more than SIG_LIMITS.viewBytes of it; otherwise the precomputed summaries.
  plan(contig, a, b, bpp, chrLen, want) {
    const c = this.chrs.get(contigKey(contig)); if (!c) return { none: true };
    let z = Math.ceil(Math.log2(Math.max(1, chrLen / (bpp * 700)))); if (!(z >= 0)) z = 0;
    const wfs = this.windowFunctions, has = f => wfs.includes(f);
    const zname = (f, zz) => `/${c}/z${zz}` + (this.version >= 2 ? '/' + f : '');
    // summaries of zoom level zz, or of the finest coarser level the file has (a file may stop at a
    // low maxZoom: its coarsest summaries are then still better than nothing)
    const summary = () => {
      for (let zz = Math.min(z, Math.max(0, this.maxZoom)); zz >= 0; zz--) {
        const fields = {};
        for (const f of ['mean', 'min', 'max']) if ((this.version >= 2 ? has(f) : f === 'mean') && this.dsIndex.has(zname(f, zz))) fields[f] = zname(f, zz);
        if (!Object.keys(fields).length && this.version >= 2) for (const f of wfs) if (this.dsIndex.has(zname(f, zz))) { fields.mean = zname(f, zz); break; }
        if (Object.keys(fields).length) return { mode: 'summary', fields, z: zz };
      }
      return null;
    };
    const raw = `/${c}/raw`;
    if (z > this.maxZoom || want === 'raw') {
      if (this.dsIndex.has(raw)) {
        const d = this.dataset(raw); if (d === undefined) return { wait: true };
        if (d) {
          const t0 = Math.max(0, Math.floor(a / d.tileWidth)), t1 = Math.min(d.pos.length - 1, Math.floor(b / d.tileWidth));
          let bytes = 0; for (let i = t0; i <= t1; i++) if (d.pos[i] >= 0) bytes += d.size[i];
          // neighbouring tiles come in one read, so the limit is on bytes; the tile count only bounds the work
          if (bytes <= sigBudget(this) && t1 - t0 < SIG_LIMITS.rawTiles) return { mode: 'raw', fields: { v: raw }, z };
        }
      }
    }
    return summary() || (this.dsIndex.has(raw) ? { tooMuch: true } : { none: true });
  }
  // Blocks over [a, b) for a plan; tiles not read yet are asked for and arrive through onData.
  blocks(plan, a, b) {
    const out = []; let missing = 0;
    const fields = Object.entries(plan.fields), ds = fields.map(([, n]) => this.dataset(n));
    if (ds.some(d => d === undefined)) return { blocks: out, missing: 1 };
    const d0 = ds.find(Boolean); if (!d0) return { blocks: out, missing: 0 };
    const t0 = Math.max(0, Math.floor(a / d0.tileWidth)), t1 = Math.min(d0.pos.length - 1, Math.floor(b / d0.tileWidth));
    if (t1 - t0 >= SIG_LIMITS.viewBlocks) return { blocks: out, missing: 0, tooMuch: true };
    const need = [];
    for (let i = t0; i <= t1; i++) {
      const parts = {}; let ready = true;
      fields.forEach(([f, n], k) => {
        const d = ds[k]; if (!d || i >= d.pos.length || d.pos[i] < 0 || d.size[i] <= 0) return;
        const key = `${this.id}|${n}|${i}`, got = sigCacheGet(key);
        if (got) parts[f] = got;
        else { ready = false; need.push({ key, pos: d.pos[i], size: d.size[i], decode: u8 => this.decodeTile(u8) }); }
      });
      if (!ready) { missing++; continue; }
      const ks = Object.keys(parts); if (!ks.length) continue;
      if (plan.mode === 'raw') out.push(parts.v);
      else { // several summary tiles of one zoom share positions
        const base = parts[ks[0]], s = {};
        for (const f of ks) if (parts[f].start.length === base.start.length) s[f] = parts[f].v;
        out.push({ start: base.start, end: base.end, s });
      }
    }
    if (need.length) sigFetch(this, need);
    return { blocks: out, missing };
  }
}
TdfFile.seq = 0;

function decodeTdfTile(u8, nTracks) {
  const r = new SigBin(u8, true, 'tdfBad'), type = r.str();
  const rows = (n, nS) => { const v = []; for (let t = 0; t < nS; t++) { r.need(n * 4); const a = new Float32Array(n); for (let k = 0; k < n; k++) a[k] = r.f32(); if (t < nTracks) v.push(a); } while (v.length < nTracks) v.push(new Float32Array(n).fill(NaN)); return v; };
  if (type === 'fixedStep') {
    const n = r.count(4, SIG_LIMITS.items, 'tile'), start = r.i32(), span = r.f32();
    if (!(span > 0 && span < 1e10)) throw r.bad('tile step');
    if (n * 4 * nTracks > r.left) throw r.bad('tile cut short');
    const st = new Int32Array(n), en = new Int32Array(n);
    for (let k = 0; k < n; k++) { st[k] = start + Math.floor(k * span); en[k] = start + Math.floor((k + 1) * span); }
    return { start: st, end: en, v: rows(n, nTracks) };
  }
  if (type === 'variableStep') {
    r.i32(); const span = r.f32(), n = r.count(4, SIG_LIMITS.items, 'tile');
    if (!(span > 0 && span < 1e10)) throw r.bad('tile span');
    const st = new Int32Array(n), en = new Int32Array(n);
    for (let k = 0; k < n; k++) { st[k] = r.i32(); en[k] = st[k] + Math.max(1, Math.floor(span)); }
    const nS = r.count(0, SIG_LIMITS.tracks, 'tracks in a tile');
    if (n * 4 * nS > r.left) throw r.bad('tile cut short');
    return { start: st, end: en, v: rows(n, nS) };
  }
  if (type === 'bed' || type === 'bedWithName') {
    const n = r.count(8, SIG_LIMITS.items, 'tile'), st = new Int32Array(n), en = new Int32Array(n);
    for (let k = 0; k < n; k++) st[k] = r.i32();
    for (let k = 0; k < n; k++) { en[k] = r.i32(); if (en[k] <= st[k]) en[k] = st[k] + 1; }
    const nS = r.count(0, SIG_LIMITS.tracks, 'tracks in a tile');
    if (n * 4 * nS > r.left) throw r.bad('tile cut short');
    const out = { start: st, end: en, v: rows(n, nS) };
    if (type === 'bedWithName') { out.names = []; for (let k = 0; k < n; k++) out.names.push(r.str().slice(0, 200)); }
    return out;
  }
  throw r.bad('unknown tile type');
}

// ======================= bigWig =======================
const BW_MAGIC = 0x888FFC26, BB_MAGIC = 0x8789F2EB, BPT_MAGIC = 0x78CA8C91, CIR_MAGIC = 0x2468ACE0;
class BigWigFile {
  constructor(src) { this.src = src; this.name = src.name; this.format = 'bigwig'; this.code = 'bwBad'; this.loader = new SigLoader(this); this.nodes = new Map(); this.leafCache = new Map(); this.id = ++TdfFile.seq; }
  bad(what) { return new CgvError('bwBad', { what }); }
  bin(u8) { return new SigBin(u8, this.le, 'bwBad'); }
  async open() {
    const h = await this.src.read(0, 64);
    if (h.length < 64) throw this.bad('cut short');
    const dv = new DataView(h.buffer, h.byteOffset, 64);
    if (dv.getUint32(0, true) === BW_MAGIC) this.le = true;
    else if (dv.getUint32(0, false) === BW_MAGIC) this.le = false;
    else if (dv.getUint32(0, true) === BB_MAGIC || dv.getUint32(0, false) === BB_MAGIC) throw new CgvError('bwBigBed', { name: this.name });
    else throw this.bad('no bigWig signature');
    const r = this.bin(h); r.p = 4;
    this.version = r.u16(); const nZoom = r.u16();
    this.chromTreeOffset = r.i64(); this.fullDataOffset = r.i64(); this.fullIndexOffset = r.i64();
    r.u16(); r.u16(); r.i64(); this.totalSummaryOffset = r.i64(); this.uncompressBufSize = r.u32();
    if (nZoom > 64) throw this.bad('zoom levels');
    if (this.uncompressBufSize > SIG_LIMITS.block) throw this.bad('the block size is too large');
    const zr = this.bin(await sigRead(this.src, 64, nZoom * 24, SIG_LIMITS.header, 'bwBad'));
    this.zooms = [];
    for (let i = 0; i < nZoom; i++) { const reduction = zr.u32(); zr.u32(); const dataOffset = zr.i64(), indexOffset = zr.i64(); if (reduction > 0 && indexOffset > 0) this.zooms.push({ reduction, indexOffset }); }
    this.zooms.sort((x, y) => x.reduction - y.reduction);
    if (this.totalSummaryOffset > 0) {
      const s = this.bin(await sigRead(this.src, this.totalSummaryOffset, 40, 40, 'bwBad'));
      const n = s.i64(), mn = s.f64(), mx = s.f64(), sum = s.f64(), sq = s.f64();
      if (n > 0 && isFinite(mn) && isFinite(mx)) this.total = { n, min: mn, max: mx, mean: sum / n, sd: Math.sqrt(Math.max(0, sq / n - (sum / n) ** 2)) };
    }
    await this.readChromTree();
    this.trackLine = ''; this.props = {}; this.trackType = '';
    this.tracks = [this.name.replace(/\.(bw|bigwig)$/i, '')];
    return this;
  }
  // chromosome B+ tree: header (magic, block size, key size, value size, item count), then nodes
  // (leaf: key, chromId, chromSize; inner: key, child offset), read whole
  async readChromTree() {
    const hr = this.bin(await sigRead(this.src, this.chromTreeOffset, 32, 32, 'bwBad'));
    if (hr.u32() !== BPT_MAGIC) throw this.bad('chromosome tree');
    const blockSize = hr.u32(), keySize = hr.u32(), valSize = hr.u32(), itemCount = hr.i64();
    if (!(keySize >= 1 && keySize <= 255) || valSize !== 8 || !(blockSize >= 1 && blockSize <= 65536) || itemCount > SIG_LIMITS.chroms) throw this.bad('chromosome tree');
    this.chrById = new Map(); this.chrs = new Map(); let nodes = 0;
    const walk = async (off, depth) => {
      if (depth > SIG_LIMITS.depth || ++nodes > SIG_LIMITS.nodes) throw this.bad('chromosome tree');
      const head = this.bin(await sigRead(this.src, off, 4, 4, 'bwBad'));
      const leaf = head.u8v(); head.u8v(); const n = head.u16();
      if (n > blockSize) throw this.bad('chromosome tree');
      const r = this.bin(await sigRead(this.src, off + 4, n * (keySize + 8), SIG_LIMITS.header, 'bwBad'));
      const kids = [];
      for (let i = 0; i < n; i++) {
        const key = r.fixedStr(keySize);
        if (leaf) { const id = r.u32(), size = r.u32(); if (this.chrById.size >= SIG_LIMITS.chroms) throw this.bad('too many contigs'); this.chrById.set(id, key); if (!this.chrs.has(contigKey(key))) this.chrs.set(contigKey(key), { id, size, name: key }); }
        else kids.push(r.i64());
      }
      for (const k of kids) await walk(k, depth + 1);
    };
    await walk(this.chromTreeOffset + 32, 0);
  }
  // R tree (cirTree) node at an offset: leaves give data blocks (offset, size); read on demand, kept
  async node(off) {
    const k = off; if (this.nodes.has(k)) return this.nodes.get(k);
    if (this.nodes.size > SIG_LIMITS.nodes) throw this.bad('index tree');
    const head = this.bin(await sigRead(this.src, off, 4, 4, 'bwBad'));
    const leaf = head.u8v(); head.u8v(); const n = head.u16();
    const r = this.bin(await sigRead(this.src, off + 4, n * (leaf ? 32 : 24), SIG_LIMITS.header, 'bwBad'));
    const items = [];
    for (let i = 0; i < n; i++) {
      const it = { c0: r.u32(), s0: r.u32(), c1: r.u32(), s1: r.u32(), off: r.i64() };
      if (leaf) it.size = r.i64();
      items.push(it);
    }
    const nd = { leaf: !!leaf, items }; this.nodes.set(k, nd); return nd;
  }
  // data blocks of an index (the full one or a zoom level's) overlapping [a, b) of contig id
  async leaves(indexOffset, id, a, b) {
    const key = `${indexOffset}|${id}|${a}|${b}`; if (this.leafCache.has(key)) return this.leafCache.get(key);
    const hr = this.bin(await sigRead(this.src, indexOffset, 48, 48, 'bwBad'));
    if (hr.u32() !== CIR_MAGIC) throw this.bad('index tree');
    const out = [];
    const over = it => (it.c0 < id || (it.c0 === id && it.s0 < b)) && (it.c1 > id || (it.c1 === id && it.s1 > a));
    const walk = async (off, depth) => {
      if (depth > SIG_LIMITS.depth) throw this.bad('index tree');
      const nd = await this.node(off);
      for (const it of nd.items) if (over(it)) { if (nd.leaf) { if (out.length < SIG_LIMITS.viewBlocks + 1) out.push({ off: it.off, size: it.size }); } else await walk(it.off, depth + 1); }
    };
    await walk(indexOffset + 48, 0);
    if (this.leafCache.size > 2000) this.leafCache.clear();
    this.leafCache.set(key, out); return out;
  }
  async readBlock(lf, zoom, id) { return this.decodeBlock(await sigRead(this.src, lf.off, lf.size, SIG_LIMITS.block, 'bwBad'), zoom, id); }
  async decodeBlock(u8, zoom, id) {
    if (this.uncompressBufSize > 0) u8 = await inflateStream(u8, 'deflate', Math.max(this.uncompressBufSize, 1 << 20), 'bwBad', { what: 'a block does not decompress' });
    return zoom ? decodeBwZoom(u8, this.le, id) : decodeBwData(u8, this.le, id);
  }
  plan(contig, a, b, bpp, chrLen, want) {
    const c = this.chrs.get(contigKey(contig)); if (!c) return { none: true };
    const lk = `L|${c.id}|${Math.floor(a)}|${Math.ceil(b)}`;
    // full-resolution data when the view needs little of it (or the file has no summaries)
    const zl = this.zooms.filter(z => z.reduction <= bpp).pop();
    const tryRaw = want === 'raw' || !zl;
    if (tryRaw) {
      const L = this.leafList(this.fullIndexOffset, c.id, a, b, lk + '|full'); if (!L) return { wait: true };
      let bytes = 0; for (const x of L) bytes += x.size;
      if (bytes <= sigBudget(this) && L.length <= SIG_LIMITS.viewBlocks) return { mode: 'raw', id: c.id, leaves: L };
      if (!zl && !this.zooms.length) return { tooMuch: true };
    }
    // summaries: the coarsest level that still gives a value per pixel (the finest one when zoomed in)
    const z = zl || this.zooms[0]; if (!z) return { tooMuch: true };
    const L = this.leafList(z.indexOffset, c.id, a, b, lk + '|' + z.reduction); if (!L) return { wait: true };
    if (L.length > SIG_LIMITS.viewBlocks) return { tooMuch: true };
    return { mode: 'summary', id: c.id, leaves: L, zoom: z.reduction };
  }
  leafList(indexOffset, id, a, b, key) {
    const hit = this.leafCache.get(key); if (hit) return hit;
    this.loader.want(key, async () => { const L = await this.leaves(indexOffset, id, Math.floor(a), Math.ceil(b)); this.leafCache.set(key, L); });
    return null;
  }
  blocks(plan) {
    const out = [], need = [];
    for (const lf of plan.leaves) {
      const key = `${this.id}|${lf.off}|${plan.id}`, got = sigCacheGet(key); // a block can hold several contigs: decoded per contig
      if (got) { out.push(got); continue; }
      need.push({ key, pos: lf.off, size: lf.size, decode: u8 => this.decodeBlock(u8, plan.mode === 'summary', plan.id) });
    }
    if (need.length) sigFetch(this, need);
    return { blocks: out, missing: need.length };
  }
}

// ======================= bigBed =======================
// Features (BED records) in the bbi container (UCSC bbiFile.h): same header, contig tree and R tree as
// bigWig; data blocks hold records chromId, start, end, then the remaining BED columns as text. Read
// by region; optional extra indexes (B+ trees on a column, e.g. the name) find a feature by name.
const BB_LIMITS = { viewBytes: 8 << 20, sql: 1 << 16, names: 64, blocks: 2000 };
class BigBedFile extends BigWigFile {
  constructor(src) { super(src); this.format = 'bigbed'; this.blockCache = new Map(); }
  bad(what) { return new CgvError('bbBad', { what }); }
  bin(u8) { return new SigBin(u8, this.le, 'bbBad'); }
  async open() {
    const h = await this.src.read(0, 64);
    if (h.length < 64) throw this.bad('cut short');
    const dv = new DataView(h.buffer, h.byteOffset, 64);
    if (dv.getUint32(0, true) === BB_MAGIC) this.le = true;
    else if (dv.getUint32(0, false) === BB_MAGIC) this.le = false;
    else throw this.bad('no bigBed signature');
    const r = this.bin(h); r.p = 4;
    this.version = r.u16(); r.u16();
    this.chromTreeOffset = r.i64(); this.fullDataOffset = r.i64(); this.fullIndexOffset = r.i64();
    this.fieldCount = r.u16(); this.definedFieldCount = r.u16();
    const sqlOff = r.i64(); r.i64(); this.uncompressBufSize = r.u32(); const extOff = this.version >= 4 ? r.i64() : 0;
    if (this.uncompressBufSize > SIG_LIMITS.block) throw this.bad('the block size is too large');
    // column names from the autoSql text (fields after the BED ones are shown in the card)
    this.fields = [];
    if (sqlOff > 0) {
      const t = new TextDecoder().decode(await this.src.read(sqlOff, BB_LIMITS.sql)).split('\0')[0];
      // "table name "comment" ( type field; "comment" ... )": the last word before each ';'
      const body = t.slice(t.indexOf('(') + 1).replace(/"[^"]*"/g, ' ');
      for (const part of body.split(';').slice(0, -1)) { const w = part.trim().split(/\s+/); this.fields.push(w[w.length - 1].replace(/\W/g, '')); }
    }
    // extra indexes: [{ field, offset }]
    this.extra = [];
    if (extOff > 0) {
      const e = this.bin(await sigRead(this.src, extOff, 12, 12, 'bbBad')); e.u16(); const n = e.u16(), listOff = e.i64();
      if (n > BB_LIMITS.names) throw this.bad('extra indexes');
      if (n && listOff > 0) {
        const l = this.bin(await sigRead(this.src, listOff, n * 20, n * 20 + 64, 'bbBad'));
        for (let i = 0; i < n; i++) { l.p = i * 20; l.u16(); const fc = l.u16(), off = l.i64(); l.u32(); const field = l.u16(); if (fc >= 1 && off > 0) this.extra.push({ field, offset: off }); }
      }
    }
    await this.readChromTree();
    this.tracks = [this.name.replace(/\.(bb|bigbed)$/i, '')];
    return this;
  }
  // records of one block, decoded once and kept (a few hundred blocks)
  async block(lf) {
    const k = lf.off; let b = this.blockCache.get(k); if (b) return b;
    let u8 = await sigRead(this.src, lf.off, lf.size, SIG_LIMITS.block, 'bbBad');
    if (this.uncompressBufSize > 0) u8 = await inflateStream(u8, 'deflate', Math.max(this.uncompressBufSize, 1 << 20), 'bbBad', { what: 'a block does not decompress' });
    const dv = new DataView(u8.buffer, u8.byteOffset, u8.byteLength), td = new TextDecoder(), out = [];
    for (let p = 0; p + 12 <= u8.length;) {
      const cid = dv.getUint32(p, this.le), s = dv.getUint32(p + 4, this.le), e = dv.getUint32(p + 8, this.le); p += 12;
      let z = p; while (z < u8.length && u8[z] !== 0) z++;
      if (z >= u8.length) throw this.bad('a record is cut short');
      out.push({ cid, s, e, rest: td.decode(u8.subarray(p, z)) }); p = z + 1;
    }
    if (this.blockCache.size > BB_LIMITS.blocks) this.blockCache.clear();
    this.blockCache.set(k, out); return out;
  }
  feature(rec) {
    const chr = this.chrById.get(rec.cid) || '';
    const f = [chr, String(rec.s), String(rec.e), ...(rec.rest ? rec.rest.split('\t') : [])];
    const ft = bedFeature(f); if (!ft) return null;
    // columns past BED 12 (or past the defined ones): named by the autoSql text
    const from = Math.min(12, this.definedFieldCount || 12);
    if (f.length > from) { ft.attrs = {}; for (let i = from; i < f.length; i++) if (f[i] !== '') ft.attrs[this.fields[i] || 'field' + (i + 1)] = f[i]; }
    // gene tracks (bigGenePred and the like): the gene symbol as the label, the transcript in the card
    const g = ft.attrs && (ft.attrs.geneName || ft.attrs.name2 || ft.attrs.geneSymbol || ft.attrs.gene_name);
    if (g) ft.label = g;
    return ft;
  }
  // features overlapping [a, b) of a contig; { tooMuch } when the view would read more than allowed
  async features(contig, a, b, noLimit) {
    const c = this.chrs.get(contigKey(contig)); if (!c) return [];
    const L = await this.leaves(this.fullIndexOffset, c.id, Math.floor(a), Math.ceil(b));
    let bytes = 0; for (const x of L) bytes += x.size;
    if (!noLimit && (bytes > BB_LIMITS.viewBytes || L.length > SIG_LIMITS.viewBlocks)) return { tooMuch: true };
    const out = [];
    for (const lf of L) for (const r of await this.block(lf)) if (r.cid === c.id && r.s < b && r.e > a) { const ft = this.feature(r); if (ft) out.push(ft); }
    return out;
  }
  // features whose indexed column equals a name (any case), through the extra B+ tree indexes
  async findName(q) {
    const want = String(q).trim(); if (!want || !this.extra.length) return [];
    const out = [];
    for (const ex of this.extra) {
      const bare = want.replace(/\.\d+$/, ''), blocks = await this.bptFind(ex.offset, [want, want.toUpperCase(), bare, bare.toUpperCase()]);
      for (const lf of blocks) for (const r of await this.block(lf)) {
        const f = [String(r.cid), String(r.s), String(r.e), ...r.rest.split('\t')];
        const v = f[ex.field]; if (v === undefined) continue;
        if (v.toUpperCase() === want.toUpperCase() || v.replace(/\.\d+$/, '').toUpperCase() === want.toUpperCase()) { const ft = this.feature(r); if (ft) out.push(ft); }
      }
      if (out.length) break;
    }
    return out;
  }
  // B+ tree (bPlusTree.h) lookup: values (block offset, size) of every item whose key is one of keys
  async bptFind(off, keys) {
    const hr = this.bin(await sigRead(this.src, off, 32, 32, 'bbBad'));
    if (hr.u32() !== BPT_MAGIC) throw this.bad('name index');
    const blockSize = hr.u32(), keySize = hr.u32(), valSize = hr.u32();
    if (!(keySize >= 1 && keySize <= 1024) || valSize !== 16 || !(blockSize >= 2 && blockSize <= 65536)) throw this.bad('name index');
    const ks = [...new Set(keys.filter(k => k.length <= keySize))]; if (!ks.length) return [];
    const found = new Map(); let nodes = 0;
    const walk = async (p, depth) => {
      if (depth > SIG_LIMITS.depth || ++nodes > 200) throw this.bad('name index');
      const head = this.bin(await sigRead(this.src, p, 4, 4, 'bbBad')); const leaf = head.u8v(); head.u8v(); const n = head.u16();
      if (n > blockSize) throw this.bad('name index');
      const r = this.bin(await sigRead(this.src, p + 4, n * (keySize + (leaf ? 16 : 8)), SIG_LIMITS.header, 'bbBad'));
      const items = [];
      for (let i = 0; i < n; i++) { const key = r.fixedStr(keySize); if (leaf) { const o = r.i64(), sz = r.i64(); items.push({ key, o, sz }); } else items.push({ key, o: r.i64() }); }
      if (leaf) { for (const it of items) if (ks.includes(it.key)) found.set(it.o, { off: it.o, size: it.sz }); return; }
      // a child may hold a key when its first key is not past it and the next child's is not before it
      for (let i = 0; i < n; i++) if (ks.some(k => items[i].key <= k && (i + 1 === n || items[i + 1].key >= k))) await walk(items[i].o, depth + 1);
    };
    await walk(off + 32, 0);
    return [...found.values()];
  }
}
// A bigBed as an annotation track's index: query() answers from the window read last and asks for a
// new one when the view leaves it (onData redraws when it arrives).
class BigBedIndex {
  constructor(file) { this.file = file; this.win = null; this.busy = null; this.onData = null; this.byChr = null; this.count = 0; this.skipped = 0; }
  query(contig, a, b) {
    const w = this.win, k = contigKey(contig);
    if (!w || w.k !== k || a < w.a || b > w.b) this.load(contig, a, b);
    if (!w || w.k !== k || w.tooMuch) { this.tooMuch = !!(w && w.k === k && w.tooMuch && a >= w.a && b <= w.b); return []; }
    this.tooMuch = false;
    return w.list.filter(f => f.start < b && f.end > a);
  }
  load(contig, a, b) {
    const k = contigKey(contig), span = b - a, A = Math.max(0, Math.floor(a - span)), B = Math.ceil(b + span);
    if (this.busy && this.busy.k === k && a >= this.busy.a && b <= this.busy.b) return;
    const job = this.busy = { k, a: A, b: B };
    this.file.features(contig, A, B).then(r => {
      if (this.busy !== job) return;
      this.busy = null;
      this.win = r && r.tooMuch ? { k, a: A, b: B, tooMuch: true, list: [] } : { k, a: A, b: B, list: r };
      if (this.onData) this.onData();
    }, e => { if (this.busy === job) this.busy = null; if (this.onData) this.onData(e); });
  }
}
async function openBigBed(src) { return new BigBedFile(src).open(); }

// bigWig data section: chromId, start, end, step, span, type (1 bedGraph, 2 variableStep, 3 fixedStep),
// reserved, item count; then the items. A block may hold several sections; only contig id is kept.
function decodeBwData(u8, le, id) {
  const r = new SigBin(u8, le, 'bwBad'), S = [], E = [], V = [];
  while (r.left >= 24) {
    const cid = r.u32(), start = r.u32(); r.u32(); const step = r.u32(), span = r.u32(), type = r.u8v(); r.u8v(); const n = r.u16();
    const each = type === 1 ? 12 : type === 2 ? 8 : type === 3 ? 4 : 0;
    if (!each) throw r.bad('unknown section type');
    if (n * each > r.left) throw r.bad('section cut short');
    for (let i = 0; i < n; i++) {
      let s, e, v;
      if (type === 1) { s = r.u32(); e = r.u32(); v = r.f32(); }
      else if (type === 2) { s = r.u32(); v = r.f32(); e = s + span; }
      else { v = r.f32(); s = start + i * step; e = s + span; }
      if (cid !== id || !Number.isFinite(v) || e <= s) continue;
      S.push(s); E.push(e); V.push(v);
    }
  }
  return { start: Int32Array.from(S), end: Int32Array.from(E), v: [Float32Array.from(V)] };
}
// bigWig zoom records: chromId, start, end, validCount, min, max, sum, sum of squares (per base)
function decodeBwZoom(u8, le, id) {
  const r = new SigBin(u8, le, 'bwBad'), n = Math.floor(u8.length / 32);
  const S = new Int32Array(n), E = new Int32Array(n), mean = new Float32Array(n), mn = new Float32Array(n), mx = new Float32Array(n), sd = new Float32Array(n), cnt = new Float32Array(n);
  let k = 0;
  for (let i = 0; i < n; i++) {
    const cid = r.u32(), s = r.u32(), e = r.u32(), c = r.u32(), a = r.f32(), b = r.f32(), sum = r.f32(), sq = r.f32();
    if (cid !== id || !c || e <= s) continue;
    S[k] = s; E[k] = e; mn[k] = a; mx[k] = b; mean[k] = sum / c; cnt[k] = c;
    sd[k] = c > 1 ? Math.sqrt(Math.max(0, (sq - sum * sum / c) / (c - 1))) : 0; // as UCSC's calcStdFromSums
    k++;
  }
  const cut = a => a.subarray(0, k);
  return { start: cut(S), end: cut(E), s: { mean: [cut(mean)], min: [cut(mn)], max: [cut(mx)], sd: [cut(sd)], n: [cut(cnt)] } };
}

// Opens a signal file by name: .tdf or .bw / .bigwig.
function signalKind(name) { const n = String(name).toLowerCase(); return /\.tdf$/.test(n) ? 'tdf' : /\.(bw|bigwig)$/.test(n) ? 'bigwig' : null; }
const isBigBed = name => /\.(bb|bigbed)$/i.test(String(name));
async function openSignal(src) {
  const k = signalKind(src.name);
  const f = k === 'tdf' ? new TdfFile(src) : new BigWigFile(src);
  await f.open();
  f.kind = 'signal';
  return f;
}

// ---- the default way to draw a file that does not say ----
// A track line's graphType wins. A TDF of type ALLELE_FREQUENCY is drawn as points. Otherwise a TDF
// whose values all lie in 0–1 and cover less than half of the stretch they span (separate positions,
// as variant allele frequencies are) is drawn as points on a 0–1 scale; everything else as bars, the
// UCSC default. The sample is the first non-empty tile of the first raw dataset: one small read.
async function signalDefaults(f) {
  const p = f.props || {}, d = { graph: p.graphType || 'bar', why: p.graphType ? 'trackLine' : 'ucsc' };
  if (p.graphType) return d;
  if (f.format === 'tdf' && /^ALLELE_FREQUENCY$/i.test(f.trackType)) return { graph: 'points', why: 'type', fixed01: true };
  if (f.format !== 'tdf') return d;
  try {
    const name = [...f.dsIndex.keys()].find(k => /\/raw$/.test(k)); if (!name) return d;
    const e = f.dsIndex.get(name);
    const r = new SigBin(await sigRead(f.src, e.pos, e.size, SIG_LIMITS.index, 'tdfBad'));
    const na = r.count(2, 100000); for (let i = 0; i < na; i++) { r.str(); r.str(); }
    r.str(); r.f32(); const n = r.count(12, SIG_LIMITS.items);
    for (let i = 0; i < n; i++) {
      const pos = r.i64(), size = r.i32(); if (pos < 0 || size <= 0) continue;
      const t = await f.readTile({ pos: [pos], size: [size] }, 0);
      const m = t.start.length; if (m < 2) continue;
      let cov = 0, lo = Infinity, hi = -Infinity, ok = 0;
      for (let k = 0; k < m; k++) { cov += t.end[k] - t.start[k]; for (const v of t.v) { const x = v[k]; if (Number.isNaN(x)) continue; ok++; if (x < lo) lo = x; if (x > hi) hi = x; } }
      const spanBp = t.end[m - 1] - t.start[0];
      if (ok && lo >= 0 && hi <= 1 && spanBp > 0 && cov / spanBp < 0.5) return { graph: 'points', why: 'rule', fixed01: true };
      return d;
    }
  } catch (e) { /* the tracks still open; a damaged tile shows when it is drawn */ }
  return d;
}

if (typeof module !== 'undefined') module.exports = { SIG_LIMITS, SigBin, parseTrackLine, TdfFile, BigWigFile, BigBedFile, BigBedIndex, openBigBed, isBigBed, decodeTdfTile, decodeBwData, decodeBwZoom, signalKind, openSignal, signalDefaults, SIG_CACHE };
