// CGV: user annotation tracks — BED, bedGraph, VCF, GTF and GFF3 files (plain or gzip), read whole
// into memory and indexed by contig. Own code, no dependencies.
//
// Safety rules (see SECURITY.md): the file size and the number of features are capped; every field
// taken from a line is checked (coordinates whole numbers, end >= start, block lists consistent);
// a line that does not parse is skipped and counted, never trusted. Text from files is only ever
// drawn on the canvas or passed through esc() in cards.
'use strict';
if (typeof contigKey === 'undefined' && typeof require !== 'undefined') Object.assign(globalThis, require('./core.js'));

const TRACK_LIMITS = { bytes: 512 << 20, features: 5000000, blocks: 10000, samples: 100000 };
Object.assign(CGV_ERR_EN, {
  trackFormat: 'Unknown annotation format: {name}. LGV reads BED, bedGraph, VCF, GTF, GFF3 and TDF',
  trackEmpty: 'No features found in {name}',
  trackBig: '{name} is larger than LGV reads at once ({size})',
  trackMany: '{name} has more than {n} features',
});

// Format from the file name (a .gz ending is looked through).
function trackKind(name) {
  const n = String(name).toLowerCase().replace(/\.(gz|bgz)$/, '');
  if (/\.(bed|bed\d+|narrowpeak|broadpeak)$/.test(n)) return 'bed';
  if (/\.(bedgraph|bg)$/.test(n)) return 'bedgraph';
  if (/\.vcf$/.test(n)) return 'vcf';
  if (/\.gtf$/.test(n)) return 'gtf';
  if (/\.tdf$/.test(n)) return 'tdf';
  if (/\.(gff|gff3)$/.test(n)) return 'gff';
  return null;
}
const isTrackFile = name => !!trackKind(name);

// Whole, plain-number coordinate (no signs, decimals or exponents); -1 when it is not one.
const intField = s => /^\d{1,12}$/.test(s) ? +s : -1;

class FeatureIndex {
  constructor() { this.byChr = new Map(); this.count = 0; this.skipped = 0; }
  add(chr, f) {
    const k = contigKey(chr); let c = this.byChr.get(k);
    if (!c) this.byChr.set(k, c = { list: [], maxLen: 0 });
    c.list.push(f); if (f.end - f.start > c.maxLen) c.maxLen = f.end - f.start;
    if (++this.count > TRACK_LIMITS.features) throw new CgvError('trackMany', { name: this.name, n: TRACK_LIMITS.features.toLocaleString('en') });
  }
  finish() { for (const c of this.byChr.values()) c.list.sort((a, b) => a.start - b.start || a.end - b.end); return this; }
  // features overlapping [a, b) on a contig (any naming: chr7, 7)
  query(contig, a, b) {
    const c = this.byChr.get(contigKey(contig)); if (!c) return [];
    const L = c.list, from = a - c.maxLen; let lo = 0, hi = L.length;
    while (lo < hi) { const m = (lo + hi) >> 1; if (L[m].start < from) lo = m + 1; else hi = m; }
    const out = [];
    for (let i = lo; i < L.length && L[i].start < b; i++) if (L[i].end > a) out.push(L[i]);
    return out;
  }
}

function* lines(text) { let p = 0; while (p < text.length) { let e = text.indexOf('\n', p); if (e < 0) e = text.length; const l = text.charCodeAt(e - 1) === 13 ? text.slice(p, e - 1) : text.slice(p, e); p = e + 1; yield l; } }

// BED 3–12: chrom start end [name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts]
function parseBed(text, name) {
  const ix = new FeatureIndex(); ix.name = name;
  for (const line of lines(text)) {
    if (!line || line[0] === '#' || /^(track|browser)\b/.test(line)) continue;
    const f = line.includes('\t') ? line.split('\t') : line.trim().split(/\s+/);
    const start = intField(f[1]), end = intField(f[2]);
    if (f.length < 3 || !f[0] || start < 0 || end < start) { ix.skipped++; continue; }
    const ft = { start, end: Math.max(end, start + 1), name: f[3] && f[3] !== '.' ? f[3].slice(0, 200) : '', strand: f[5] === '-' ? -1 : f[5] === '+' ? 1 : 0, chr: f[0] };
    if (f[4] !== undefined && f[4] !== '.' && isFinite(+f[4])) ft.score = +f[4];
    const ts = intField(f[6]), te = intField(f[7]);
    if (ts >= start && te >= ts && te <= end) { ft.cdsStart = ts; ft.cdsEnd = te; }
    if (f[8] && /^\d{1,3},\d{1,3},\d{1,3}$/.test(f[8])) { const c = f[8].split(',').map(Number); if (c.every(v => v <= 255)) ft.rgb = c; }
    const n = intField(f[9]);
    if (n > 0 && n <= TRACK_LIMITS.blocks && f[10] && f[11]) {
      const sz = f[10].split(',').filter(x => x !== '').map(intField), st = f[11].split(',').filter(x => x !== '').map(intField);
      // blocks: in order, inside the feature, not overlapping
      let ok = sz.length === n && st.length === n, prev = 0;
      const ex = new Int32Array(n * 2);
      for (let i = 0; ok && i < n; i++) {
        const a = start + st[i], b = a + sz[i];
        if (sz[i] < 0 || st[i] < 0 || a < start + prev || b > end) ok = false;
        ex[2 * i] = a; ex[2 * i + 1] = b; prev = st[i] + sz[i];
      }
      if (ok) ft.exons = ex;
    }
    ix.add(f[0], ft);
  }
  return ix.finish();
}

// bedGraph: chrom start end value
function parseBedGraph(text, name) {
  const ix = new FeatureIndex(); ix.name = name; ix.min = Infinity; ix.max = -Infinity;
  for (const line of lines(text)) {
    if (!line || line[0] === '#' || /^(track|browser)\b/.test(line)) continue;
    const f = line.includes('\t') ? line.split('\t') : line.trim().split(/\s+/);
    const start = intField(f[1]), end = intField(f[2]), v = +f[3];
    if (f.length < 4 || !f[0] || start < 0 || end <= start || !isFinite(v) || f[3] === '') { ix.skipped++; continue; }
    ix.add(f[0], { start, end, value: v });
    if (v < ix.min) ix.min = v; if (v > ix.max) ix.max = v;
  }
  return ix.finish();
}

// VCF: sites with REF, ALT, QUAL, FILTER, INFO and the samples' fields (kept as text, read on demand).
function parseVcf(text, name) {
  const ix = new FeatureIndex(); ix.name = name; ix.samples = []; ix.header = [];
  for (const line of lines(text)) {
    if (!line) continue;
    if (line.startsWith('##')) { if (ix.header.length < 5000) ix.header.push(line.slice(0, 2000)); continue; }
    if (line[0] === '#') { const f = line.split('\t'); ix.samples = f.slice(9, 9 + TRACK_LIMITS.samples).map(s => s.slice(0, 200)); continue; }
    const f = line.split('\t'), pos = intField(f[1]);
    if (f.length < 8 || !f[0] || pos < 1 || !/^[A-Za-z*.]+$/.test(f[3] || '')) { ix.skipped++; continue; }
    const ref = f[3].toUpperCase(), start = pos - 1;
    ix.add(f[0], {
      start, end: start + Math.max(1, ref.length), id: f[2] === '.' ? '' : f[2].slice(0, 200), ref: ref.slice(0, 1000), alt: f[4].slice(0, 1000),
      qual: f[5], filter: f[6], info: f[7], format: f[8] || '', gts: f.length > 9 ? f.slice(9) : null,
    });
  }
  return ix.finish();
}
// Genotype class of sample i at a site: 0 no call, 1 hom ref, 2 het, 3 hom alt.
function vcfGenotype(site, i) {
  if (!site.gts || i >= site.gts.length) return 0;
  const k = site.format.split(':').indexOf('GT'); if (k < 0) return 0;
  const gt = (site.gts[i].split(':')[k] || '').split(/[\/|]/);
  if (!gt.length || gt.some(a => a === '.' || a === '')) return 0;
  const alt = gt.filter(a => a !== '0').length;
  return alt === 0 ? 1 : alt === gt.length ? 3 : 2;
}
function vcfInfo(site) {
  const out = [];
  if (!site.info || site.info === '.') return out;
  for (const kv of site.info.split(';')) { const e = kv.indexOf('='); out.push(e < 0 ? [kv, ''] : [kv.slice(0, e), kv.slice(e + 1)]); if (out.length > 200) break; }
  return out;
}

// GTF / GFF3: exons (and CDS) grouped into transcripts; files without exons show every line.
function parseGtfGff(text, name, kind) {
  const ix = new FeatureIndex(); ix.name = name;
  const tx = new Map(), plain = [], names = new Map();
  const attr = s => {
    const a = {};
    if (kind === 'gtf') { for (const m of s.matchAll(/(\w+)\s+"([^"]*)"/g)) if (!(m[1] in a)) a[m[1]] = m[2]; }
    else for (const kv of s.split(';')) { const e = kv.indexOf('='); if (e > 0) a[kv.slice(0, e).trim()] = decodeGff(kv.slice(e + 1).trim()); }
    return a;
  };
  for (const line of lines(text)) {
    if (!line || line[0] === '#') { if (line === '##FASTA') break; continue; }
    const f = line.split('\t');
    const s1 = intField(f[3]), e1 = intField(f[4]);
    if (f.length < 9 || !f[0] || s1 < 1 || e1 < s1) { ix.skipped++; continue; }
    const type = f[2], a = attr(f[8]), strand = f[6] === '-' ? -1 : f[6] === '+' ? 1 : 0, start = s1 - 1, end = e1;
    const child = /^(exon|CDS|five_prime_UTR|three_prime_UTR|UTR|start_codon|stop_codon)$/.test(type);
    if (child) {
      const parents = kind === 'gtf' ? [a.transcript_id] : String(a.Parent || '').split(',');
      for (const pid of parents) {
        if (!pid) continue;
        let t = tx.get(pid);
        if (!t) tx.set(pid, t = { chr: f[0], strand, exons: [], cds: [], start, end, name: '', id: pid, attrs: a, type: 'transcript' });
        if (type === 'exon') t.exons.push(start, end); else if (type === 'CDS') t.cds.push(start, end);
        t.start = Math.min(t.start, start); t.end = Math.max(t.end, end);
        if (kind === 'gtf' && !t.name) t.name = a.transcript_name || a.gene_name || '';
      }
    } else {
      plain.push({ chr: f[0], start, end, strand, type, name: (a.Name || a.gene_name || a.Name || a.ID || a.gene_id || '').slice(0, 200), attrs: a });
      if (a.ID) names.set(a.ID, (a.Name || a.gene_name || '').slice(0, 200));
    }
    if (tx.size + plain.length > TRACK_LIMITS.features) throw new CgvError('trackMany', { name, n: TRACK_LIMITS.features.toLocaleString('en') });
  }
  if (tx.size) {
    const inTx = new Set(tx.keys());
    for (const t of tx.values()) {
      if (!t.exons.length) t.exons = t.cds.slice();
      if (!t.exons.length) continue;
      const pairs = []; for (let i = 0; i < t.exons.length; i += 2) pairs.push([t.exons[i], t.exons[i + 1]]);
      pairs.sort((x, y) => x[0] - y[0]);
      const ex = new Int32Array(pairs.length * 2); pairs.forEach((p, i) => { ex[2 * i] = p[0]; ex[2 * i + 1] = p[1]; });
      const ft = { start: t.start, end: t.end, strand: t.strand, exons: ex, name: (t.name || names.get(t.id) || t.id).slice(0, 200), attrs: t.attrs, type: 'transcript' };
      if (t.cds.length) { ft.cdsStart = Math.min(...t.cds.filter((_, i) => i % 2 === 0)); ft.cdsEnd = Math.max(...t.cds.filter((_, i) => i % 2 === 1)); }
      ix.add(t.chr, ft);
    }
    // transcript and gene lines are shown through their exons; other lines (repeats, regions) as they are
    for (const p of plain) if (!inTx.has(p.attrs.ID) && !/^(gene|mRNA|transcript|ncRNA|lnc_RNA|pseudogene|primary_transcript|miRNA|snRNA|snoRNA|rRNA|tRNA)$/.test(p.type)) ix.add(p.chr, p);
  } else for (const p of plain) ix.add(p.chr, p);
  return ix.finish();
}
function decodeGff(s) { try { return decodeURIComponent(s); } catch (e) { return s; } }

// gzip or bgzip (blocks of gzip one after another, as tabix-ready files are), capped in size.
async function gunzipTrack(bytes, name) {
  const big = { name, size: '512 MB' };
  if (!isBgzfHeader(bytes)) return inflateStream(bytes, 'gzip', TRACK_LIMITS.bytes, 'inflate', big, 'trackBig');
  const parts = []; let total = 0;
  for (const b of parseBgzfBlocks(bytes, 0)) {
    const u = await inflateBgzf(bytes, b); total += u.length;
    if (total > TRACK_LIMITS.bytes) throw new CgvError('trackBig', big);
    parts.push(u);
  }
  const out = new Uint8Array(total); let o = 0; for (const u of parts) { out.set(u, o); o += u.length; }
  return out;
}
// ---- TDF (tiled data format: coverage at several zoom levels), read by region with Range reads ----
// Layout (little-endian): magic "TDF", version, index position and size, header (window functions,
// track type and line, track names, genome, flags: bit 0 = tiles compressed with zlib); the index
// lists datasets ("/chr1/z3/mean", "/chr1/raw") and groups ("/" with "chromosomes", "maxZoom").
// A dataset lists tiles; a tile holds fixedStep, variableStep or bed records with float values.
// Every count and size is checked against the bytes actually there before use.
const TDF_LIMITS = { header: 64 << 20, index: 64 << 20, dataset: 64 << 20, tile: 64 << 20, string: 100000, tiles: 10000000 };
Object.assign(CGV_ERR_EN, { tdfBad: 'The TDF file is damaged or not TDF: {what}' });
class BinReader {
  constructor(u8) { this.u8 = u8; this.dv = new DataView(u8.buffer, u8.byteOffset, u8.byteLength); this.p = 0; }
  need(n, what) { if (n < 0 || this.p + n > this.u8.length) throw new CgvError('tdfBad', { what: what || 'cut short' }); }
  int() { this.need(4); const v = this.dv.getInt32(this.p, true); this.p += 4; return v; }
  float() { this.need(4); const v = this.dv.getFloat32(this.p, true); this.p += 4; return v; }
  long() { this.need(8); const lo = this.dv.getUint32(this.p, true), hi = this.dv.getUint32(this.p + 4, true); this.p += 8; if (hi & 0x80000000) return -1; /* empty tile */ if (hi > 0x1fffff) throw new CgvError('tdfBad', { what: 'an offset is too large' }); return hi * 4294967296 + lo; }
  str() { let e = this.p; while (e < this.u8.length && this.u8[e] !== 0) e++; if (e >= this.u8.length || e - this.p > TDF_LIMITS.string) throw new CgvError('tdfBad', { what: 'a text field does not end' }); const s = new TextDecoder().decode(this.u8.subarray(this.p, e)); this.p = e + 1; return s; }
  count(each, max, what) { const n = this.int(); if (n < 0 || n > max || (each && this.p + n * each > this.u8.length)) throw new CgvError('tdfBad', { what }); return n; }
}
class TdfFile {
  constructor(src) { this.src = src; this.name = src.name; this.datasets = new Map(); this.tiles = new Map(); this.pending = new Set(); this.count = 1; this.skipped = 0; this.min = 0; }
  async read(pos, size, max) { if (size < 0 || size > max) throw new CgvError('tdfBad', { what: 'a block is too large' }); const b = await this.src.read(pos, size); if (b.length < size) throw new CgvError('tdfBad', { what: 'cut short' }); return b; }
  async open() {
    const h = await this.src.read(0, 64000), r = new BinReader(h);
    if (h.length < 24 || h[0] !== 84 || h[1] !== 68 || h[2] !== 70) throw new CgvError('tdfBad', { what: 'no TDF signature' });
    r.int(); this.version = r.int();
    const indexPos = r.long(), indexSize = r.int(); r.int();
    this.windowFunctions = [];
    if (this.version >= 2) { const n = r.count(1, 1000, 'window functions'); for (let i = 0; i < n; i++) this.windowFunctions.push(r.str()); }
    this.trackType = r.str(); this.trackLine = r.str();
    const nt = r.count(1, 100000, 'track names'); this.trackNames = []; for (let i = 0; i < nt; i++) this.trackNames.push(r.str());
    if (this.version >= 3) { this.genomeId = r.str(); this.compressed = (r.int() & 1) !== 0; } else this.compressed = false;
    const ix = new BinReader(await this.read(indexPos, indexSize, TDF_LIMITS.index));
    this.dsIndex = new Map(); this.grIndex = new Map();
    for (const map of [this.dsIndex, this.grIndex]) { const n = ix.count(13, TDF_LIMITS.tiles, 'index'); for (let i = 0; i < n; i++) { const k = ix.str(); map.set(k, { pos: ix.long(), size: ix.int() }); } }
    const root = await this.group('/');
    this.maxZoom = root && isFinite(+root.maxZoom) ? +root.maxZoom : -1;
    this.chrs = new Map();
    for (const c of String((root && root.chromosomes) || '').split(',')) if (c) this.chrs.set(contigKey(c), c);
    for (const k of this.dsIndex.keys()) { const c = k.split('/')[1]; if (c && !this.chrs.has(contigKey(c))) this.chrs.set(contigKey(c), c); }
    return this;
  }
  async group(name) {
    const e = this.grIndex.get(name); if (!e) return null;
    const r = new BinReader(await this.read(e.pos, e.size, TDF_LIMITS.header)), g = {};
    const n = r.count(2, 100000, 'group'); for (let i = 0; i < n; i++) { const k = r.str(); g[k] = r.str(); }
    return g;
  }
  async dataset(name) {
    if (this.datasets.has(name)) return this.datasets.get(name);
    const e = this.dsIndex.get(name); let ds = null;
    if (e) {
      const r = new BinReader(await this.read(e.pos, e.size, TDF_LIMITS.dataset));
      const na = r.count(2, 100000, 'dataset'); for (let i = 0; i < na; i++) { r.str(); r.str(); }
      r.str(); const tileWidth = r.float(), n = r.count(12, TDF_LIMITS.tiles, 'tiles');
      if (!(tileWidth > 0)) throw new CgvError('tdfBad', { what: 'tile width' });
      const tiles = []; for (let i = 0; i < n; i++) tiles.push({ pos: r.long(), size: r.int() });
      ds = { tileWidth, tiles };
    }
    this.datasets.set(name, ds); return ds;
  }
  async tile(dsName, i, t) {
    const raw = await this.read(t.pos, t.size, TDF_LIMITS.tile);
    let u8 = raw;
    if (this.compressed) u8 = await inflateStream(raw, 'deflate', TDF_LIMITS.tile, 'tdfBad', { what: 'a tile does not decompress' });
    const r = new BinReader(u8), type = r.str(), out = [];
    const floats = n => { const a = new Float32Array(n); for (let k = 0; k < n; k++) a[k] = r.float(); return a; };
    if (type === 'fixedStep') {
      const n = r.count(4, 1e8, 'tile'), start = r.int(), span = r.float(), v = floats(n);
      if (!(span > 0 && span < 1e9)) throw new CgvError('tdfBad', { what: 'tile step' });
      for (let k = 0; k < n; k++) if (!Number.isNaN(v[k])) out.push({ start: Math.round(start + k * span), end: Math.round(start + (k + 1) * span), value: v[k] });
    } else if (type === 'variableStep') {
      r.int(); const span = r.float(), n = r.count(8, 1e8, 'tile'), st = new Int32Array(n);
      if (!(span > 0 && span < 1e9)) throw new CgvError('tdfBad', { what: 'tile span' });
      for (let k = 0; k < n; k++) st[k] = r.int();
      r.int(); const v = floats(n);
      for (let k = 0; k < n; k++) if (!Number.isNaN(v[k])) out.push({ start: st[k], end: Math.round(st[k] + span), value: v[k] });
    } else if (type === 'bed' || type === 'bedWithName') {
      const n = r.count(12, 1e8, 'tile'), st = new Int32Array(n), en = new Int32Array(n);
      for (let k = 0; k < n; k++) st[k] = r.int(); for (let k = 0; k < n; k++) en[k] = r.int();
      r.int(); const v = floats(n);
      for (let k = 0; k < n; k++) if (!Number.isNaN(v[k]) && en[k] >= st[k]) out.push({ start: st[k], end: Math.max(en[k], st[k] + 1), value: v[k] });
    } else throw new CgvError('tdfBad', { what: 'unknown tile type' });
    return out;
  }
  // Values over [a, b) of a contig at the zoom fitting bpPerPx; tiles not read yet are fetched in the
  // background and onLoad() is called when they arrive. Returns what is already there.
  values(contig, a, b, bpPerPx, chrLen, onLoad) {
    const c = this.chrs.get(contigKey(contig)); if (!c) return [];
    let z = Math.ceil(Math.log2(Math.max(1, chrLen / (bpPerPx * 700))));
    if (!(z >= 0)) z = 0;
    const wf = this.windowFunctions.includes('mean') ? 'mean' : this.windowFunctions[0] || '';
    const name = z > this.maxZoom ? `/${c}/raw` : `/${c}/z${z}` + (this.version >= 2 ? '/' + wf : '');
    const ds = this.datasets.get(name);
    if (ds === undefined) { this.load(name, null, onLoad); return []; }
    if (!ds) return [];
    const out = [], t0 = Math.max(0, Math.floor(a / ds.tileWidth)), t1 = Math.min(ds.tiles.length - 1, Math.floor(b / ds.tileWidth));
    if (t1 - t0 > 2000) return out; // far too many tiles for one view: a coarser zoom should be used
    for (let i = t0; i <= t1; i++) {
      const t = ds.tiles[i]; if (!t || t.size <= 0 || t.pos < 0) continue;
      const v = this.tiles.get(name + '#' + i);
      if (v === undefined) { this.load(name, i, onLoad); continue; }
      for (const f of v) if (f.end > a && f.start < b) out.push(f);
    }
    return out;
  }
  load(name, i, onLoad) {
    const key = name + '#' + (i == null ? 'ds' : i);
    if (this.pending.has(key) || this.failed) return;
    this.pending.add(key);
    (i == null ? this.dataset(name) : this.dataset(name).then(ds => this.tile(name, i, ds.tiles[i])).then(v => { this.tiles.set(name + '#' + i, v); if (this.tiles.size > 4000) this.tiles.delete(this.tiles.keys().next().value); }))
      .catch(e => { this.failed = e; if (onLoad) onLoad(e); })
      .then(() => { this.pending.delete(key); if (onLoad && !this.failed) onLoad(); });
  }
}

// Reads a track file from a source (file or link): gzip / bgzip is unpacked, the size capped.
async function loadTrack(src) {
  const name = src.name || 'track', kind = trackKind(name);
  if (!kind) throw new CgvError('trackFormat', { name });
  if (kind === 'tdf') { const f = await new TdfFile(src).open(); f.kind = 'tdf'; return f; } // read by region, not whole
  let bytes = await src.readAll();
  if (bytes.length > TRACK_LIMITS.bytes) throw new CgvError('trackBig', { name, size: '512 MB' });
  if (bytes.length >= 2 && bytes[0] === 31 && bytes[1] === 139) bytes = await gunzipTrack(bytes, name);
  const text = new TextDecoder().decode(bytes);
  const ix = kind === 'bed' ? parseBed(text, name) : kind === 'bedgraph' ? parseBedGraph(text, name) : kind === 'vcf' ? parseVcf(text, name) : parseGtfGff(text, name, kind);
  if (!ix.count) throw new CgvError('trackEmpty', { name });
  ix.kind = kind;
  return ix;
}

if (typeof module !== 'undefined') module.exports = { TdfFile, gunzipTrack, TRACK_LIMITS, trackKind, isTrackFile, FeatureIndex, parseBed, parseBedGraph, parseVcf, vcfGenotype, vcfInfo, parseGtfGff, loadTrack };
