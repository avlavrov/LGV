// CGV: CRAM reading (CRAM 3.0 and 3.1). Written from the CRAM specification (samtools/hts-specs,
// CRAMv3.tex); no third-party code in this file.
//
// File definition, containers, blocks, the SAM header, the .crai index and records (CRAM 3.0 data
// series and encodings; record decoding follows the specification and is checked against htslib's
// output on the same data). Block codecs: raw and gzip so far.
//
// Safety rules for this file (see SECURITY.md):
//  - every read is bounds-checked through CramReader; a short or damaged file gives a CgvError,
//    never an exception from outside the reader or an endless loop;
//  - sizes taken from the file are checked against CRAM_LIMITS before anything is allocated;
//  - CRC32 checksums of container headers and blocks are verified (CRAM 3);
//  - nothing here touches the page, the network or storage: data comes in through a source object
//    (file or URL) chosen by the user.
'use strict';
if (typeof CgvError === 'undefined' && typeof require !== 'undefined') Object.assign(globalThis, require('./core.js'));
if (typeof ransDecode4x8 === 'undefined' && typeof require !== 'undefined') Object.assign(globalThis, require('./htscodecs.js'));
if (typeof bzip2Decode === 'undefined' && typeof require !== 'undefined') Object.assign(globalThis, require('./compress.js'));

const CRAM_LIMITS = {
  block: 256 << 20,       // decoded size of one block
  container: 1024 << 20,  // compressed size of one container
  blocks: 100000,         // blocks in one container
  header: 64 << 20,       // SAM header text
  index: 512 << 20,       // decoded .crai
  readLength: 1 << 23,    // bases in one read (ultra-long nanopore reads stay below)
  sliceBases: 256 << 20,  // bases in one slice: values that cost no input bits cannot inflate memory
};
// Checksums are verified always; the tests switch them off to send damaged data deeper into the decoder.
const CRAM_CHECK = { crc: true, md5: true };
const CRAM_METHODS = ['raw', 'gzip', 'bzip2', 'lzma', 'rANS 4x8', 'rANS 4x16', 'arithmetic', 'fqzcomp', 'name tokenizer'];
const CRAM_CT = { FILE_HEADER: 0, COMPRESSION_HEADER: 1, SLICE_HEADER: 2, EXTERNAL: 4, CORE: 5 };

Object.assign(CGV_ERR_EN, {
  notCram: 'This is not a CRAM file (no CRAM signature)',
  cramVersion: 'CRAM {v} is not supported; CGV reads CRAM 3.0 and 3.1',
  cramBad: 'The CRAM file is damaged: {what}',
  cramCrc: 'CRAM checksum mismatch in {what}: the file is damaged',
  cramCodec: 'A CRAM block is compressed with {codec}, which CGV does not read yet',
  cramLimit: 'A CRAM {what} is larger than CGV allows ({size})',
  notCrai: 'This is not a CRAM index (.crai)',
  cramRefMd5: 'The loaded reference is not the one this CRAM was made with: the MD5 of contig {name} differs. Load the reference the file was made with',
  cramNeedRef: 'This CRAM stores reads as differences from its reference: load the reference it was made with (contig {name})',
  cramEncoding: 'CRAM data series {ds} uses encoding {enc}, which CGV does not read',
});

// CRC32: crc32() from core.js (the same check as for BGZF blocks).

// ---------- MD5 (RFC 1321), for the reference check of each slice ----------
const MD5_K = new Uint32Array(64).map((_, i) => Math.floor(Math.abs(Math.sin(i + 1)) * 4294967296));
const MD5_S = [7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22, 5, 9, 14, 20, 5, 9, 14, 20, 5, 9, 14, 20, 5, 9, 14, 20, 4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23, 6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21];
function md5(data) {
  const n = data.length, total = ((n + 8) >> 6) + 1 << 6, m = new Uint8Array(total);
  m.set(data); m[n] = 0x80;
  const bits = n * 8; for (let i = 0; i < 8; i++) m[total - 8 + i] = Math.floor(bits / 2 ** (8 * i)) & 255;
  const w = new Uint32Array(16);
  let a0 = 0x67452301, b0 = 0xefcdab89, c0 = 0x98badcfe, d0 = 0x10325476;
  for (let off = 0; off < total; off += 64) {
    for (let i = 0; i < 16; i++) w[i] = m[off + 4 * i] | (m[off + 4 * i + 1] << 8) | (m[off + 4 * i + 2] << 16) | (m[off + 4 * i + 3] << 24);
    let A = a0, B = b0, C = c0, D = d0;
    for (let i = 0; i < 64; i++) {
      let F, g;
      if (i < 16) { F = (B & C) | (~B & D); g = i; }
      else if (i < 32) { F = (D & B) | (~D & C); g = (5 * i + 1) & 15; }
      else if (i < 48) { F = B ^ C ^ D; g = (3 * i + 5) & 15; }
      else { F = C ^ (B | ~D); g = (7 * i) & 15; }
      F = (F + A + MD5_K[i] + w[g]) | 0;
      A = D; D = C; C = B;
      B = (B + ((F << MD5_S[i]) | (F >>> (32 - MD5_S[i])))) | 0;
    }
    a0 = (a0 + A) | 0; b0 = (b0 + B) | 0; c0 = (c0 + C) | 0; d0 = (d0 + D) | 0;
  }
  const out = new Uint8Array(16);
  [a0, b0, c0, d0].forEach((v, i) => { for (let k = 0; k < 4; k++) out[4 * i + k] = (v >>> (8 * k)) & 255; });
  return out;
}

// ---------- Bounds-checked reader ----------
class CramReader {
  constructor(u8, pos, end, what) { this.u8 = u8; this.pos = pos || 0; this.end = end === undefined ? u8.length : end; this.what = what || 'data'; }
  need(n) { if (n < 0 || this.pos + n > this.end) throw new CgvError('cramBad', { what: this.what + ' is cut short' }); }
  byte() { this.need(1); return this.u8[this.pos++]; }
  i32() { this.need(4); const u = this.u8, p = this.pos; this.pos += 4; return (u[p] | u[p + 1] << 8 | u[p + 2] << 16 | u[p + 3] << 24); }
  u32() { return this.i32() >>> 0; }
  // ITF8: 1–5 bytes, the count given by leading 1 bits of the first byte; a signed 32-bit value
  itf8() {
    const b0 = this.byte();
    if (b0 < 0x80) return b0;
    if (b0 < 0xC0) return ((b0 & 0x3F) << 8) | this.byte();
    if (b0 < 0xE0) { this.need(2); const u = this.u8, p = this.pos; this.pos += 2; return ((b0 & 0x1F) << 16) | (u[p] << 8) | u[p + 1]; }
    if (b0 < 0xF0) { this.need(3); const u = this.u8, p = this.pos; this.pos += 3; return ((b0 & 0x0F) << 24) | (u[p] << 16) | (u[p + 1] << 8) | u[p + 2]; }
    this.need(4); const u = this.u8, p = this.pos; this.pos += 4;
    return ((b0 & 0x0F) << 28) | (u[p] << 20) | (u[p + 1] << 12) | (u[p + 2] << 4) | (u[p + 3] & 0x0F);
  }
  // LTF8: 1–9 bytes, up to 64 bits. Values beyond 2^53 cannot occur in a sane file and are refused.
  ltf8() {
    const b0 = this.byte();
    let n = 0; while (n < 8 && (b0 << n) & 0x80) n++;
    if (n === 0) return b0;
    this.need(n);
    let v = n < 8 ? b0 & (0xFF >> (n + 1)) : 0;
    for (let i = 0; i < n; i++) { v = v * 256 + this.u8[this.pos++]; if (v > Number.MAX_SAFE_INTEGER) throw new CgvError('cramBad', { what: this.what + ' has an impossible number' }); }
    return v;
  }
  itf8Array(max) {
    const n = this.itf8();
    if (n < 0 || n > max) throw new CgvError('cramBad', { what: this.what + ' has an impossible array length' });
    const a = new Array(n); for (let i = 0; i < n; i++) a[i] = this.itf8(); return a;
  }
  bytes(n) { this.need(n); const s = this.u8.subarray(this.pos, this.pos + n); this.pos += n; return s; }
}

// ---------- Containers and blocks ----------
// Container header; `off` is where it starts in `u8`. Returns null when `u8` ends before the header
// does (the caller reads more and tries again).
function parseContainerHeader(u8, off, major) {
  const r = new CramReader(u8, off, u8.length, 'a container header');
  try {
    const h = { offset: off };
    h.length = r.i32();
    h.refId = r.itf8(); h.start = r.itf8(); h.span = r.itf8(); h.nRecords = r.itf8();
    h.counter = major >= 3 ? r.ltf8() : r.itf8();
    h.bases = r.ltf8();
    h.nBlocks = r.itf8();
    h.landmarks = r.itf8Array(CRAM_LIMITS.blocks);
    if (major >= 3) { const end = r.pos; const crc = r.u32(); if (CRAM_CHECK.crc && crc32(u8, off, end) !== crc) throw new CgvError('cramCrc', { what: 'a container header' }); }
    h.headerSize = r.pos - off;
    if (h.length < 0 || h.length > CRAM_LIMITS.container) throw new CgvError('cramLimit', { what: 'container', size: fmtSize(h.length) });
    if (h.nBlocks < 0 || h.nBlocks > CRAM_LIMITS.blocks) throw new CgvError('cramBad', { what: 'a container has an impossible number of blocks' });
    // the end-of-file marker: an empty container on reference -1 at position 4542278
    h.eof = h.nRecords === 0 && h.refId === -1 && h.start === 4542278 && h.nBlocks === 1;
    return h;
  } catch (e) { if (e.code === 'cramBad' && /cut short/.test(e.message)) return null; throw e; }
}
// One block at `off` in `u8` (all of it must be there). Data stays compressed until decodeBlock.
function parseBlock(u8, off, major, end) {
  const r = new CramReader(u8, off, end === undefined ? u8.length : end, 'a block');
  const b = { offset: off };
  b.method = r.byte(); b.contentType = r.byte(); b.contentId = r.itf8();
  b.size = r.itf8(); b.rawSize = r.itf8();
  if (b.size < 0 || b.rawSize < 0) throw new CgvError('cramBad', { what: 'a block has a negative size' });
  if (b.rawSize > CRAM_LIMITS.block) throw new CgvError('cramLimit', { what: 'block', size: fmtSize(b.rawSize) });
  b.data = r.bytes(b.size);
  if (major >= 3) { const e = r.pos; const crc = r.u32(); if (CRAM_CHECK.crc && crc32(u8, off, e) !== crc) throw new CgvError('cramCrc', { what: 'a block' }); }
  b.next = r.pos;
  return b;
}
// gzip: the browser's own inflate, stopped as soon as the output passes the size the block declares.
function inflateGzip(u8, size) { return inflateStream(u8, 'gzip', size, 'cramBad', { what: 'a gzip block does not inflate to its declared size' }); }
// Decoded contents of a block. Codecs other than raw and gzip come in later stages.
// Where decoding time goes (milliseconds by codec, records, reference); shown in tests and diagnostics.
const CRAM_STATS = { codec: {}, records: 0, ref: 0, read: 0 };
const cramNow = () => (typeof performance !== 'undefined' ? performance.now() : Date.now());
async function decodeBlock(b) {
  const t0 = cramNow();
  const out = await decodeBlockRaw(b);
  const k = CRAM_METHODS[b.method] || b.method; CRAM_STATS.codec[k] = (CRAM_STATS.codec[k] || 0) + cramNow() - t0;
  return out;
}
async function decodeBlockRaw(b) {
  let out;
  if (b.method === 0) out = b.data;
  else if (b.method === 1) out = await inflateGzip(b.data, b.rawSize);
  else if (b.method === 2) out = bzip2Decode(b.data, b.rawSize);         // our own, src/compress.js
  else if (b.method === 3) out = xzDecode(b.data, b.rawSize);
  else if (b.method === 4) out = ransDecode4x8(b.data, b.rawSize); // htscodecs, see src/htscodecs.js
  else if (b.method === 5) out = ransDecode4x16(b.data, b.rawSize);   // CRAM 3.1 codecs, also from htscodecs
  else if (b.method === 6) out = arithDecode(b.data, b.rawSize);
  else if (b.method === 7) out = fqzDecode(b.data, b.rawSize);
  else if (b.method === 8) out = tok3Decode(b.data, b.rawSize);
  else throw new CgvError('cramCodec', { codec: CRAM_METHODS[b.method] || 'method ' + b.method });
  if (out.length !== b.rawSize) throw new CgvError('cramBad', { what: 'a block decodes to ' + out.length + ' bytes instead of ' + b.rawSize });
  return out;
}
function fmtSize(n) { return n >= 1 << 20 ? Math.round(n / (1 << 20)) + ' MB' : n + ' bytes'; }

// ---------- SAM header text ----------
function parseSamHeaderText(text) {
  const refs = [], refIndex = new Map(), samples = [];
  for (const line of text.split('\n')) {
    if (line.startsWith('@SQ')) {
      const sn = line.match(/\tSN:([^\t\r]+)/), ln = line.match(/\tLN:(\d+)/);
      if (!sn) continue;
      refIndex.set(sn[1], refs.length); refs.push({ name: sn[1], len: ln ? +ln[1] : 0 });
    } else if (line.startsWith('@RG')) {
      const m = line.match(/\tSM:([^\t\r]+)/); if (m && !samples.includes(m[1])) samples.push(m[1]);
    }
  }
  return { headerText: text, refs, refIndex, samples };
}

// ---------- .crai ----------
// Gzipped text, one line per slice: reference id, alignment start (1-based), span, container offset,
// slice offset within the container's data, slice size. Positions are kept 0-based here.
async function parseCrai(raw) {
  let u8 = raw;
  if (raw.length >= 2 && raw[0] === 0x1F && raw[1] === 0x8B) u8 = await inflateGzip(raw, CRAM_LIMITS.index);
  if (u8.length > CRAM_LIMITS.index) throw new CgvError('cramLimit', { what: 'index', size: fmtSize(u8.length) });
  const text = new TextDecoder().decode(u8);
  const byRef = new Map(); let n = 0;
  for (const line of text.split('\n')) {
    if (!line.trim()) continue;
    const f = line.split('\t');
    if (f.length < 6 || !f.every(x => /^-?\d+$/.test(x.trim()))) throw new CgvError('notCrai');
    const e = { refId: +f[0], start: Math.max(0, +f[1] - 1), span: +f[2], container: +f[3], slice: +f[4], size: +f[5] };
    if (!byRef.has(e.refId)) byRef.set(e.refId, []);
    byRef.get(e.refId).push(e); n++;
  }
  for (const list of byRef.values()) list.sort((a, b) => a.start - b.start || a.container - b.container);
  return { byRef, slices: n };
}

// ---------- CRAM file ----------
// Same outward shape as BamFile: refs, refIndex, samples, headerText, maxSpan, refId(),
// estimateBytes() and fetch().
class CramFile {
  constructor(src, craiSrc) { this.src = src; this.baiSrc = craiSrc; this.maxSpan = 0; this.isCram = true; }
  async open() {
    const [def, idx] = await Promise.all([this.src.read(0, 26), this._readIndex()]);
    if (def.length < 26 || def[0] !== 67 || def[1] !== 82 || def[2] !== 65 || def[3] !== 77) throw new CgvError('notCram');
    this.major = def[4]; this.minor = def[5]; this.version = this.major + '.' + this.minor;
    if (this.major !== 3) throw new CgvError('cramVersion', { v: this.version });
    this.index = await parseCrai(idx);
    Object.assign(this, await this._readHeader(26));
    return this;
  }
  async _readIndex() {
    const list = Array.isArray(this.baiSrc) ? this.baiSrc : [this.baiSrc];
    let first = null;
    for (const s of list) {
      try { const b = await s.readAll(); this.baiSrc = s; return b; }
      catch (e) {
        if (list.length === 1 || !['notfound', 'auth', 'http', 'network'].includes(e.code)) throw e;
        if (!first || first.code === 'notfound') first = e;
      }
    }
    if (first.code === 'notfound' || first.code === 'auth') throw new CgvError('noIndex', { name: this.src.name, tried: list.map(s => s.name).join(', ') });
    throw first;
  }
  // The whole container at `off`: { header, data (its blocks' bytes), dataOffset }.
  async readContainer(off) {
    let size = 512, raw, h;
    for (;;) { // the header is small; read a bit more if it did not fit
      raw = await this.src.read(off, size);
      h = parseContainerHeader(raw, 0, this.major);
      if (h) break;
      if (raw.length < size) throw new CgvError('cramBad', { what: 'the file ends inside a container header' });
      size *= 4; if (size > 1 << 20) throw new CgvError('cramBad', { what: 'a container header is impossibly long' });
    }
    const dataOffset = off + h.headerSize;
    let data = raw.subarray(h.headerSize, h.headerSize + h.length);
    if (data.length < h.length) { const t = cramNow(); data = await this.src.read(dataOffset, h.length); CRAM_STATS.read += cramNow() - t; }
    if (data.length < h.length) throw new CgvError('cramBad', { what: 'the file ends inside a container' });
    h.offset = off;
    return { header: h, data, dataOffset };
  }
  // Blocks of a container's data, in order.
  blocksOf(c, from, to) {
    const out = []; let p = from || 0; const end = to === undefined ? c.data.length : to;
    while (p < end) { const b = parseBlock(c.data, p, this.major, end); out.push(b); p = b.next; if (out.length > CRAM_LIMITS.blocks) throw new CgvError('cramBad', { what: 'too many blocks' }); }
    return out;
  }
  async _readHeader(off) {
    const c = await this.readContainer(off);
    const blocks = this.blocksOf(c, 0, c.data.length);
    const hb = blocks.find(b => b.contentType === CRAM_CT.FILE_HEADER);
    if (!hb) throw new CgvError('cramBad', { what: 'no SAM header block' });
    const d = await decodeBlock(hb);
    const r = new CramReader(d, 0, d.length, 'the SAM header');
    const len = r.i32();
    if (len < 0 || len > CRAM_LIMITS.header) throw new CgvError('cramLimit', { what: 'SAM header', size: fmtSize(len) });
    const text = new TextDecoder().decode(r.bytes(len)).replace(/\0+$/, '');
    this.firstContainer = off + c.header.headerSize + c.header.length;
    return parseSamHeaderText(text);
  }
  refId(name) {
    if (this.refIndex.has(name)) return this.refIndex.get(name);
    if (!this._keys) { this._keys = new Map(); this.refs.forEach((r, i) => { const k = contigKey(r.name); if (!this._keys.has(k)) this._keys.set(k, i); }); }
    const id = this._keys.get(contigKey(name));
    return id === undefined ? -1 : id;
  }
  // Index entries of slices overlapping [beg, end) on refId (multi-reference slices, id -2, included).
  slices(refId, beg, end) {
    const out = [];
    for (const id of [refId, -2]) for (const e of this.index.byRef.get(id) || []) if (id === -2 || (e.start < end && e.start + e.span > beg)) out.push(e);
    return out;
  }
  estimateBytes(refId, beg, end) { let n = 0; for (const e of this.slices(refId, beg, end)) n += e.size; return n; }
  // Walks every container from the first one to the end-of-file marker, checking headers, blocks and
  // checksums. Used by the tests; a check of the file's structure, not needed for viewing.
  async scan(limit) {
    const out = []; let off = this.firstContainer;
    for (;;) {
      const c = await this.readContainer(off);
      const blocks = this.blocksOf(c);
      if (blocks.length !== c.header.nBlocks) throw new CgvError('cramBad', { what: 'a container holds ' + blocks.length + ' blocks instead of ' + c.header.nBlocks });
      out.push({ header: c.header, methods: [...new Set(blocks.map(b => CRAM_METHODS[b.method] || b.method))] });
      if (c.header.eof || (limit && out.length >= limit)) break;
      off += c.header.headerSize + c.header.length;
    }
    return out;
  }
  // Reference bases come from the page: refProvider(refId, start0, end0) → uppercase bases or null.
  // Decoded slices are kept for a while, so panning around reuses them.
  _cacheGet(k) { const c = this._cache || (this._cache = new Map()); const v = c.get(k); if (v) { c.delete(k); c.set(k, v); } return v; }
  _cachePut(k, v) { const c = this._cache || (this._cache = new Map()); c.set(k, v); while (c.size > 12) c.delete(c.keys().next().value); }
  async _container(off) {
    let e = this._cacheGet('c' + off);
    if (!e) {
      const c = await this.readContainer(off);
      const blocks = this.blocksOf(c, 0, c.header.landmarks.length ? c.header.landmarks[0] : c.data.length);
      const chb = blocks.find(b => b.contentType === CRAM_CT.COMPRESSION_HEADER);
      if (!chb) throw new CgvError('cramBad', { what: 'a container has no compression header' });
      e = { c, ch: parseCompressionHeader(await decodeBlock(chb)) };
      this._cachePut('c' + off, e);
    }
    return e;
  }
  // All records of slice `si` of the container at `off`.
  async slice(off, si) {
    const key = 's' + off + ':' + si;
    let recs = this._cacheGet(key);
    if (recs) return recs;
    const { c, ch } = await this._container(off);
    const lm = c.header.landmarks;
    if (si < 0 || si >= lm.length) throw new CgvError('cramBad', { what: 'the index points to a slice the container does not have' });
    const end = si + 1 < lm.length ? lm[si + 1] : c.data.length;
    recs = await decodeSlice(this, c, ch, lm[si], end, off);
    this._cachePut(key, recs);
    return recs;
  }
  // Same contract as BamFile.fetch: reads overlapping [beg, end) on refId, sorted by position, with
  // load-time filters and downsampling decided the same way, and coverage counted from every read.
  async fetch(refId, beg, end, opts) {
    opts = opts || {};
    const flt = opts.filters || null;
    const ds = opts.downsample && opts.downsample.max > 0 && opts.downsample.window > 0 ? opts.downsample : null;
    const cov = new SparseCoverage(beg, end);
    const depth = opts.depthOnly ? new DepthDiff(beg, end) : null;
    const out = { reads: [], cov, total: 0, sampled: [] };
    let win = -1, n = 0, kept = [], rng = 1;
    const flush = () => {
      if (kept.length) { kept.sort((a, b) => a.pos - b.pos); for (const r of kept) out.reads.push(r); }
      if (ds && n > ds.max) out.sampled.push([win * ds.window, (win + 1) * ds.window]);
      kept = []; n = 0;
    };
    // the slices the index lists for this window, in file order, each once
    const seen = new Set(), list = [];
    for (const e of this.slices(refId, beg, end)) { const k = e.container + ':' + e.slice; if (!seen.has(k)) { seen.add(k); list.push(e); } }
    list.sort((a, b) => a.container - b.container || a.slice - b.slice);
    for (const e of list) {
      const { c } = await this._container(e.container);
      const si = c.header.landmarks.indexOf(e.slice);
      for (const rec of await this.slice(e.container, si)) {
        if (rec.refID !== refId || rec.flag & 4 && rec.pos < 0) continue;
        if (rec.end - rec.pos > this.maxSpan) this.maxSpan = rec.end - rec.pos;
        if (rec.end <= beg || rec.pos >= end) continue;
        if (rec.flag & 4) continue; // unmapped reads placed next to their mate are not drawn
        if (flt && !passFlags(rec.flag, rec.mapq, flt)) continue;
        if (depth) depth.addRec(rec); else addCoverageRec(rec, cov);
        if (!ds || rec.pos >= beg) out.total++;
        if (!ds) { out.reads.push(rec); continue; }
        if (rec.pos < beg) continue; // as in BamFile.fetch: the left edge belongs to a wider query
        const w = Math.floor(rec.pos / ds.window);
        if (w !== win) { flush(); win = w; rng = (Math.imul(w + 1, 2654435761) >>> 0) || 1; }
        n++;
        if (kept.length < ds.max) kept.push(rec);
        else { rng ^= rng << 13; rng >>>= 0; rng ^= rng >>> 17; rng ^= rng << 5; rng >>>= 0; const j = rng % n; if (j < ds.max) kept[j] = rec; }
      }
    }
    flush();
    if (depth) depth.finish(cov);
    if (!ds) out.reads.sort((a, b) => a.pos - b.pos);
    return out;
  }
}

// ---------- Encodings (how each data series is stored) ----------
// Core data is a bit stream, read most significant bit first.
class BitReader {
  constructor(u8) { this.u8 = u8; this.byte = 0; this.bit = 7; }
  one() {
    if (this.byte >= this.u8.length) throw new CgvError('cramBad', { what: 'core data is cut short' });
    const v = (this.u8[this.byte] >> this.bit) & 1;
    if (--this.bit < 0) { this.bit = 7; this.byte++; }
    return v;
  }
  bits(n) { let v = 0; for (let i = 0; i < n; i++) v = v * 2 + this.one(); return v; }
}
const DS_INT = new Set(['BF', 'CF', 'RI', 'RL', 'AP', 'RG', 'MF', 'NS', 'NP', 'TS', 'NF', 'TL', 'FN', 'FP', 'DL', 'RS', 'PD', 'HC', 'MQ']);
const DS_BYTE = new Set(['FC', 'BA', 'QS', 'BS']);
const DS_ARRAY = new Set(['RN', 'IN', 'SC', 'BB', 'QQ']);
const ENC_NAMES = ['NULL', 'EXTERNAL', 'GOLOMB', 'HUFFMAN', 'BYTE_ARRAY_LEN', 'BYTE_ARRAY_STOP', 'BETA', 'SUBEXP', 'GOLOMB_RICE', 'GAMMA'];

// A decoder for one data series. Each has int(s), byte(s) and bytes(s, n) where s is the slice state
// { core: BitReader, ext: Map(content id → CramReader) }; n is the count for fixed-length arrays and
// is omitted for arrays that carry their own length.
function readEncoding(r, ds) {
  const id = r.itf8(), len = r.itf8();
  if (len < 0 || len > 1 << 20) throw new CgvError('cramBad', { what: 'an encoding has an impossible size' });
  const p = new CramReader(r.bytes(len), 0, len, 'encoding parameters');
  const ext = (s, cid) => { const b = s.ext.get(cid); if (!b) throw new CgvError('cramBad', { what: 'data series ' + ds + ' points to a missing block ' + cid }); return b; };
  const many = f => (s, n) => { if (!(n >= 0) || n > CRAM_LIMITS.block) throw new CgvError('cramBad', { what: 'an impossible array length in ' + ds }); const o = new Uint8Array(n); for (let i = 0; i < n; i++) o[i] = f(s); return o; };
  switch (id) {
    case 0: return { int: () => 0, byte: () => 0, bytes: (s, n) => new Uint8Array(n || 0) };
    case 1: { // EXTERNAL: integers as ITF8, bytes as bytes, in a block of its own
      const cid = p.itf8();
      return { int: s => ext(s, cid).itf8(), byte: s => ext(s, cid).byte(), bytes: (s, n) => { if (n === undefined) throw new CgvError('cramBad', { what: ds + ' has no length' }); return n === 0 ? new Uint8Array(0) : ext(s, cid).bytes(n); } };
    }
    case 3: { // HUFFMAN, canonical codes
      const syms = p.itf8Array(1 << 16), lens = p.itf8Array(1 << 16);
      if (syms.length !== lens.length || !syms.length) throw new CgvError('cramBad', { what: 'a Huffman table of ' + ds + ' is inconsistent' });
      if (syms.length === 1 && lens[0] === 0) { const v = syms[0]; return { int: () => v, byte: () => v, bytes: (s, n) => new Uint8Array(n).fill(v) }; }
      const codes = syms.map((sym, i) => ({ sym, len: lens[i] })).sort((a, b) => a.len - b.len || a.sym - b.sym);
      if (codes.some(c => c.len < 1 || c.len > 31)) throw new CgvError('cramBad', { what: 'a Huffman code of ' + ds + ' has an impossible length' });
      // per length: first code, its index in `codes`, count
      const first = new Map(); let code = 0, prevLen = codes[0].len;
      codes.forEach((c, i) => { if (i) { code++; if (c.len > prevLen) { code *= 2 ** (c.len - prevLen); prevLen = c.len; } } if (!first.has(c.len)) first.set(c.len, { code, index: i, count: 0 }); first.get(c.len).count++; });
      const maxLen = codes[codes.length - 1].len;
      const dec = s => {
        let v = 0;
        for (let len = 1; len <= maxLen; len++) {
          v = v * 2 + s.core.one();
          const f = first.get(len);
          if (f && v >= f.code && v < f.code + f.count) return codes[f.index + v - f.code].sym;
        }
        throw new CgvError('cramBad', { what: 'a Huffman code of ' + ds + ' is not in its table' });
      };
      return { int: dec, byte: dec, bytes: many(dec) };
    }
    case 4: { // BYTE_ARRAY_LEN: a length, then that many bytes
      const lenEnc = readEncoding(p, ds), valEnc = readEncoding(p, ds);
      return { bytes: s => { const n = lenEnc.int(s); if (n < 0 || n > CRAM_LIMITS.block) throw new CgvError('cramBad', { what: 'an impossible length in ' + ds }); return valEnc.bytes(s, n); } };
    }
    case 5: { // BYTE_ARRAY_STOP: bytes up to a stop byte
      const stop = p.byte(), cid = p.itf8();
      return { bytes: s => { const b = ext(s, cid), st = b.pos; let e = st; while (e < b.end && b.u8[e] !== stop) e++; if (e >= b.end) throw new CgvError('cramBad', { what: ds + ' has no stop byte' }); b.pos = e + 1; return b.u8.subarray(st, e); } };
    }
    case 6: { // BETA: fixed number of bits
      const off = p.itf8(), nb = p.itf8();
      if (nb < 0 || nb > 32) throw new CgvError('cramBad', { what: 'BETA width of ' + ds });
      const dec = s => s.core.bits(nb) - off;
      return { int: dec, byte: dec, bytes: many(dec) };
    }
    case 7: { // SUBEXP
      const off = p.itf8(), k = p.itf8();
      if (k < 0 || k > 31) throw new CgvError('cramBad', { what: 'SUBEXP of ' + ds });
      const dec = s => { let u = 0; while (s.core.one()) { if (++u > 31 - k) throw new CgvError('cramBad', { what: 'SUBEXP value of ' + ds }); } return (u ? 2 ** (u + k - 1) + s.core.bits(u + k - 1) : s.core.bits(k)) - off; };
      return { int: dec, byte: dec, bytes: many(dec) };
    }
    case 9: { // GAMMA
      const off = p.itf8();
      const dec = s => { let nz = 0; while (!s.core.one()) { if (++nz > 31) throw new CgvError('cramBad', { what: 'GAMMA value of ' + ds }); } return 2 ** nz + s.core.bits(nz) - off; };
      return { int: dec, byte: dec, bytes: many(dec) };
    }
    default: throw new CgvError('cramEncoding', { ds, enc: ENC_NAMES[id] || String(id) });
  }
}

// ---------- Compression header ----------
const SUBST_BASES = [65, 67, 71, 84, 78]; // A C G T N
const CRAM_BASE_IDX = new Uint8Array(256).fill(4); CRAM_BASE_IDX[65] = 0; CRAM_BASE_IDX[67] = 1; CRAM_BASE_IDX[71] = 2; CRAM_BASE_IDX[84] = 3; CRAM_BASE_IDX[97] = 0; CRAM_BASE_IDX[99] = 1; CRAM_BASE_IDX[103] = 2; CRAM_BASE_IDX[116] = 3;
function parseCompressionHeader(d) {
  const r = new CramReader(d, 0, d.length, 'the compression header');
  const h = { RN: true, AP: true, RR: true, QO: true, subst: null, TD: [], ds: {}, tags: new Map() };
  r.itf8(); // size of the preservation map
  let n = r.itf8();
  if (n < 0 || n > 64) throw new CgvError('cramBad', { what: 'the preservation map' });
  for (let i = 0; i < n; i++) {
    const key = String.fromCharCode(r.byte(), r.byte());
    if (key === 'RN' || key === 'AP' || key === 'RR' || key === 'QO') h[key] = !!r.byte();
    else if (key === 'SM') {
      const sm = r.bytes(5);
      // for each reference base, the 2-bit code of every other base in A C G T N order
      h.subst = SUBST_BASES.map((rb, ri) => { const t = new Uint8Array(4), others = SUBST_BASES.filter((_, j) => j !== ri); others.forEach((b, j) => { t[(sm[ri] >> (6 - 2 * j)) & 3] = b; }); return t; });
    } else if (key === 'TD') {
      const len = r.itf8(); if (len < 0 || len > CRAM_LIMITS.block) throw new CgvError('cramBad', { what: 'the tag dictionary' });
      const td = r.bytes(len); let line = [], st = 0;
      for (let j = 0; j < td.length; j++) if (td[j] === 0) { for (let k = st; k + 3 <= j; k += 3) line.push([td[k], td[k + 1], td[k + 2]]); h.TD.push(line); line = []; st = j + 1; }
    } else throw new CgvError('cramBad', { what: 'an unknown preservation key ' + key });
  }
  r.itf8(); n = r.itf8();
  if (n < 0 || n > 64) throw new CgvError('cramBad', { what: 'the data series map' });
  for (let i = 0; i < n; i++) { const key = String.fromCharCode(r.byte(), r.byte()); h.ds[key] = readEncoding(r, key); }
  r.itf8(); n = r.itf8();
  if (n < 0 || n > 1 << 16) throw new CgvError('cramBad', { what: 'the tag encoding map' });
  for (let i = 0; i < n; i++) { const key = r.itf8(); h.tags.set(key, readEncoding(r, 'tag ' + String.fromCharCode((key >> 16) & 255, (key >> 8) & 255))); }
  if (!h.subst) h.subst = SUBST_BASES.map(() => new Uint8Array([78, 78, 78, 78]));
  return h;
}

// ---------- Slices and records ----------
const NT16 = new Uint8Array(256).fill(15);
'=ACMGRSVTWYHKDBN'.split('').forEach((c, i) => { NT16[c.charCodeAt(0)] = i; NT16[c.toLowerCase().charCodeAt(0)] = i; });
const TAG_SIZE = { 65: 1, 99: 1, 67: 1, 115: 2, 83: 2, 105: 4, 73: 4, 102: 4 }; // A c C s S i I f

// Reads one slice into records shaped like parseBamRecord's.
async function decodeSlice(file, c, ch, from, to, containerOffset) {
  const major = file.major;
  const blocks = file.blocksOf(c, from, to);
  if (!blocks.length || blocks[0].contentType !== CRAM_CT.SLICE_HEADER) throw new CgvError('cramBad', { what: 'a slice does not start with its header' });
  const r = new CramReader(await decodeBlock(blocks[0]), 0, undefined, 'a slice header');
  const sh = { refId: r.itf8(), start: r.itf8(), span: r.itf8(), nRecords: r.itf8(), counter: r.ltf8(), nBlocks: r.itf8() };
  sh.ids = r.itf8Array(CRAM_LIMITS.blocks); sh.embedded = r.itf8(); sh.md5 = r.bytes(16);
  if (sh.nRecords < 0 || sh.nRecords > 10000000) throw new CgvError('cramBad', { what: 'a slice has an impossible number of records' });
  const s = { core: new BitReader(new Uint8Array(0)), ext: new Map() };
  for (const b of blocks.slice(1, 1 + sh.nBlocks)) {
    const d = await decodeBlock(b);
    if (b.contentType === CRAM_CT.CORE) s.core = new BitReader(d);
    else if (b.contentType === CRAM_CT.EXTERNAL) s.ext.set(b.contentId, new CramReader(d, 0, d.length, 'block ' + b.contentId));
  }
  // reference: embedded in the slice, or from the page for each reference the slice covers
  const refs = new Map(); // refId → { start0, bases }
  const needRef = ch.RR;
  if (sh.embedded >= 0) {
    const b = s.ext.get(sh.embedded); if (!b) throw new CgvError('cramBad', { what: 'the embedded reference block is missing' });
    refs.set(sh.refId, { start0: sh.start - 1, bases: b.u8.slice(0, b.end).map(x => x & ~32) });
  } else if (needRef) {
    const spans = sh.refId >= 0 ? [[sh.refId, sh.start - 1, sh.start - 1 + sh.span]]
      : sh.refId === -2 ? (file.index.byRef ? [...file.index.byRef.entries()].flatMap(([id, list]) => list.filter(e => e.container === containerOffset && e.slice === from).map(e => [id, e.start, e.start + e.span])).filter(x => x[0] >= 0) : []) : [];
    for (const [id, a, b] of spans) {
      const tr = cramNow(); const bases = file.refProvider ? await file.refProvider(id, a, b) : null; CRAM_STATS.ref += cramNow() - tr;
      if (!bases) throw new CgvError('cramNeedRef', { name: (file.refs[id] || {}).name || id });
      // the slice records the MD5 of the reference stretch it was written against: a different
      // reference (another build, another patch release) would silently give wrong bases
      if (sh.refId >= 0 && CRAM_CHECK.md5 && sh.md5.some(x => x) && bases.length === b - a) {
        const got = md5(bases); CRAM_STATS.md5 = (CRAM_STATS.md5 || 0) + 1;
        if (!got.every((x, k) => x === sh.md5[k])) throw new CgvError('cramRefMd5', { name: (file.refs[id] || {}).name || id });
      }
      refs.set(id, { start0: a, bases });
    }
  }
  // one reference in the slice (the usual case) is read without a lookup per base
  const only = refs.size === 1 ? [...refs.values()][0] : null, onlyId = only ? [...refs.keys()][0] : -9;
  const refAt = (id, p0) => {
    const R = id === onlyId ? only : refs.get(id); if (!R) return 78;
    const i = p0 - R.start0; return i >= 0 && i < R.bases.length ? R.bases[i] : 78;
  };
  // decoders of the data series, looked up once; a missing one fails only when a record needs it
  const E = {};
  for (const k of ['BF', 'CF', 'RI', 'RL', 'AP', 'RG', 'RN', 'MF', 'NS', 'NP', 'TS', 'NF', 'TL', 'FN', 'FC', 'FP', 'SC', 'BS', 'DL', 'IN', 'BA', 'QS', 'HC', 'PD', 'RS', 'BB', 'QQ', 'MQ']) {
    E[k] = ch.ds[k] || { int: () => { throw new CgvError('cramBad', { what: 'data series ' + k + ' has no encoding' }); } };
    if (!E[k].byte) E[k].byte = E[k].int; if (!E[k].bytes) E[k].bytes = E[k].int;
  }
  const rgBytes = file._rgBytes || (file._rgBytes = [...file.headerText.matchAll(/^@RG\t(?:.*\t)?ID:([^\t\r\n]+)/gm)].map(m => { const id = new TextEncoder().encode(m[1]); const b = new Uint8Array(id.length + 4); b.set([82, 71, 90]); b.set(id, 3); return b; }));
  // tag lines resolved to their encodings once per slice
  const tagLines = ch.TD.map(line => line.map(([a, b, t]) => {
    const enc = ch.tags.get((a << 16) | (b << 8) | t);
    return { a, b, t, enc, size: TAG_SIZE[t] || 0, self: !!enc && !!enc.bytes && !enc.int, keep: !(a === 99 && b === 70 && t === 67) && t !== 42 };
  }));
  const bad = what => new CgvError('cramBad', { what });
  const N = sh.nRecords, recs = new Array(N);
  const mateLine = new Int32Array(N).fill(-1), tlenOf = new Float64Array(N).fill(NaN), explicitTlen = new Float64Array(N).fill(NaN);
  let seq = new Uint8Array(512); // scratch for bases as letters, reused by every record
  let lastPos = sh.start, bases = 0; const tRec = cramNow();
  for (let i = 0; i < N; i++) {
    const rec = { refID: sh.refId, pos: -1, end: 0, name: '', mapq: 0, flag: 0, lseq: 0, mateRef: -1, matePos: -1, tlen: 0, cigar: null, seq: null, qual: null, tagBytes: null, _tags: null, feats: null };
    const flag = E.BF.int(s), cf = E.CF.int(s);
    rec.flag = flag & 0xFFFF;
    if (sh.refId === -2) rec.refID = E.RI.int(s);
    const rl = E.RL.int(s);
    if (rl < 0 || rl > CRAM_LIMITS.readLength) throw bad('a record has an impossible read length');
    if ((bases += rl) > CRAM_LIMITS.sliceBases) throw new CgvError('cramLimit', { what: 'slice', size: fmtSize(bases) + ' of bases' });
    let apos = E.AP.int(s);
    if (ch.AP) { apos += lastPos; lastPos = apos; }
    rec.pos = apos - 1;
    const rg = E.RG.int(s);
    if (ch.RN) rec.name = latin1Str(E.RN.bytes(s));
    if (cf & 2) { // mate kept separately (other slice, other reference, or not reproducible)
      const mf = E.MF.int(s);
      if (!ch.RN) rec.name = latin1Str(E.RN.bytes(s));
      rec.mateRef = E.NS.int(s); rec.matePos = E.NP.int(s) - 1; tlenOf[i] = E.TS.int(s);
      if (mf & 1) rec.flag |= 0x1 | 0x20;
      if (mf & 2) rec.flag |= 0x8;
      if (!(rec.flag & 1)) rec.mateRef = -1;
    } else if (cf & 4) { // mate further down in this slice
      mateLine[i] = i + 1 + E.NF.int(s);
      if (cf & 16) explicitTlen[i] = E.TS.int(s);
    } else if (cf & 16) explicitTlen[i] = E.TS.int(s);
    // tags, in the order of the tag line, as BAM writes them
    const tl = E.TL.int(s);
    if (tl < 0 || tl >= tagLines.length) throw bad('a record points to a missing tag line');
    const line = tagLines[tl], vals = new Array(line.length);
    let tn = 0;
    for (let k = 0; k < line.length; k++) {
      const T = line[k];
      if (!T.enc) throw bad('tag ' + String.fromCharCode(T.a, T.b) + ' has no encoding');
      const v = T.self ? T.enc.bytes(s) : T.size ? T.enc.bytes(s, T.size) : T.enc.bytes(s); // Z, H, B: self-delimiting
      if (!T.keep) continue; // cF (CRAM's own flags) and MD* / NM* (regenerated by htslib on request)
      vals[k] = v; tn += 3 + v.length;
    }
    const rgb = rg >= 0 && rg < rgBytes.length ? rgBytes[rg] : null;
    if (rgb) tn += rgb.length;
    const tagBytes = new Uint8Array(tn); let to = 0;
    for (let k = 0; k < line.length; k++) { const v = vals[k]; if (!v) continue; const T = line[k]; tagBytes[to] = T.a; tagBytes[to + 1] = T.b; tagBytes[to + 2] = T.t; tagBytes.set(v, to + 3); to += 3 + v.length; }
    if (rgb) tagBytes.set(rgb, to);
    rec.tagBytes = tagBytes;
    // sequence, qualities and CIGAR
    if (seq.length < rl) seq = new Uint8Array(rl * 2);
    const qual = new Uint8Array(rl).fill(255);
    // a CIGAR operation holds at most 2^28 - 1: longer values from a damaged file would wrap around
    const cig = []; const op = (o, n) => { if (n <= 0) return; if (n > 0xFFFFFFF) throw bad('an alignment operation is too long'); const L = cig.length; if (L && (cig[L - 1] & 15) === o && (cig[L - 1] >>> 4) + n <= 0xFFFFFFF) cig[L - 1] += n * 16; else cig.push(n * 16 + o); };
    const put = (at, b) => { if (at < 0 || at + b.length > rl) throw bad('a feature runs past the end of its read'); seq.set(b, at); };
    if (!(rec.flag & 4)) {
      const fn = E.FN.int(s);
      if (fn < 0 || fn > rl + 1000) throw bad('a record has an impossible number of features');
      let seqPos = 1, refPos = apos - 1, prev = 0;
      const id = rec.refID;
      const qFill = () => { if (!(cf & 1) && rl && qual[0] === 255) qual.fill(30); };
      for (let f = 0; f < fn; f++) {
        const fc = E.FC.byte(s), fp = E.FP.int(s) + prev; prev = fp;
        if (fp <= 0 || fp > rl + 1) throw bad('a feature lies outside its read');
        if (fp > seqPos) { for (let k = seqPos; k < fp; k++) seq[k - 1] = refAt(id, refPos + k - seqPos); op(0, fp - seqPos); refPos += fp - seqPos; seqPos = fp; }
        switch (fc) {
          case 83: { const b = E.SC.bytes(s); put(fp - 1, b); op(4, b.length); seqPos += b.length; break; }            // S soft clip
          case 88: { const code = E.BS.byte(s) & 3; seq[fp - 1] = ch.subst[CRAM_BASE_IDX[refAt(id, refPos)]][code]; op(0, 1); seqPos++; refPos++; break; } // X substitution
          case 68: { const n = E.DL.int(s); if (n < 0) throw bad('a negative deletion'); op(2, n); refPos += n; break; } // D
          case 73: { const b = E.IN.bytes(s); put(fp - 1, b); op(1, b.length); seqPos += b.length; break; }            // I insertion
          case 105: { if (fp > rl) throw bad('a feature runs past the end of its read'); seq[fp - 1] = E.BA.byte(s); op(1, 1); seqPos++; break; } // i one inserted base
          case 66: { if (fp > rl) throw bad('a feature runs past the end of its read'); seq[fp - 1] = E.BA.byte(s); qFill(); qual[fp - 1] = E.QS.byte(s); op(0, 1); seqPos++; refPos++; break; } // B base + quality
          case 81: { qFill(); if (fp - 1 < rl) qual[fp - 1] = E.QS.byte(s); else E.QS.byte(s); break; }          // Q quality
          case 72: { const n = E.HC.int(s); if (n < 0) throw bad('a negative hard clip'); op(5, n); break; }           // H
          case 80: { const n = E.PD.int(s); if (n < 0) throw bad('a negative padding'); op(6, n); break; }             // P
          case 78: { const n = E.RS.int(s); if (n < 0) throw bad('a negative skip'); op(3, n); refPos += n; break; }   // N
          case 98: { const b = E.BB.bytes(s); put(fp - 1, b); op(0, b.length); seqPos += b.length; refPos += b.length; break; } // b bases
          case 113: { const b = E.QQ.bytes(s); qFill(); if (fp - 1 + b.length > rl) throw bad('qualities run past the read'); qual.set(b, fp - 1); break; } // q qualities
          default: throw bad('an unknown read feature ' + String.fromCharCode(fc));
        }
      }
      if (rl >= seqPos) {
        const R = id === onlyId ? only : null;
        if (R) { // the common tail: bases straight from the reference
          for (let k = seqPos, j = refPos - R.start0; k <= rl; k++, j++) seq[k - 1] = j >= 0 && j < R.bases.length ? R.bases[j] : 78;
        } else for (let k = seqPos; k <= rl; k++) seq[k - 1] = refAt(id, refPos + k - seqPos);
        op(0, rl - seqPos + 1); refPos += rl - seqPos + 1;
      }
      if (refPos > 0x7FFFFFFF) throw bad('an alignment ends past the largest position');
      rec.end = Math.max(refPos, rec.pos + 1);
      rec.mapq = E.MQ.int(s);
      if (cf & 1) qual.set(E.QS.bytes(s, rl));
    } else {
      if (rl) seq.set(E.BA.bytes(s, rl));
      if (cf & 1) qual.set(E.QS.bytes(s, rl));
      rec.end = rec.pos + 1;
    }
    if (!ch.QO && (rec.flag & 16)) qual.reverse(); // qualities stored in sequencing order
    rec.cigar = Uint32Array.from(cig);
    const noSeq = !!(cf & 8);
    const ls = rec.lseq = noSeq ? 0 : rl;
    const packed = rec.seq = new Uint8Array((ls + 1) >> 1);
    for (let k = 0; k + 1 < ls; k += 2) packed[k >> 1] = (NT16[seq[k]] << 4) | NT16[seq[k + 1]];
    if (ls & 1) packed[ls >> 1] = NT16[seq[ls - 1]] << 4;
    rec.qual = noSeq ? new Uint8Array(0) : qual;
    if (!ch.RN && !rec.name) rec.name = file.src.name + ':' + (sh.counter + i);
    recs[i] = rec;
  }
  resolveMates(recs, mateLine, tlenOf, explicitTlen); CRAM_STATS.records += cramNow() - tRec;
  return recs;
}
function latin1Str(u8) { let s = ''; for (let i = 0; i < u8.length; i++) s += String.fromCharCode(u8[i]); return s; }

// Mates inside a slice: positions, flags and template length, as htslib computes them
// (cram_decode_slice_xref), so the result matches what samtools writes back to BAM.
function resolveMates(recs, mateLine, tlenOf, explicitTlen) {
  const N = recs.length;
  for (let i = 0; i < N; i++) if (mateLine[i] >= 0 && (mateLine[i] <= i || mateLine[i] >= N)) throw new CgvError('cramBad', { what: 'a record points to a mate outside its slice' });
  const ap = r => r.pos + 1, ae = r => (r.flag & 4) ? r.pos + 1 : r.end; // 1-based start, inclusive end
  for (let i = 0; i < N; i++) {
    const cr = recs[i];
    if (mateLine[i] >= 0) {
      if (Number.isNaN(tlenOf[i])) {
        let id2 = i, aleft = ap(cr), aright = ae(cr), leftCnt = 0, rightCnt = 0, ref = cr.refID, guard = 0;
        do {
          const x = recs[id2];
          if (aleft > ap(x)) { aleft = ap(x); leftCnt = 1; } else if (aleft === ap(x)) leftCnt++;
          if (aright < ae(x)) { aright = ae(x); rightCnt = 1; } else if (aright === ae(x)) rightCnt++;
          if (mateLine[id2] === -1) { mateLine[id2] = i; break; }
          id2 = mateLine[id2];
          if (recs[id2].refID !== ref) ref = -1;
          if (++guard > N) throw new CgvError('cramBad', { what: 'a loop in mate links' });
        } while (id2 !== i);
        if (ref !== -1) {
          let tlen = aright - aleft + 1;
          const x = recs[i];
          if (ap(x) === aleft && (ae(x) < aright || leftCnt <= 1)) { tlenOf[i] = tlen; tlen = -tlen; }
          else if (ap(x) === aleft && ae(x) === aright && leftCnt > 1 && rightCnt > 1) { if (x.flag & 64) { tlenOf[i] = tlen; tlen = -tlen; } else tlenOf[i] = -tlen; }
          else tlenOf[i] = -tlen;
          for (let j = mateLine[i], g = 0; j !== i; j = mateLine[j]) { tlenOf[j] = tlen; if (++g > N) break; }
        } else {
          tlenOf[i] = 0;
          for (let j = mateLine[i], g = 0; j !== i; j = mateLine[j]) { tlenOf[j] = 0; if (++g > N) break; }
        }
      }
      const m = recs[mateLine[i]];
      cr.matePos = m.pos; cr.mateRef = m.refID;
      cr.flag |= 1;
      if (m.flag & 4) { cr.flag |= 8; tlenOf[i] = 0; }
      if (cr.flag & 4) tlenOf[i] = 0;
      if (m.flag & 16) cr.flag |= 32;
    }
    if (Number.isNaN(tlenOf[i])) tlenOf[i] = 0;
  }
  for (let i = 0; i < N; i++) recs[i].tlen = Number.isNaN(explicitTlen[i]) ? tlenOf[i] : explicitTlen[i];
}

// Coverage counts from a parsed record (BamFile does the same from raw bytes in addCoverage).
function addCoverageRec(rec, cov) {
  const lo = cov.start, hi = cov.end;
  let r = rec.pos, q = 0, bi = -1, blk = null;
  const slot = pos => { const b = Math.floor(pos / COV_BLOCK); if (b !== bi) { bi = b; blk = cov.block(b, true); } return (pos - b * COV_BLOCK) * 6; };
  for (const c of rec.cigar) {
    const op = c & 15, n = c >>> 4;
    const k0 = Math.max(0, lo - r), k1 = Math.min(n, hi - r); // only the part inside the window
    if (op === 0 || op === 7 || op === 8) {
      for (let k = k0; k < k1; k++) { const qi = q + k; const code = rec.lseq ? (qi & 1 ? rec.seq[qi >> 1] & 15 : rec.seq[qi >> 1] >> 4) : 15; const o = slot(r + k); blk[o + NIB_IDX[code]]++; }
      r += n; q += n;
    } else if (op === 2) { for (let k = k0; k < k1; k++) { const o = slot(r + k); blk[o + 5]++; } r += n; }
    else if (op === 3) r += n;
    else if (op === 1 || op === 4) q += n;
  }
}

if (typeof module !== 'undefined') module.exports = { md5, CRAM_STATS, CRAM_CHECK, decodeSlice, parseCompressionHeader, readEncoding, BitReader, CramFile, CramReader, parseContainerHeader, parseBlock, decodeBlock, parseCrai, parseSamHeaderText, CRAM_LIMITS };
