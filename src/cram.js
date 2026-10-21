// CGV: CRAM reading (CRAM 3.0 and 3.1). Written from the CRAM specification (samtools/hts-specs,
// CRAMv3.tex); no third-party code in this file.
//
// Stage 1 of CRAM support: file definition, containers, blocks (raw and gzip), the SAM header and
// the .crai index. Records are decoded in the next stage.
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

const CRAM_LIMITS = {
  block: 256 << 20,       // decoded size of one block
  container: 1024 << 20,  // compressed size of one container
  blocks: 100000,         // blocks in one container
  header: 64 << 20,       // SAM header text
  index: 512 << 20,       // decoded .crai
};
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
  cramLater: 'CRAM {v}: header and index read ({slices} slices). Showing reads from CRAM comes in the next step.',
});

// ---------- CRC32 (IEEE, as in zlib) ----------
const CRC_TABLE = (() => { const t = new Uint32Array(256); for (let n = 0; n < 256; n++) { let c = n; for (let k = 0; k < 8; k++) c = c & 1 ? 0xEDB88320 ^ (c >>> 1) : c >>> 1; t[n] = c >>> 0; } return t; })();
function crc32(u8, start, end) {
  let c = 0xFFFFFFFF;
  for (let i = start; i < end; i++) c = CRC_TABLE[(c ^ u8[i]) & 0xFF] ^ (c >>> 8);
  return (c ^ 0xFFFFFFFF) >>> 0;
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
    if (major >= 3) { const end = r.pos; const crc = r.u32(); if (crc32(u8, off, end) !== crc) throw new CgvError('cramCrc', { what: 'a container header' }); }
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
  if (major >= 3) { const e = r.pos; const crc = r.u32(); if (crc32(u8, off, e) !== crc) throw new CgvError('cramCrc', { what: 'a block' }); }
  b.next = r.pos;
  return b;
}
async function inflateGzip(u8) {
  const ds = new DecompressionStream('gzip');
  const w = ds.writable.getWriter();
  w.write(u8).catch(() => {}); w.close().catch(() => {});
  try { return new Uint8Array(await new Response(ds.readable).arrayBuffer()); }
  catch (e) { throw new CgvError('cramBad', { what: 'a gzip block does not inflate' }); }
}
// Decoded contents of a block. Codecs other than raw and gzip come in later stages.
async function decodeBlock(b) {
  let out;
  if (b.method === 0) out = b.data;
  else if (b.method === 1) out = await inflateGzip(b.data);
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
  if (raw.length >= 2 && raw[0] === 0x1F && raw[1] === 0x8B) u8 = await inflateGzip(raw);
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
    if (data.length < h.length) data = await this.src.read(dataOffset, h.length);
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
  async fetch() { throw new CgvError('cramLater', { v: this.version, slices: this.index.slices }); }
}

if (typeof module !== 'undefined') module.exports = { CramFile, CramReader, parseContainerHeader, parseBlock, decodeBlock, parseCrai, parseSamHeaderText, crc32, CRAM_LIMITS };
