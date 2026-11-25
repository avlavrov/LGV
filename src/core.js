// CGV core: data sources, BGZF, BAI, BAM, FASTA, model (coverage, layout).
// No external dependencies. Decompression uses the browser's built-in DecompressionStream('deflate-raw').
'use strict';

// ======================= Data sources =======================
// One interface: read(start, len) -> Uint8Array, readAll() -> Uint8Array.
class FileSource {
  constructor(blob, name) { this.blob = blob; this.name = name || blob.name || 'file'; this.bytes = 0; }
  async read(start, len) {
    const end = Math.min(start + len, this.blob.size);
    if (end <= start) return new Uint8Array(0);
    const b = new Uint8Array(await this.blob.slice(start, end).arrayBuffer());
    this.bytes += b.length; return b;
  }
  async readAll() { const b = new Uint8Array(await this.blob.arrayBuffer()); this.bytes += b.length; return b; }
}

class UrlSource {
  constructor(url, headers) {
    this.url = url; this.headers = headers || {}; this.bytes = 0;
    try { this.name = decodeURIComponent(new URL(url).pathname.split('/').pop()) || url; } catch (e) { this.name = url; }
  }
  async _fetch(range) {
    const h = Object.assign({}, this.headers);
    if (range) h.Range = range;
    let r;
    try { r = await fetch(this.url, { headers: h, credentials: 'omit', mode: 'cors' }); }
    catch (e) { throw new CgvError('network', { name: this.name }); }
    if (r.status === 401 || r.status === 403) throw new CgvError('auth', { status: r.status, name: this.name });
    if (r.status === 404) throw new CgvError('notfound', { name: this.name });
    if (!r.ok) throw new CgvError('http', { status: r.status, name: this.name });
    return r;
  }
  // Reads the body of a response, at most `max` bytes: a server that sends more than was asked for
  // (or ignores Range on a large file) cannot make CGV download or hold more than it needs.
  async _body(r, max) {
    if (!r.body) { const b = new Uint8Array(await r.arrayBuffer()); return b.length > max ? b.subarray(0, max) : b; }
    const rd = r.body.getReader(), parts = []; let n = 0;
    try {
      while (n < max) { const { value, done } = await rd.read(); if (done) break; parts.push(value); n += value.length; }
    } catch (e) { throw new CgvError('network', { name: this.name }); }
    if (n >= max) rd.cancel().catch(() => {});
    const out = new Uint8Array(Math.min(n, max)); let o = 0;
    for (const p of parts) { if (o >= out.length) break; const k = Math.min(p.length, out.length - o); out.set(k === p.length ? p : p.subarray(0, k), o); o += k; }
    this.bytes += out.length;
    return out;
  }
  async read(start, len) {
    const r = await this._fetch('bytes=' + start + '-' + (start + len - 1));
    if (r.status === 200) { // the server ignored Range and sends the whole file from the start
      if (start + len > URL_NO_RANGE_MAX) { r.body && r.body.cancel().catch(() => {}); throw new CgvError('noRange', { name: this.name }); }
      this.noRange = true;
      const buf = await this._body(r, start + len);
      return buf.subarray(Math.min(start, buf.length));
    }
    return this._body(r, len);
  }
  async readAll() { const r = await this._fetch(null); const b = await this._body(r, URL_READ_ALL_MAX + 1); if (b.length > URL_READ_ALL_MAX) throw new CgvError('tooBig', { name: this.name }); return b; }
}
// Without Range support only the start of a file can be read (each read downloads it from byte 0).
const URL_NO_RANGE_MAX = 256 << 20;
// Whole files read at once are indexes (.bai, .crai, .fai, .gzi): far smaller than this.
const URL_READ_ALL_MAX = 512 << 20;

// Errors carry a code and values; the UI turns them into text in the chosen language.
const CGV_ERR_EN = {
  network: 'Cannot reach {name}. Check the link and that the server allows CORS and Range requests.',
  auth: 'Access denied ({status}): {name}. The link has expired or a token is needed.',
  noRange: 'The server of {name} does not support Range requests, so only the start of the file can be read. Use a server or bucket with Range support',
  tooBig: '{name} is too large to read as a whole: is it really an index file?',
  notfound: 'File not found (404): {name}', noIndex: 'No index found for {name} (tried {tried})', http: 'The server answered {status} for {name}',
  bgzf: 'The file is damaged or not BGZF (offset {offset})', bsize: 'A BGZF block has no BSIZE field',
  bgzfCrc: 'A compressed block is damaged: its checksum or length does not match (offset {offset})', inflate: 'A compressed block is damaged and does not decompress',
  baiBad: 'The BAI index is damaged or cut short', bamRecord: 'A BAM record is damaged (offset {offset})', faiBad: 'The FASTA index (.fai) is damaged at line {line}',
  faGzipBig: 'This gzip FASTA unpacks to more than {size}. Compress it with bgzip and index it (samtools faidx), then load it with its .fai and .gzi',
  notBai: 'This is not a BAI index (.bai)', notBam: 'This is not a BAM file (no BAM signature)', header: 'Could not read the BAM header',
  headerCut: 'The BAM header is cut short — the file is incomplete', headerBig: 'The BAM header is too large',
  fai: 'Empty or invalid FASTA index (.fai)', gzi: 'Invalid .gzi index', needGzi: 'A compressed .fa.gz needs its .gzi index',
  faIrregular: 'Sequence {name} has lines of different lengths, so it cannot be indexed. Rewrap it to a fixed width, e.g. with seqkit seq -w 60.', faEmpty: 'No sequences found — is this a FASTA file?',
};
class CgvError extends Error {
  constructor(code, vars) {
    super((CGV_ERR_EN[code] || code).replace(/\{(\w+)\}/g, (m, k) => (vars && k in vars ? vars[k] : m)));
    this.code = code; this.vars = vars || {};
  }
}

// ======================= BGZF =======================
// CRC32 (IEEE, as in zlib and gzip), slicing by 8: fast enough to check every BGZF block.
const CRC_T = (() => {
  const t = new Uint32Array(256 * 8);
  for (let n = 0; n < 256; n++) { let c = n; for (let k = 0; k < 8; k++) c = c & 1 ? 0xEDB88320 ^ (c >>> 1) : c >>> 1; t[n] = c >>> 0; }
  for (let n = 0; n < 256; n++) for (let k = 1; k < 8; k++) t[k * 256 + n] = (t[(k - 1) * 256 + n] >>> 8) ^ t[t[(k - 1) * 256 + n] & 255];
  return t;
})();
function crc32(u8, start, end) {
  if (start === undefined) { start = 0; end = u8.length; }
  const t = CRC_T; let c = 0xFFFFFFFF, i = start;
  for (; i + 8 <= end; i += 8) {
    const a = (c ^ (u8[i] | u8[i + 1] << 8 | u8[i + 2] << 16 | u8[i + 3] << 24)) >>> 0;
    c = t[1792 + (a & 255)] ^ t[1536 + (a >>> 8 & 255)] ^ t[1280 + (a >>> 16 & 255)] ^ t[1024 + (a >>> 24)] ^
      t[768 + u8[i + 4]] ^ t[512 + u8[i + 5]] ^ t[256 + u8[i + 6]] ^ t[u8[i + 7]];
  }
  for (; i < end; i++) c = t[(c ^ u8[i]) & 255] ^ (c >>> 8);
  return (c ^ 0xFFFFFFFF) >>> 0;
}
// Built-in deflate (DecompressionStream). Reads the output in pieces and stops as soon as it grows
// past `max`, so a small damaged or hostile block cannot inflate into gigabytes ("zip bomb").
// Any failure is a CgvError.
async function inflateStream(u8, format, max, code, vars, bigCode) {
  const ds = new DecompressionStream(format);
  const w = ds.writable.getWriter();
  w.write(u8).catch(() => {}); w.close().catch(() => {});
  const r = ds.readable.getReader(), parts = []; let total = 0;
  try {
    for (;;) {
      const { value, done } = await r.read();
      if (done) break;
      total += value.length;
      if (total > max) { r.cancel().catch(() => {}); throw new CgvError(bigCode || code, vars); }
      parts.push(value);
    }
  } catch (e) { if (e instanceof CgvError) throw e; throw new CgvError(code, vars); }
  if (parts.length === 1) return parts[0];
  const out = new Uint8Array(total); let o = 0;
  for (const p of parts) { out.set(p, o); o += p.length; }
  return out;
}
async function inflateRaw(u8, max) { return inflateStream(u8, 'deflate-raw', max || 65536, 'inflate', { offset: '?' }); }

// BGZF blocks in buf (a slab read at file offset baseCoff). A block cut off by the end of the slab
// is left for the next read. Every field is checked: a damaged block gives a CgvError, never a
// negative length or a read outside the slab.
function parseBgzfBlocks(buf, baseCoff) {
  const blocks = []; let p = 0;
  while (p + 18 <= buf.length) {
    if (buf[p] !== 31 || buf[p + 1] !== 139 || buf[p + 2] !== 8 || !(buf[p + 3] & 4)) throw new CgvError('bgzf', { offset: baseCoff + p });
    const xlen = buf[p + 10] | (buf[p + 11] << 8);
    let q = p + 12, bsize = -1; const xend = p + 12 + xlen;
    if (xend > buf.length) break;
    while (q + 4 <= xend) {
      const slen = buf[q + 2] | (buf[q + 3] << 8);
      if (buf[q] === 66 && buf[q + 1] === 67 && slen === 2 && q + 6 <= xend) bsize = buf[q + 4] | (buf[q + 5] << 8);
      q += 4 + slen;
    }
    if (bsize < 0) throw new CgvError('bsize');
    const total = bsize + 1;
    // header (12 + xlen) + deflate data + CRC32 + ISIZE (8)
    if (total < xlen + 20) throw new CgvError('bgzf', { offset: baseCoff + p });
    if (p + total > buf.length) break;
    const t = p + total - 8;
    const crc = (buf[t] | buf[t + 1] << 8 | buf[t + 2] << 16 | buf[t + 3] << 24) >>> 0;
    const isize = (buf[t + 4] | buf[t + 5] << 8 | buf[t + 6] << 16 | buf[t + 7] << 24) >>> 0;
    if (isize > 65536) throw new CgvError('bgzf', { offset: baseCoff + p });
    blocks.push({ coff: baseCoff + p, cstart: p + 12 + xlen, clen: total - xlen - 20, crc, isize });
    p += total;
  }
  return blocks;
}
// Checksums are verified always; the tests switch them off to send damaged data deeper into the parsers.
const BGZF_CHECK = { crc: true };
// Inflates one block from parseBgzfBlocks and checks its length and CRC32 against the block trailer.
async function inflateBgzf(raw, b) {
  const out = await inflateStream(raw.subarray(b.cstart, b.cstart + b.clen), 'deflate-raw', b.isize, 'bgzfCrc', { offset: b.coff });
  if (out.length !== b.isize || (BGZF_CHECK.crc && crc32(out, 0, out.length) !== b.crc)) throw new CgvError('bgzfCrc', { offset: b.coff });
  return out;
}

// Reads and inflates the blocks from cStart through the block starting at cLast.
async function readBgzf(source, cStart, cLast) {
  const raw = await source.read(cStart, cLast - cStart + 65536);
  const blocks = parseBgzfBlocks(raw, cStart).filter(b => b.coff <= cLast);
  const parts = await Promise.all(blocks.map(b => inflateBgzf(raw, b)));
  let total = 0; const offs = new Map();
  blocks.forEach((b, i) => { offs.set(b.coff, total); total += parts[i].length; });
  const data = new Uint8Array(total); let o = 0;
  for (const p of parts) { data.set(p, o); o += p.length; }
  return { data, offs };
}

// ======================= BAI =======================
// BGZF virtual offset as a Number: coffset * 65536 + uoffset (exact for files up to ~137 GB).
function voffNum(dv, o) { const lo = dv.getUint32(o, true), hi = dv.getUint32(o + 4, true); return (hi * 65536 + (lo >>> 16)) * 65536 + (lo & 0xffff); }

function reg2bins(beg, end) {
  const list = [0]; end -= 1;
  for (let k = 1 + (beg >> 26); k <= 1 + (end >> 26); k++) list.push(k);
  for (let k = 9 + (beg >> 23); k <= 9 + (end >> 23); k++) list.push(k);
  for (let k = 73 + (beg >> 20); k <= 73 + (end >> 20); k++) list.push(k);
  for (let k = 585 + (beg >> 17); k <= 585 + (end >> 17); k++) list.push(k);
  for (let k = 4681 + (beg >> 14); k <= 4681 + (end >> 14); k++) list.push(k);
  return list;
}

class BaiIndex {
  constructor(u8) {
    if (u8.length < 8 || u8[0] !== 66 || u8[1] !== 65 || u8[2] !== 73 || u8[3] !== 1) throw new CgvError('notBai');
    this.dv = new DataView(u8.buffer, u8.byteOffset, u8.byteLength);
    const dv = this.dv, L = u8.length;
    // every count is checked against the bytes that are left before it is used, so a damaged index
    // gives a clear error instead of reading outside the file or allocating a huge table
    const count = (p, size) => { if (p + 4 > L) throw new CgvError('baiBad'); const n = dv.getInt32(p, true); if (n < 0 || p + 4 + n * size > L) throw new CgvError('baiBad'); return n; };
    const nref = count(4, 8);
    this.refs = []; let p = 8;
    for (let r = 0; r < nref; r++) {
      const start = p; const nbin = count(p, 8); p += 4;
      for (let b = 0; b < nbin; b++) { const nchunk = count(p + 4, 16); p += 8 + nchunk * 16; }
      const nintv = count(p, 8); p += 4 + nintv * 8;
      this.refs.push({ start, parsed: null });
    }
  }
  ref(r) {
    const R = this.refs[r]; if (!R) return null;
    if (R.parsed) return R.parsed;
    const dv = this.dv; let p = R.start;
    const nbin = dv.getInt32(p, true); p += 4; const bins = new Map();
    for (let b = 0; b < nbin; b++) {
      const bin = dv.getUint32(p, true), nchunk = dv.getInt32(p + 4, true); p += 8;
      if (bin === 37450) { p += nchunk * 16; continue; }
      const list = [];
      for (let c = 0; c < nchunk; c++) { list.push([voffNum(dv, p), voffNum(dv, p + 8)]); p += 16; }
      bins.set(bin, list);
    }
    const nintv = dv.getInt32(p, true); p += 4;
    const linear = new Float64Array(nintv);
    for (let i = 0; i < nintv; i++) linear[i] = voffNum(dv, p + i * 8);
    return (R.parsed = { bins, linear });
  }
  chunks(r, beg, end) {
    const R = this.ref(r); if (!R) return [];
    let minOff = 0;
    if (R.linear.length) minOff = R.linear[Math.min(beg >> 14, R.linear.length - 1)];
    const cs = [];
    for (const b of reg2bins(beg, end)) {
      const list = R.bins.get(b); if (!list) continue;
      for (const c of list) if (c[1] > minOff) cs.push([Math.max(c[0], minOff), c[1]]);
    }
    cs.sort((a, b) => a[0] - b[0]);
    const merged = [];
    for (const c of cs) {
      const last = merged[merged.length - 1];
      // merge when the gap is under 64 KB of compressed data: one request is cheaper than two
      if (last && Math.floor(c[0] / 65536) <= Math.floor(last[1] / 65536) + 65536) { if (c[1] > last[1]) last[1] = c[1]; }
      else merged.push([c[0], c[1]]);
    }
    return merged;
  }
}

// ======================= BAM =======================
function latin1(u8, start, len) { let s = ''; for (let i = 0; i < len; i++) s += String.fromCharCode(u8[start + i]); return s; }
const SEQ_CODES = '=ACMGRSVTWYHKDBN';
const CIGAR_OPS = 'MIDNSHP=X';

// Sizes taken from a BAM file are checked against these before anything is allocated or looped over.
const BAM_LIMITS = {
  header: 64 << 20,     // SAM header text
  refs: 10000000,       // reference sequences in the header
  record: 64 << 20,     // one alignment record (an 8 Mb nanopore read with qualities and tags fits)
  readLength: 1 << 23,  // bases in one read, as for CRAM
};
// Checks that the fields of the record at p (block_size at p) fit inside it and agree with each
// other, as htslib does when it reads a record. Returns false for a damaged record.
function bamRecordOk(buf, dv, p, bs) {
  const s = p + 4, lrn = buf[s + 8], ncig = dv.getUint16(s + 12, true), lseq = dv.getInt32(s + 16, true);
  if (lrn < 1 || lseq < 0 || lseq > BAM_LIMITS.readLength) return false;
  const cig = s + 32 + lrn;
  if (32 + lrn + 4 * ncig + ((lseq + 1) >> 1) + lseq > bs) return false;
  if (buf[cig - 1] !== 0) return false; // the read name ends with NUL
  let qlen = 0;
  for (let i = 0; i < ncig; i++) { const c = dv.getUint32(cig + 4 * i, true), op = c & 15; if (op > 8) return false; if (op === 0 || op === 1 || op === 4 || op === 7 || op === 8) qlen += c >>> 4; }
  // a stored sequence must be as long as the CIGAR says; without one, the CIGAR is still bounded
  return lseq ? qlen === lseq || ncig === 0 : qlen <= BAM_LIMITS.readLength;
}

// FNV-1a hash of a read name (ASCII), the key that downsampling keeps reads by (BAM and CRAM alike).
function nameHash(name) { let h = 2166136261; for (let i = 0; i < name.length; i++) { h ^= name.charCodeAt(i) & 255; h = Math.imul(h, 16777619); } return h >>> 0; }

// Downsampling, shared by BAM and CRAM: at most `max` reads starting in each `window` bp, chosen by
// the hash of the read name (the smallest hashes), so the choice is the same on every load. Pairs are
// kept together: a read whose mate starts in a later window of the same query decides for both, and
// its mate is kept or left out with it when it comes (kept mates may take a window up to 2 × max).
// Names in `keep` are always kept, with their mates. What is kept for a window does not depend on
// where the query began, except for reads whose mate starts before that (within an insert size).
class Downsampler {
  constructor(ds, keep, out) {
    this.ds = ds; this.keep = keep && keep.length ? new Set(keep) : null; this.out = out;
    this.win = -1; this.n = 0; this.lost = false; this.kept = []; this.kh = []; this.top = -1; this.extra = []; this.xh = []; this.pend = [];
    this.open = new Map(); // name hash → kept? for mates still to come
  }
  // pos: read start; h: name hash; ahead: its mate starts in a later window of this query;
  // name(), rec(): the read name and record, made only when needed
  add(pos, h, ahead, name, rec) {
    const w = Math.floor(pos / this.ds.window);
    if (w !== this.win) { this.flush(); this.win = w; }
    this.n++;
    if (ahead) this.pend.push(h);
    if (this.keep && this.keep.has(name())) { this.extra.push(rec()); this.xh.push(h); return; }
    if (this.open.has(h)) { // the second mate: as the first one went
      const k = this.open.get(h); this.open.delete(h);
      if (k && this.extra.length < this.ds.max) { this.extra.push(rec()); this.xh.push(h); } else this.lost = true; // a window holds at most 2 × max
      return;
    }
    if (this.kept.length < this.ds.max) { this.kept.push(rec()); this.kh.push(h); this.top = -1; return; }
    this.lost = true;
    if (this.top < 0) { this.top = 0; for (let i = 1; i < this.kh.length; i++) if (this.kh[i] > this.kh[this.top]) this.top = i; } // the largest hash kept
    if (h < this.kh[this.top]) { this.kept[this.top] = rec(); this.kh[this.top] = h; this.top = -1; }
  }
  flush() {
    const all = this.extra.length ? this.kept.concat(this.extra) : this.kept;
    if (all.length) { all.sort((a, b) => a.pos - b.pos); for (const r of all) this.out.reads.push(r); }
    if (this.lost) this.out.sampled.push([this.win * this.ds.window, (this.win + 1) * this.ds.window]);
    if (this.pend.length) {
      const yes = new Set(this.kh); for (const h of this.xh) yes.add(h);
      for (const h of this.pend) this.open.set(h, yes.has(h));
      while (this.open.size > 500000) this.open.delete(this.open.keys().next().value); // mates that never came
    }
    this.kept = []; this.kh = []; this.extra = []; this.xh = []; this.pend = []; this.n = 0; this.top = -1; this.lost = false;
  }
}
class BamFile {
  constructor(bamSource, baiSource) { this.src = bamSource; this.baiSrc = baiSource; this.maxSpan = 0; }
  // baiSource may be a list of candidates (e.g. sample.bam.bai, then sample.bai): the first one that
  // exists is used.
  async open() {
    const [idx] = await Promise.all([this._readIndex(), this._readHeader()]);
    this.index = new BaiIndex(idx);
    return this;
  }
  async _readIndex() {
    const list = Array.isArray(this.baiSrc) ? this.baiSrc : [this.baiSrc];
    let first = null;
    for (const s of list) {
      try { const b = await s.readAll(); this.baiSrc = s; return b; }
      catch (e) {
        // a missing object answers 404, or 403 on buckets without listing rights, or fails CORS
        if (list.length === 1 || !['notfound', 'auth', 'http', 'network'].includes(e.code)) throw e;
        if (!first || first.code === 'notfound') first = e;
      }
    }
    if (first.code === 'notfound' || first.code === 'auth') throw new CgvError('noIndex', { name: this.src.name, tried: list.map(s => s.name).join(', ') });
    throw first;
  }
  async _readHeader() {
    let size = 1 << 16;
    for (let attempt = 0; attempt < 12; attempt++) {
      const raw = await this.src.read(0, size);
      const blocks = parseBgzfBlocks(raw, 0);
      if (!blocks.length) throw new CgvError('header');
      const parts = await Promise.all(blocks.map(b => inflateBgzf(raw, b)));
      let total = 0; parts.forEach(p => total += p.length);
      const d = new Uint8Array(total); let o = 0; for (const p of parts) { d.set(p, o); o += p.length; }
      const h = tryParseBamHeader(d);
      if (h) { Object.assign(this, h); return; }
      if (raw.length < size) throw new CgvError('headerCut');
      size *= 2;
    }
    throw new CgvError('headerBig');
  }
  refId(name) {
    if (this.refIndex.has(name)) return this.refIndex.get(name);
    if (!this._keys) { this._keys = new Map(); this.refs.forEach((r, i) => { const k = contigKey(r.name); if (!this._keys.has(k)) this._keys.set(k, i); }); }
    const id = this._keys.get(contigKey(name));
    return id === undefined ? -1 : id;
  }
  // Reads overlapping [beg, end) on refId.
  // opts.filters: load-time filters (flag bits, MAPQ). opts.downsample: { max, window } — keep at most
  // `max` reads starting in each `window` bp, as IGV does; coverage still counts every read that passes
  // the filters. Which reads are kept is decided by a hash of the read name: the reads of a window
  // with the smallest hashes. Both mates of a pair have the same name, so where the two windows are
  // about as deep, a pair is kept or left out together (a random choice keeps the mate of a kept read
  // only by chance); the choice is also the same on every load. opts.keep: read names always kept
  // (a mate being opened). Returns { reads, cov (6 counters per bp over [beg, end)), total, sampled: [[s, e]...] }.
  async fetch(refId, beg, end, opts) {
    opts = opts || {};
    const flt = opts.filters || null;
    const ds = opts.downsample && opts.downsample.max > 0 && opts.downsample.window > 0 ? opts.downsample : null;
    const cov = new SparseCoverage(beg, end);
    // opts.depthOnly (wide windows): depth from CIGAR ranges, without counting bases one by one;
    // several times cheaper, and allele colors cannot be seen at that scale anyway
    const depth = opts.depthOnly ? new DepthDiff(beg, end) : null;
    const out = { reads: [], cov, total: 0, sampled: [] };
    const dsm = ds ? new Downsampler(ds, opts.keep, out) : null;
    const visit = (buf, dv, p, pos) => {
      const s = p + 4, flag = dv.getUint16(s + 14, true), mapq = buf[s + 9];
      if (flag & 4) return;
      const lrn = buf[s + 8], ncig = dv.getUint16(s + 12, true), cig = s + 32 + lrn;
      let refLen = 0;
      for (let i = 0; i < ncig; i++) { const c = dv.getUint32(cig + 4 * i, true), op = c & 15; if (op === 0 || op === 2 || op === 3 || op === 7 || op === 8) refLen += c >>> 4; }
      if (refLen > this.maxSpan) this.maxSpan = refLen;
      if (pos + Math.max(refLen, 1) <= beg) return;
      if (flt && !passFlags(flag, mapq, flt)) return;
      if (depth) depth.add(buf, dv, p); else addCoverage(buf, dv, p, cov);
      if (!ds || pos >= beg) out.total++;
      if (!ds) { out.reads.push(parseBamRecord(buf, dv, p)); return; }
      // Downsampling decides per window of read starts. A read starting before `beg` belongs to a window
      // this query only partly sees, so the choice there could differ from a wider query. Such reads
      // count for coverage but are not returned: what is returned for a window does not depend on where
      // the query began (see Downsampler for mates). Callers start the query a read length early.
      if (pos < beg) return;
      // the hash of the read name (without its closing NUL), as nameHash
      const ln = buf[s + 8], q0 = s + 32;
      let h = 2166136261;
      for (let i = 0; i < ln - 1; i++) { h ^= buf[q0 + i]; h = Math.imul(h, 16777619); }
      h >>>= 0;
      const nRef = dv.getInt32(s + 20, true), nPos = dv.getInt32(s + 24, true);
      const ahead = (flag & 9) === 1 && nRef === refId && nPos < end && Math.floor(nPos / ds.window) > Math.floor(pos / ds.window);
      dsm.add(pos, h, ahead, () => String.fromCharCode.apply(null, buf.subarray(q0, q0 + ln - 1)), () => parseBamRecord(buf, dv, p));
    };
    // chunks are disjoint and in file order, so reads arrive sorted by position
    for (const c of this.index.chunks(refId, beg, end)) await this._scanChunk(c, refId, beg, end, visit);
    if (dsm) dsm.flush();
    if (depth) depth.finish(cov);
    return out;
  }
  // Streams one index chunk in slabs, calling visit(buf, dv, p, pos) for records of refId that start
  // before `end`. Stops as soon as reads pass `end`, so the rest of the chunk is never read.
  // Remembers where records start in BGZF blocks, so later queries can skip ahead.
  async _scanChunk(c, refId, beg, end, visit) {
    let c0 = Math.floor(c[0] / 65536), u0 = c[0] % 65536;
    const cEnd = Math.floor(c[1] / 65536), uEnd = c[1] % 65536;
    const jump = this._jump(refId, beg, c0, u0, cEnd);
    if (jump) { c0 = jump.coff; u0 = jump.uoff; }
    if (!this.blocks) this.blocks = new Map();
    let carry = new Uint8Array(0), coff = c0, first = true;
    const SLAB = 1 << 20;
    for (;;) {
      const raw = await this.src.read(coff, SLAB);
      const blocks = raw.length ? parseBgzfBlocks(raw, coff) : [];
      if (!blocks.length) return;
      const use = blocks.filter(b => b.coff <= cEnd);
      const parts = await Promise.all(use.map(b => inflateBgzf(raw, b)));
      for (let i = 0; i < use.length; i++) {
        const bcoff = use[i].coff, shift = first ? u0 : 0;
        let part = shift ? parts[i].subarray(shift) : parts[i];
        first = false;
        let buf = part;
        if (carry.length) { buf = new Uint8Array(carry.length + part.length); buf.set(carry); buf.set(part, carry.length); }
        const dv = new DataView(buf.buffer, buf.byteOffset, buf.byteLength), cl = carry.length;
        let p = 0, cached = this.blocks.has(bcoff);
        while (p + 4 <= buf.length) {
          const bs = dv.getInt32(p, true);
          if (bs < 32 || bs > BAM_LIMITS.record) throw new CgvError('bamRecord', { offset: bcoff });
          if (p + 4 + bs > buf.length) break;
          const inBlock = p - cl; // >= 0: the record starts in this block, at uoff inBlock + shift
          if (bcoff === cEnd && inBlock >= 0 && inBlock + shift >= uEnd) return;
          const refID = dv.getInt32(p + 4, true), pos = dv.getInt32(p + 8, true);
          if (!cached && inBlock >= 0) { this.blocks.set(bcoff, { refID, pos, uoff: inBlock + shift }); cached = true; }
          if (refID !== refId) { if (refID > refId || refID < 0) return; p += 4 + bs; continue; }
          if (pos >= end) return;
          if (!bamRecordOk(buf, dv, p, bs)) throw new CgvError('bamRecord', { offset: bcoff });
          visit(buf, dv, p, pos);
          p += 4 + bs;
        }
        carry = buf.slice(p);
      }
      if (use.length < blocks.length) return;
      const last = blocks[blocks.length - 1];
      coff += last.cstart + last.clen + 8;
    }
  }
  // Compressed bytes the index points to for [beg, end): a cheap estimate of the work before reading.
  estimateBytes(refId, beg, end) {
    let n = 0;
    for (const c of this.index.chunks(refId, beg, end)) n += Math.floor(c[1] / 65536) - Math.floor(c[0] / 65536) + 65536;
    return n;
  }
  // Latest remembered record start in [c0, cEnd] that is safely left of `beg`
  // (further left than the longest alignment seen in this file).
  _jump(refId, beg, c0, u0, cEnd) {
    if (!this.blocks) return null;
    const limit = beg - Math.max(this.maxSpan, 1000) - 1;
    let best = null;
    for (const [coff, b] of this.blocks) {
      if (coff < c0 || coff > cEnd || b.refID !== refId || b.pos > limit) continue;
      if (coff === c0 && b.uoff <= u0) continue;
      if (!best || coff > best.coff) best = { coff, uoff: b.uoff };
    }
    return best;
  }
}

// Case-insensitive contig key: Chr1, chr1, CHR1 and 1 match; M and MT are the same contig.
function contigKey(name) {
  const k = String(name).toLowerCase().replace(/^chr/, '');
  return k === 'm' ? 'mt' : k;
}

// Name variants of one sequence: the whole name, the name without a version suffix (NG_012345.1),
// and each field of an NCBI-style header (gi|224589800|ref|NC_000001.10|).
function contigKeys(name) {
  const out = new Set(), add = n => { if (!n) return; const k = contigKey(n); out.add(k); out.add(k.replace(/\.\d+$/, '')); };
  add(name);
  if (name.includes('|')) for (const f of name.split('|')) if (f && !/^(gi|ref|gb|emb|dbj|lcl|\d+)$/i.test(f)) add(f);
  return [...out];
}
// Matches BAM contigs to reference sequences: by name (case, chr prefix, version and NCBI header
// fields ignored), then, for what is left, by a length that only one sequence on each side has.
// A name match with a different length is not a match but a conflict: the same name in another
// genome build (chr1 is 249,250,621 bp in GRCh37 and 248,956,422 bp in GRCh38). Coordinates would be
// wrong, so such contigs are never paired, and the caller refuses the sample when a large one differs.
function matchContigs(bamRefs, faRefs) {
  const byKey = new Map();
  for (const f of faRefs) for (const k of contigKeys(f.name)) if (!byKey.has(k)) byKey.set(k, f);
  const matched = [], used = new Set(), left = [], conflicts = [];
  for (const b of bamRefs) {
    let hit = null;
    for (const k of contigKeys(b.name)) if (byKey.has(k)) { hit = byKey.get(k); break; }
    if (hit && b.len && hit.len && b.len !== hit.len) { conflicts.push([b.name, b.len, hit.name, hit.len]); used.add(hit.name); continue; }
    if (hit) { matched.push([b.name, hit.name, 'name']); used.add(hit.name); } else left.push(b);
  }
  const count = (arr, len) => arr.filter(x => x.len === len).length;
  const faLeft = faRefs.filter(f => !used.has(f.name));
  for (const b of left) {
    if (count(left, b.len) !== 1) continue;
    const cand = faLeft.filter(f => f.len === b.len);
    if (cand.length === 1) { matched.push([b.name, cand[0].name, 'length']); used.add(cand[0].name); }
  }
  const m = new Set(matched.map(x => x[0]));
  return { matched, conflicts, unmatched: bamRefs.filter(b => !m.has(b.name)).map(b => b.name) };
}

function passFlags(flag, mapq, flt) {
  if (mapq < flt.minMapq) return false;
  if (flt.hideDuplicates && flag & 1024) return false;
  if (flt.hideSecondary && flag & 256) return false;
  if (flt.hideSupplementary && flag & 2048) return false;
  if (flt.hideQcFail && flag & 512) return false;
  return true;
}

// Coverage over [start, end) stored in 64 kb blocks that exist only where reads are, so a
// 50 Mb window with a few reads costs a few blocks, not gigabytes. 6 counters per base: A C G T N del.
const COV_BLOCK = 65536;
class SparseCoverage {
  constructor(start, end) { this.start = start; this.end = end; this.blocks = new Map(); }
  block(bi, create) {
    let b = this.blocks.get(bi);
    if (!b && create) { b = new Uint32Array(COV_BLOCK * 6); this.blocks.set(bi, b); }
    return b;
  }
  // counter k (0..5) at a position
  get(pos, k) { const b = this.blocks.get(Math.floor(pos / COV_BLOCK)); return b ? b[(pos % COV_BLOCK) * 6 + k] : 0; }
  // all six counters at a position into out (an array of 6), returns out
  counts(pos, out) {
    const b = pos >= this.start && pos < this.end ? this.blocks.get(Math.floor(pos / COV_BLOCK)) : null;
    if (!b) { out.fill(0); return out; }
    const o = (pos % COV_BLOCK) * 6; for (let k = 0; k < 6; k++) out[k] = b[o + k];
    return out;
  }
  // calls f(pos, block, offset) for every position in [a, b) that has a block; skips empty stretches fast
  forEach(a, b, f) {
    a = Math.max(a, this.start); b = Math.min(b, this.end);
    for (let bi = Math.floor(a / COV_BLOCK); bi * COV_BLOCK < b; bi++) {
      const blk = this.blocks.get(bi); if (!blk) continue;
      const s0 = Math.max(a, bi * COV_BLOCK), s1 = Math.min(b, (bi + 1) * COV_BLOCK);
      for (let pos = s0; pos < s1; pos++) f(pos, blk, (pos - bi * COV_BLOCK) * 6);
    }
  }
  hasData(a, b) { a = Math.max(a, this.start); b = Math.min(b, this.end); for (let bi = Math.floor(a / COV_BLOCK); bi * COV_BLOCK < b; bi++) if (this.blocks.has(bi)) return true; return false; }
}

// Adds one record to coverage straight from the BAM buffer, without building an object.
const NIB_IDX = new Uint8Array([4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4]);
function addCoverage(buf, dv, p, cov) {
  const s = p + 4, lrn = buf[s + 8], ncig = dv.getUint16(s + 12, true), lseq = dv.getInt32(s + 16, true);
  const cig = s + 32 + lrn, seq = cig + 4 * ncig, lo = cov.start, hi = cov.end;
  let r = dv.getInt32(s + 4, true), q = 0, bi = -1, blk = null;
  const slot = pos => { const b = Math.floor(pos / COV_BLOCK); if (b !== bi) { bi = b; blk = cov.block(b, true); } return (pos - b * COV_BLOCK) * 6; };
  for (let i = 0; i < ncig; i++) {
    const c = dv.getUint32(cig + 4 * i, true), op = c & 15, n = c >>> 4;
    // only the part inside [lo, hi) is visited, so a long operation costs nothing outside the window
    const k0 = Math.max(0, lo - r), k1 = Math.min(n, hi - r);
    if (op === 0 || op === 7 || op === 8) {
      for (let k = k0; k < k1; k++) { const qi = q + k; const code = lseq ? (qi & 1 ? buf[seq + (qi >> 1)] & 15 : buf[seq + (qi >> 1)] >> 4) : 15; const o = slot(r + k); blk[o + NIB_IDX[code]]++; }
      r += n; q += n;
    } else if (op === 2) { for (let k = k0; k < k1; k++) { const o = slot(r + k); blk[o + 5]++; } r += n; }
    else if (op === 3) r += n;
    else if (op === 1 || op === 4) q += n;
  }
}

// Depth as a difference array: +1 where an aligned block starts, -1 where it ends; summed at the end
// into the N counter of SparseCoverage, which is then marked depth-only (no allele counts).
class DepthDiff {
  constructor(start, end) { this.start = start; this.end = end; this.blocks = new Map(); }
  _blk(b) { let x = this.blocks.get(b); if (!x) { x = new Int32Array(COV_BLOCK); this.blocks.set(b, x); } return x; }
  _mark(pos, v) { if (pos < this.start) pos = this.start; if (pos >= this.end) return; const b = Math.floor(pos / COV_BLOCK); this._blk(b)[pos - b * COV_BLOCK] += v; }
  add(buf, dv, p) {
    const s = p + 4, lrn = buf[s + 8], ncig = dv.getUint16(s + 12, true), cig = s + 32 + lrn;
    let r = dv.getInt32(s + 4, true);
    for (let i = 0; i < ncig; i++) {
      const c = dv.getUint32(cig + 4 * i, true), op = c & 15, n = c >>> 4;
      if (op === 0 || op === 7 || op === 8) { if (r + n > this.start && r < this.end) { this._mark(r, 1); this._mark(r + n, -1); } r += n; }
      else if (op === 2 || op === 3) r += n;
    }
  }
  // the same from a parsed record (CRAM records arrive parsed)
  addRec(rec) {
    let r = rec.pos;
    for (const c of rec.cigar) {
      const op = c & 15, n = c >>> 4;
      if (op === 0 || op === 7 || op === 8) { if (r + n > this.start && r < this.end) { this._mark(r, 1); this._mark(r + n, -1); } r += n; }
      else if (op === 2 || op === 3) r += n;
    }
  }
  finish(cov) {
    let run = 0, prev = null;
    for (const b of [...this.blocks.keys()].sort((x, y) => x - y)) {
      if (prev !== null && b !== prev + 1) run = 0; // aligned blocks never span a whole empty 64 kb block
      const d = this.blocks.get(b), out = cov.block(b, true);
      for (let i = 0; i < COV_BLOCK; i++) { run += d[i]; if (run) out[i * 6 + 4] = run; }
      prev = b;
    }
    cov.depthOnly = true;
  }
}

// Reference bases kept only in stretches that are needed (under reads, or the view when zoomed in).
class SparseRef {
  constructor() { this.segs = []; this.last = null; } // segs: { s, a: Uint8Array } sorted by s, disjoint
  add(s, a) { this.segs.push({ s, a }); this.segs.sort((x, y) => x.s - y.s); }
  at(pos) {
    const L = this.last;
    if (L && pos >= L.s && pos < L.s + L.a.length) return L.a[pos - L.s];
    let lo = 0, hi = this.segs.length - 1;
    while (lo <= hi) { const m = (lo + hi) >> 1, g = this.segs[m]; if (pos < g.s) hi = m - 1; else if (pos >= g.s + g.a.length) lo = m + 1; else { this.last = g; return g.a[pos - g.s]; } }
    return 0;
  }
}
// Merges [s, e) intervals that are closer than `gap`.
function mergeIntervals(list, gap) {
  list.sort((a, b) => a[0] - b[0]);
  const out = [];
  for (const [s, e] of list) { const l = out[out.length - 1]; if (l && s <= l[1] + gap) { if (e > l[1]) l[1] = e; } else out.push([s, e]); }
  return out;
}
// Position of a read base at a reference position: { qi } for a base, { del: true } inside a
// deletion, null outside the read or in a skipped region.
function queryAt(rec, pos) {
  if (pos < rec.pos || pos >= rec.end) return null;
  let r = rec.pos, q = 0;
  for (const c of rec.cigar) {
    const op = c & 15, n = c >>> 4;
    if (op === 0 || op === 7 || op === 8) { if (pos < r + n) return { qi: q + pos - r }; r += n; q += n; }
    else if (op === 2) { if (pos < r + n) return { del: true }; r += n; }
    else if (op === 3) { if (pos < r + n) return null; r += n; }
    else if (op === 1 || op === 4) q += n;
  }
  return null;
}

function tryParseBamHeader(d) {
  if (d.length < 12) return null;
  if (d[0] !== 66 || d[1] !== 65 || d[2] !== 77 || d[3] !== 1) throw new CgvError('notBam');
  const dv = new DataView(d.buffer, d.byteOffset, d.byteLength);
  const lText = dv.getInt32(4, true);
  if (lText < 0) throw new CgvError('header');
  if (lText > BAM_LIMITS.header) throw new CgvError('headerBig');
  let p = 8 + lText;
  if (p + 4 > d.length) return null;
  const text = new TextDecoder().decode(d.subarray(8, 8 + lText)).replace(/\0+$/, '');
  const nRef = dv.getInt32(p, true); p += 4;
  if (nRef < 0 || nRef > BAM_LIMITS.refs) throw new CgvError('header');
  const refs = [], refIndex = new Map();
  for (let i = 0; i < nRef; i++) {
    if (p + 4 > d.length) return null;
    const ln = dv.getInt32(p, true); p += 4;
    if (ln < 1 || ln > 65536) throw new CgvError('header');
    if (p + ln + 4 > d.length) return null;
    const name = latin1(d, p, ln - 1); p += ln;
    const len = dv.getInt32(p, true); p += 4;
    if (len < 0) throw new CgvError('header');
    refs.push({ name, len }); refIndex.set(name, i);
  }
  const samples = [];
  for (const line of text.split('\n')) if (line.startsWith('@RG')) { const m = line.match(/\tSM:([^\t]+)/); if (m && !samples.includes(m[1])) samples.push(m[1]); }
  return { headerText: text, refs, refIndex, samples };
}

function parseBamRecord(data, dv, p) {
  const bs = dv.getInt32(p, true), s = p + 4;
  const lrn = data[s + 8], ncig = dv.getUint16(s + 12, true), lseq = dv.getInt32(s + 16, true);
  let q = s + 32;
  const name = latin1(data, q, lrn - 1); q += lrn;
  const cigar = new Uint32Array(ncig);
  for (let i = 0; i < ncig; i++) cigar[i] = dv.getUint32(q + 4 * i, true);
  q += 4 * ncig;
  const seqLen = (lseq + 1) >> 1;
  const seq = data.slice(q, q + seqLen); q += seqLen;
  const qual = data.slice(q, q + lseq); q += lseq;
  const tags = data.slice(q, s + bs);
  const pos = dv.getInt32(s + 4, true);
  let refLen = 0;
  for (let i = 0; i < ncig; i++) { const op = cigar[i] & 15; if (op === 0 || op === 2 || op === 3 || op === 7 || op === 8) refLen += cigar[i] >>> 4; }
  return {
    refID: dv.getInt32(s, true), pos, end: pos + Math.max(refLen, 1), name,
    mapq: data[s + 9], flag: dv.getUint16(s + 14, true), lseq,
    mateRef: dv.getInt32(s + 20, true), matePos: dv.getInt32(s + 24, true), tlen: dv.getInt32(s + 28, true),
    cigar, seq, qual, tagBytes: tags, _tags: null, feats: null,
  };
}

function readBase(rec, i) { const b = rec.seq[i >> 1]; return SEQ_CODES.charCodeAt(i & 1 ? b & 15 : b >> 4); }
function cigarString(rec) { let s = ''; for (const c of rec.cigar) s += (c >>> 4) + CIGAR_OPS[c & 15]; return s || '*'; }

function recTags(rec) {
  if (rec._tags) return rec._tags;
  const t = rec.tagBytes, out = {}, dv = new DataView(t.buffer, t.byteOffset, t.byteLength);
  let p = 0;
  const SZ = { A: 1, c: 1, C: 1, s: 2, S: 2, i: 4, I: 4, f: 4 };
  // a damaged tag area ends the list: what was read before it is kept
  while (p + 3 <= t.length) {
    const key = String.fromCharCode(t[p], t[p + 1]), type = String.fromCharCode(t[p + 2]); p += 3;
    let v;
    if (SZ[type] && p + SZ[type] > t.length) break;
    if (type === 'B' && p + 5 > t.length) break;
    switch (type) {
      case 'A': v = String.fromCharCode(t[p]); p += 1; break;
      case 'c': v = dv.getInt8(p); p += 1; break;
      case 'C': v = t[p]; p += 1; break;
      case 's': v = dv.getInt16(p, true); p += 2; break;
      case 'S': v = dv.getUint16(p, true); p += 2; break;
      case 'i': v = dv.getInt32(p, true); p += 4; break;
      case 'I': v = dv.getUint32(p, true); p += 4; break;
      case 'f': v = dv.getFloat32(p, true); p += 4; break;
      case 'Z': case 'H': { let e = p; while (e < t.length && t[e] !== 0) e++; v = latin1(t, p, e - p); p = e + 1; break; }
      case 'B': {
        const sub = String.fromCharCode(t[p]), n = dv.getInt32(p + 1, true); p += 5;
        const sz = { c: 1, C: 1, s: 2, S: 2, i: 4, I: 4, f: 4 }[sub]; if (!sz || n < 0 || p + n * sz > t.length) return (rec._tags = out);
        const arr = new Array(n);
        for (let i = 0; i < n; i++) {
          const o = p + i * sz;
          arr[i] = sub === 'c' ? dv.getInt8(o) : sub === 'C' ? t[o] : sub === 's' ? dv.getInt16(o, true) : sub === 'S' ? dv.getUint16(o, true) : sub === 'i' ? dv.getInt32(o, true) : sub === 'I' ? dv.getUint32(o, true) : dv.getFloat32(o, true);
        }
        p += n * sz; v = arr; break;
      }
      default: return (rec._tags = out);
    }
    out[key] = v;
  }
  return (rec._tags = out);
}

// Rebuilds reference bases from the MD tag into refArr (indexed from winStart).
function fillRefFromMD(rec, refArr, winStart) {
  const md = recTags(rec).MD; if (typeof md !== 'string') return;
  if (rec.end - rec.pos > 1 << 22) return; // per-base lists below: skip implausibly long alignments
  // aligned positions in order: refPos, and queryPos or -1 for deletions
  const rp = [], qp = []; let r = rec.pos, q = 0;
  for (const c of rec.cigar) {
    const op = c & 15, n = c >>> 4;
    if (op === 0 || op === 7 || op === 8) { for (let i = 0; i < n; i++) { rp.push(r + i); qp.push(q + i); } r += n; q += n; }
    else if (op === 2) { for (let i = 0; i < n; i++) { rp.push(r + i); qp.push(-1); } r += n; }
    else if (op === 3) r += n;
    else if (op === 1 || op === 4) q += n;
  }
  let k = 0, i = 0;
  const put = (pos, base) => { const j = pos - winStart; if (j >= 0 && j < refArr.length) refArr[j] = base; };
  while (i < md.length && k < rp.length) {
    const ch = md.charCodeAt(i);
    if (ch >= 48 && ch <= 57) {
      let n = 0; while (i < md.length && md.charCodeAt(i) >= 48 && md.charCodeAt(i) <= 57) { n = n * 10 + md.charCodeAt(i) - 48; i++; }
      for (let j = 0; j < n && k < rp.length; j++, k++) { if (qp[k] >= 0) { const b = readBase(rec, qp[k]); if (b !== 78) put(rp[k], b); } }
    } else if (ch === 94) { i++; while (i < md.length && md.charCodeAt(i) >= 65) { put(rp[k], md.charCodeAt(i) & ~32); i++; k++; } }
    else { put(rp[k], ch & ~32); i++; k++; }
  }
}

// Read features to draw: mismatches, insertions, deletions, skips, soft clips.
// t: 0=mismatch 1=ins 2=del 3=skip(N) 4=softclip-left 5=softclip-right
function computeFeatures(rec, refAt) {
  const f = []; let r = rec.pos, q = 0; const nc = rec.cigar.length;
  for (let ci = 0; ci < nc; ci++) {
    const c = rec.cigar[ci], op = c & 15, n = c >>> 4;
    if (op === 0 || op === 7 || op === 8) {
      // a record without a stored sequence (SEQ "*") has no bases to compare
      if (rec.lseq) for (let i = 0; i < n; i++) {
        const qb = readBase(rec, q + i);
        if (op === 8) { f.push({ t: 0, pos: r + i, b: qb, q: rec.qual[q + i] }); continue; }
        if (op === 7) continue;
        const rb = refAt(r + i);
        if (rb && qb !== rb && qb !== 78) f.push({ t: 0, pos: r + i, b: qb, q: rec.qual[q + i] });
      }
      r += n; q += n;
    } else if (op === 1) { f.push({ t: 1, pos: r, len: n, q0: q }); q += n; }
    else if (op === 2) { f.push({ t: 2, pos: r, len: n }); r += n; }
    else if (op === 3) { f.push({ t: 3, pos: r, len: n }); r += n; }
    else if (op === 4) { f.push({ t: ci === 0 || (ci === 1 && (rec.cigar[0] & 15) === 5) ? 4 : 5, pos: r, len: n, q0: q }); q += n; }
  }
  return f;
}

// Read base at a reference position (char code); 0 = not covered, 45 ('-') = deletion.
function baseAtRefPos(rec, pos) {
  if (pos < rec.pos || pos >= rec.end) return 0;
  let r = rec.pos, q = 0;
  for (const c of rec.cigar) {
    const op = c & 15, n = c >>> 4;
    if (op === 0 || op === 7 || op === 8) { if (pos < r + n) return readBase(rec, q + pos - r); r += n; q += n; }
    else if (op === 2 || op === 3) { if (pos < r + n) return op === 2 ? 45 : 0; r += n; }
    else if (op === 1 || op === 4) q += n;
  }
  return 0;
}

// ======================= FASTA =======================
// Plain FASTA + .fai, or bgzip-compressed FASTA + .fai + .gzi.
// ---- Building a FASTA index in the browser (no .fai given) ----
// One pass over the uncompressed bytes gives what samtools faidx writes:
// name, length, offset of the first base, bases per line, bytes per line.
class FaiBuilder {
  constructor() { this.entries = []; this.pos = 0; this.cur = null; this.inHeader = false; this.hdr = []; this.lb = 0; this.lbytes = 0; this.atStart = true; }
  _endLine() {
    const c = this.cur;
    if (!c) { this.lb = this.lbytes = 0; return; }
    if (this.lb === 0) { if (this.lbytes === 0) c.blank = true; this.lb = this.lbytes = 0; return; }
    if (c.blank || c.short) throw new CgvError('faIrregular', { name: c.name });
    if (c.lb === 0) { c.lb = this.lb; c.lw = this.lbytes + 1; }
    else if (this.lb !== c.lb || this.lbytes + 1 !== c.lw) { if (this.lb < c.lb) c.short = true; else throw new CgvError('faIrregular', { name: c.name }); }
    c.len += this.lb; this.lb = this.lbytes = 0;
  }
  push(u8) {
    for (let i = 0; i < u8.length; i++, this.pos++) {
      const ch = u8[i];
      if (this.inHeader) {
        if (ch === 10) {
          this.inHeader = false; this.atStart = true;
          let n = ''; for (const x of this.hdr) { if (x === 32 || x === 9 || x === 13) break; n += String.fromCharCode(x); }
          this.cur = { name: n, len: 0, off: this.pos + 1, lb: 0, lw: 0 }; this.entries.push(this.cur);
        } else if (this.hdr.length < 512) this.hdr.push(ch);
        continue;
      }
      if (ch === 10) { this._endLine(); this.atStart = true; continue; }
      if (this.atStart && ch === 62) { this.inHeader = true; this.hdr = []; continue; }
      this.atStart = false;
      this.lbytes++; if (ch !== 13) this.lb++;
    }
  }
  finish() {
    if (this.lbytes) { const c = this.cur; if (c && c.lb === 0) { c.lb = this.lb; c.lw = this.lbytes + 1; } if (c) c.len += this.lb; }
    const idx = new Map();
    for (const e of this.entries) if (!idx.has(e.name)) idx.set(e.name, { len: e.len, off: e.off, lb: e.lb || 1, lw: e.lw || 2 });
    if (!idx.size) throw new CgvError('faEmpty');
    return idx;
  }
}

// An ordinary (not bgzip) .fa.gz is unpacked into memory; this caps it (a whole human genome fits).
const FASTA_GZIP_MAX = 4 * 1024 * 1024 * 1024;
function isBgzfHeader(b) { return b.length >= 18 && b[0] === 31 && b[1] === 139 && (b[3] & 4) && b[12] === 66 && b[13] === 67; }

// Indexes a FASTA source with no .fai. Handles plain FASTA, BGZF (.fa.gz from bgzip,
// also builds the .gzi table) and ordinary gzip (decompressed into memory).
// Returns { idx, gzi, src } where src is the source to read sequence from.
async function indexFasta(src, size, onProgress) {
  const head = await src.read(0, 18);
  const b = new FaiBuilder(), SLAB = 8 << 20;
  const report = done => onProgress && onProgress(done, size);
  if (isBgzfHeader(head)) {
    const gzi = [[0, 0]]; let coff = 0, uoff = 0;
    for (;;) {
      const raw = await src.read(coff, SLAB);
      if (!raw.length) break;
      const blocks = parseBgzfBlocks(raw, coff);
      if (!blocks.length) break;
      const parts = await Promise.all(blocks.map(bl => inflateBgzf(raw, bl)));
      for (let i = 0; i < blocks.length; i++) {
        if (blocks[i].coff && parts[i].length) gzi.push([blocks[i].coff, uoff]);
        b.push(parts[i]); uoff += parts[i].length;
      }
      const last = blocks[blocks.length - 1];
      // parseBgzfBlocks gives cstart relative to raw; the block ends 8 bytes after its deflate data
      const consumed = last.cstart + last.clen + 8;
      coff += consumed; report(coff);
      if (raw.length < SLAB && consumed >= raw.length) break;
    }
    return { idx: b.finish(), gzi, src };
  }
  if (head[0] === 31 && head[1] === 139) { // ordinary gzip: no random access, keep the text in memory
    const all = await src.readAll();
    const ds = new DecompressionStream('gzip'); const w = ds.writable.getWriter(); w.write(all).catch(() => {}); w.close().catch(() => {});
    const r = ds.readable.getReader(), chunks = []; let done = 0;
    try {
      for (;;) {
        const { value, done: end } = await r.read(); if (end) break;
        done += value.length;
        if (done > FASTA_GZIP_MAX) { r.cancel().catch(() => {}); throw new CgvError('faGzipBig', { size: '4 GB' }); }
        b.push(value); chunks.push(value); report(done);
      }
    } catch (e) { if (e instanceof CgvError) throw e; throw new CgvError('inflate'); }
    return { idx: b.finish(), gzi: null, src: new FileSource(new Blob(chunks), src.name.replace(/\.gz$/i, '')) };
  }
  for (let off = 0; ; off += SLAB) {
    const raw = await src.read(off, SLAB);
    if (raw.length) b.push(raw);
    report(off + raw.length);
    if (raw.length < SLAB) break;
  }
  return { idx: b.finish(), gzi: null, src };
}

class IndexedFasta {
  constructor(faSource, faiSource, gziSource) { this.src = faSource; this.faiSrc = faiSource || null; this.gziSrc = gziSource || null; this.name = faSource.name; }
  // Index the FASTA in the browser when no .fai is available.
  async openUnindexed(size, onProgress) {
    const r = await indexFasta(this.src, size, onProgress);
    this.idx = r.idx; this.src = r.src; this.built = true;
    if (r.gzi) this.gzi = r.gzi;
    return this;
  }
  async open() {
    if (/\.gz$/i.test(this.src.name || '') && !this.gziSrc) throw new CgvError('needGzi');
    const t = new TextDecoder().decode(await this.faiSrc.readAll());
    this.idx = new Map();
    const lines = t.split('\n');
    for (let i = 0; i < lines.length; i++) {
      const line = lines[i].replace(/\r$/, ''); if (!line) continue;
      const f = line.split('\t'), int = x => /^\d{1,15}$/.test(x) ? +x : -1;
      const e = f.length >= 5 ? { len: int(f[1]), off: int(f[2]), lb: int(f[3]), lw: int(f[4]) } : null;
      // bases per line > 0, bytes per line = bases + 1 or 2 (\n or \r\n)
      if (!e || !f[0] || e.len < 0 || e.off < 0 || e.lb < 1 || e.lw < e.lb + 1 || e.lw > e.lb + 2) throw new CgvError('faiBad', { line: i + 1 });
      if (!this.idx.has(f[0])) this.idx.set(f[0], e);
    }
    if (!this.idx.size) throw new CgvError('fai');
    if (this.gziSrc) {
      const g = await this.gziSrc.readAll(), dv = new DataView(g.buffer, g.byteOffset, g.byteLength);
      if (g.length < 8) throw new CgvError('gzi');
      const n = dv.getUint32(0, true) + dv.getUint32(4, true) * 4294967296;
      if (g.length < 8 + n * 16) throw new CgvError('gzi');
      this.gzi = [[0, 0]];
      for (let i = 0; i < n; i++) {
        const o = 8 + i * 16, c = dv.getUint32(o, true) + dv.getUint32(o + 4, true) * 4294967296, u = dv.getUint32(o + 8, true) + dv.getUint32(o + 12, true) * 4294967296;
        // entries go up in both offsets, and a block holds at most 64 KB
        const last = this.gzi[this.gzi.length - 1];
        if (c <= last[0] || u <= last[1] || u - last[1] > 65536 || c > Number.MAX_SAFE_INTEGER || u > Number.MAX_SAFE_INTEGER) throw new CgvError('gzi');
        this.gzi.push([c, u]);
      }
    }
    return this;
  }
  entry(name) {
    if (this.idx.has(name)) return this.idx.get(name);
    if (!this._keys) { this._keys = new Map(); for (const [n, e] of this.idx) { const k = contigKey(n); if (!this._keys.has(k)) this._keys.set(k, e); } }
    return this._keys.get(contigKey(name)) || null;
  }
  has(name) { return !!this.entry(name); }
  // Uncompressed bytes [u0, u1) of the FASTA file.
  async _bytes(u0, u1) {
    if (!this.gzi) return this.src.read(u0, u1 - u0);
    const g = this.gzi; let lo = 0, hi = g.length - 1;
    const find = u => { let a = 0, b = g.length - 1, r = 0; while (a <= b) { const m = (a + b) >> 1; if (g[m][1] <= u) { r = m; a = m + 1; } else b = m - 1; } return r; };
    lo = find(u0); hi = find(Math.max(u0, u1 - 1));
    const { data } = await readBgzf(this.src, g[lo][0], g[hi][0]);
    return data.subarray(u0 - g[lo][1], u1 - g[lo][1]);
  }
  // Char codes (upper case) for [beg, end), or null if the contig is absent.
  async fetch(name, beg, end) {
    const e = this.entry(name); if (!e) return null;
    beg = Math.max(0, beg); end = Math.min(e.len, end); if (end <= beg) return new Uint8Array(0);
    const o0 = e.off + Math.floor(beg / e.lb) * e.lw + (beg % e.lb);
    const o1 = e.off + Math.floor(end / e.lb) * e.lw + (end % e.lb);
    const raw = await this._bytes(o0, o1);
    const out = new Uint8Array(end - beg); let j = 0;
    for (let i = 0; i < raw.length && j < out.length; i++) { const c = raw[i]; if (c === 10 || c === 13) continue; out[j++] = c & ~32; }
    return out;
  }
}

// ======================= Model =======================
const F_PAIRED = 1, F_PROPER = 2, F_REV = 16, F_MREV = 32, F_R1 = 64, F_R2 = 128, F_SECONDARY = 256, F_QCFAIL = 512, F_DUP = 1024, F_SUPP = 2048;

function passFilters(rec, flt) {
  if (rec.mapq < flt.minMapq) return false;
  if (flt.hideDuplicates && rec.flag & F_DUP) return false;
  if (flt.hideSecondary && rec.flag & F_SECONDARY) return false;
  if (flt.hideSupplementary && rec.flag & F_SUPP) return false;
  if (flt.hideQcFail && rec.flag & F_QCFAIL) return false;
  return true;
}

// Coverage: 6 counters per position (A C G T N del).
const BASE_IDX = new Int8Array(128).fill(4);
BASE_IDX[65] = 0; BASE_IDX[67] = 1; BASE_IDX[71] = 2; BASE_IDX[84] = 3;
function computeCoverage(reads, winStart, winLen) {
  const cov = new Uint32Array(winLen * 6);
  for (const rec of reads) {
    let r = rec.pos, q = 0;
    for (const c of rec.cigar) {
      const op = c & 15, n = c >>> 4;
      const i0 = Math.max(0, winStart - r), i1 = Math.min(n, winStart + winLen - r);
      if (op === 0 || op === 7 || op === 8) {
        for (let i = i0; i < i1; i++) cov[(r + i - winStart) * 6 + (rec.lseq ? BASE_IDX[readBase(rec, q + i)] : 4)]++;
        r += n; q += n;
      } else if (op === 2) { for (let i = i0; i < i1; i++) cov[(r + i - winStart) * 6 + 5]++; r += n; }
      else if (op === 3) r += n;
      else if (op === 1 || op === 4) q += n;
    }
  }
  return cov;
}

// Pair orientation for an FR library: 'FR' normal, 'RF', 'FF', 'RR', 'other' (mate on another contig), 'unpaired', 'mateUnmapped'
function pairOrientation(rec) {
  if (!(rec.flag & F_PAIRED)) return 'unpaired';
  if (rec.flag & 8) return 'mateUnmapped';
  if (rec.mateRef !== rec.refID) return 'other';
  const rev = !!(rec.flag & F_REV), mrev = !!(rec.flag & F_MREV);
  const leftRev = rec.pos <= rec.matePos ? rev : mrev, rightRev = rec.pos <= rec.matePos ? mrev : rev;
  if (!leftRev && rightRev) return 'FR';
  if (leftRev && !rightRev) return 'RF';
  return leftRev ? 'RR' : 'FF';
}

// Bounds of an unusual insert size from the loaded data (percentiles, as in IGV).
function insertSizeBounds(reads) {
  const v = [];
  for (const r of reads) if ((r.flag & F_PROPER) && r.tlen > 0) v.push(r.tlen);
  if (v.length < 20) return { lo: 50, hi: 1000 };
  v.sort((a, b) => a - b);
  return { lo: v[Math.floor(v.length * 0.005)], hi: v[Math.floor(v.length * 0.995)] };
}

// ======================= Read layout =======================
// Stable read layout. A read keeps its row for as long as the layout lives (panning, reloads,
// redraws, recoloring); only reads seen for the first time are placed, each into the lowest
// row where it fits. A new layout is made only when what decides placement changes
// (contig, grouping, sort, filters, downsampling), so the picture never reshuffles on its own.

// Identity of an alignment record: name + which segment (first/second, secondary, supplementary) + position.
function readKey(r) { return r.name + '\u0001' + (r.flag & 0x9C0) + '\u0001' + r.pos; }

// One row: disjoint intervals sorted by start.
class LayoutRow {
  constructor() { this.s = []; this.e = []; }
  // index of the first interval starting at or after x
  _at(x) { let lo = 0, hi = this.s.length; while (lo < hi) { const m = (lo + hi) >> 1; if (this.s[m] < x) lo = m + 1; else hi = m; } return lo; }
  fits(s, e, pad) {
    const i = this._at(s);
    if (i > 0 && this.e[i - 1] + pad > s) return false;
    if (i < this.s.length && this.s[i] < e + pad) return false;
    return true;
  }
  add(s, e) {
    const n = this.s.length;
    if (!n || this.s[n - 1] <= s) { this.s.push(s); this.e.push(e); return; }
    const i = this._at(s); this.s.splice(i, 0, s); this.e.splice(i, 0, e);
  }
  // does any interval overlap [a, b)?
  overlaps(a, b) { const i = this._at(b); return i > 0 && this.e[i - 1] > a; }
}

class ReadLayout {
  // sig: anything that, when changed, must produce a fresh layout. pad: minimum gap in bp.
  constructor(sig, pad) { this.sig = sig; this.pad = pad == null ? 2 : pad; this.place = new Map(); this.order = []; this.size = new Map(); }

  // reads: sorted by position. groupOf(r) -> group label. prio: { pos, key(r) -> [a, b] } places the
  // reads covering `pos` first, in key order, one per row from the top (sorting). Used only when the
  // layout is fresh; afterwards the rows it made are kept like any others.
  // Returns { groups: [{ key, rows: [[rec...]], start, count, sets }], totalRows }.
  // keyOf(item) gives a stable identity (readKey for reads; a pair key for pairs). Items need pos and end.
  update(reads, groupOf, prio, keyOf) {
    keyOf = keyOf || readKey;
    const pad = this.pad, fresh = this.place.size === 0;
    const byGroup = new Map();
    for (const r of reads) { const k = groupOf(r); let a = byGroup.get(k); if (!a) byGroup.set(k, a = []); a.push(r); }
    // group order is kept; a new group goes into its sorted place
    for (const k of byGroup.keys()) if (!this.order.includes(k)) this.order.push(k);
    this.order.sort((a, b) => String(a).localeCompare(String(b)));
    const place = new Map(), groups = [];
    let start = 0;
    for (const k of this.order) {
      const list = byGroup.get(k) || [];
      const sets = [], rowOf = new Array(list.length);
      const ensure = n => { while (sets.length <= n) sets.push(new LayoutRow()); };
      const pending = [];
      // 1. reads already placed keep their rows
      for (let i = 0; i < list.length; i++) {
        const r = list[i], key = keyOf(r), old = this.place.get(key);
        // an item keeps its row unless it grew (a pair whose mate arrived) and no longer fits there
        if (old && old.g === k) { ensure(old.r); if (sets[old.r].fits(r.pos, r.end, pad)) { sets[old.r].add(r.pos, r.end); rowOf[i] = old.r; place.set(key, old); continue; } }
        pending.push(i);
      }
      // 2. sorting: reads covering the sort position take rows from the top, in key order
      if (fresh && prio && pending.length) {
        const cover = pending.filter(i => list[i].pos <= prio.pos && list[i].end > prio.pos);
        const keyed = cover.map(i => [prio.key(list[i]), i]);
        keyed.sort((a, b) => a[0][0] - b[0][0] || a[0][1] - b[0][1] || list[a[1]].pos - list[b[1]].pos);
        const taken = new Set();
        keyed.forEach(([, i], n) => { ensure(n); sets[n].add(list[i].pos, list[i].end); rowOf[i] = n; taken.add(i); });
        for (let n = pending.length - 1; n >= 0; n--) if (taken.has(pending[n])) pending.splice(n, 1);
      }
      // 3. everything else goes into the lowest row where it fits
      if (sets.length === 0 && pending.length === list.length) {
        // empty group: sweep with a heap of rows in use and a heap of free row numbers (lowest first)
        const act = new MinHeap((a, b) => a[0] - b[0] || a[1] - b[1]), free = new MinHeap((a, b) => a - b);
        for (const i of pending) {
          const r = list[i];
          while (act.size && act.peek()[0] + pad <= r.pos) free.push(act.pop()[1]);
          const row = free.size ? free.pop() : sets.length;
          ensure(row); sets[row].add(r.pos, r.end); rowOf[i] = row; act.push([r.end, row]);
        }
      } else {
        for (const i of pending) {
          const r = list[i]; let row = 0;
          while (row < sets.length && !sets[row].fits(r.pos, r.end, pad)) row++;
          ensure(row); sets[row].add(r.pos, r.end); rowOf[i] = row;
        }
      }
      for (let i = 0; i < list.length; i++) { const key = keyOf(list[i]); if (!place.has(key)) place.set(key, { g: k, r: rowOf[i] }); }
      // rows never shrink while the layout lives, so groups below never jump up and down
      const n = Math.max(sets.length, this.size.get(k) || 0);
      this.size.set(k, n); ensure(n - 1);
      const rows = Array.from({ length: n }, () => []);
      for (let i = 0; i < list.length; i++) rows[rowOf[i]].push(list[i]);
      groups.push({ key: k, rows, start, count: list.length, sets });
      start += n;
    }
    this.place = place;
    return { groups, totalRows: start };
  }
}

// Reference length of a CIGAR string (for the MC tag: the mate's CIGAR).
function cigarStrRefLen(str) {
  let n = 0, num = 0;
  for (let i = 0; i < str.length; i++) {
    const c = str.charCodeAt(i);
    if (c >= 48 && c <= 57) { num = num * 10 + c - 48; continue; }
    if (c === 77 || c === 68 || c === 78 || c === 61 || c === 88) n += num; // M D N = X
    num = 0;
  }
  return n;
}
// Soft-clipped lengths at the left and right end of a record (hard clips are not shown).
function clipLens(r) {
  if (r._clips) return r._clips;
  const c = r.cigar, n = c.length; let L = 0, R = 0;
  for (let i = 0; i < n; i++) { const op = c[i] & 15; if (op === 5) continue; if (op === 4) L = c[i] >>> 4; break; }
  for (let i = n - 1; i >= 0; i--) { const op = c[i] & 15; if (op === 5) continue; if (op === 4) R = c[i] >>> 4; break; }
  return (r._clips = [L, R]);
}
// Hard-clipped lengths at both ends.
function hardClipLens(r) {
  const c = r.cigar, n = c.length;
  return [n && (c[0] & 15) === 5 ? c[0] >>> 4 : 0, n > 1 && (c[n - 1] & 15) === 5 ? c[n - 1] >>> 4 : 0];
}
// Drawn extent of a read: the aligned part, plus soft clips when they are shown.
function readExtent(r, clips) { if (!clips) return [r.pos, r.end]; const [L, R] = clipLens(r); return [r.pos - L, r.end + R]; }
// Soft clips of a CIGAR string (the MC tag), for the mate's drawn extent.
function cigarStrClips(str) {
  const m = str.match(/^(?:\d+H)?(\d+)S/), e = str.match(/(\d+)S(?:\d+H)?$/);
  const L = m ? +m[1] : 0, R = e && e.index > 0 ? +e[1] : 0; // an all-S CIGAR counts once, on the left
  return [L, R];
}

// Can this read be shown joined with its mate? Primary, paired, mate mapped on the same contig.
function pairable(r) { return (r.flag & 1) && !(r.flag & 4) && !(r.flag & 8) && !(r.flag & 0x900) && r.mateRef === r.refID; }
// Where the mate lies: exact with the MC tag, otherwise assumed as long as this read.
function mateExtent(r, clips) {
  const mc = recTags(r).MC, has = typeof mc === 'string';
  const len = has ? cigarStrRefLen(mc) : r.end - r.pos;
  let s = r.matePos, e = r.matePos + Math.max(1, len);
  if (clips && has) { const [L, R] = cigarStrClips(mc); s -= L; e += R; }
  return [s, e];
}
// Groups reads into layout units for pair view. A unit spans both mates (one taken from the mate
// fields when the mate itself is not loaded), so its size does not change when the mate arrives.
// Pairs spanning more than maxSpan, and reads that cannot pair, stay single.
// Returns units sorted by start: { key, pos, end, reads: [rec, rec?], mate: [start, end] | null }.
function pairUnits(reads, maxSpan, clips) {
  const open = new Map(), units = [];
  for (const r of reads) {
    const [rs, re] = readExtent(r, clips);
    if (pairable(r)) {
      const u = open.get(r.name);
      if (u) {
        u.reads.push(r); open.delete(r.name);
        const [s0, e0] = readExtent(u.reads[0], clips);
        u.pos = Math.min(s0, rs); u.end = Math.max(e0, re); u.mate = null;
        continue;
      }
      const m = mateExtent(r, clips), pos = Math.min(rs, m[0]), end = Math.max(re, m[1]);
      if (end - pos <= maxSpan) { const nu = { key: 'P\u0001' + r.name, pos, end, reads: [r], mate: m }; open.set(r.name, nu); units.push(nu); continue; }
    }
    units.push({ key: readKey(r), pos: rs, end: re, reads: [r], mate: null });
  }
  units.sort((a, b) => a.pos - b.pos || (a.key < b.key ? -1 : a.key > b.key ? 1 : 0));
  return units;
}

// Number of rows of a group that hold a read overlapping [a, b): what the visible part needs.
function rowsInUse(group, a, b) {
  for (let r = group.sets.length - 1; r >= 0; r--) if (group.sets[r].overlaps(a, b)) return r + 1;
  return 0;
}

class MinHeap {
  constructor(less) { this.a = []; this.cmp = less; }
  get size() { return this.a.length; }
  peek() { return this.a[0]; }
  push(x) { const a = this.a; a.push(x); let i = a.length - 1; while (i) { const p = (i - 1) >> 1; if (this.cmp(a[i], a[p]) < 0) { [a[i], a[p]] = [a[p], a[i]]; i = p; } else break; } }
  pop() {
    const a = this.a, top = a[0], last = a.pop();
    if (a.length) { a[0] = last; let i = 0; for (;;) { const l = 2 * i + 1, r = l + 1; let m = i; if (l < a.length && this.cmp(a[l], a[m]) < 0) m = l; if (r < a.length && this.cmp(a[r], a[m]) < 0) m = r; if (m === i) break; [a[i], a[m]] = [a[m], a[i]]; i = m; } }
    return top;
  }
}

// Plain packing of a list (no persistence), lowest free row first.
function packRows(reads, padBp) { const g = new ReadLayout(null, padBp).update(reads, () => '').groups[0]; return g ? g.rows : []; }


// ======================= Worker =======================
// The same file runs in a Web Worker (started from this script's own text): samples are read and
// inflated there, in parallel and off the page's thread, so the view stays responsive while loading.
// Messages: { id, op: 'open', key, bam, bai, cram } and { id, op: 'fetch', key, args } → { id, ok, res | err }.
// A CRAM in a worker asks the page for reference bases: { op: 'ref', key, rid, refId, a, b }, answered
// by { op: 'refResult', rid, bases } — the page owns the reference and the contig matching.
if (typeof WorkerGlobalScope !== 'undefined' && typeof self !== 'undefined' && self instanceof WorkerGlobalScope) {
  const files = new Map(), refWait = new Map();
  let refSeq = 0;
  const src = d => d.url ? new UrlSource(d.url, d.headers) : new FileSource(d.blob, d.name);
  self.onmessage = async e => {
    const { id, op, key } = e.data;
    if (op === 'refResult') { const cb = refWait.get(e.data.rid); if (cb) { refWait.delete(e.data.rid); cb(e.data.bases); } return; }
    try {
      if (op === 'open') {
        let f;
        if (e.data.cram) {
          f = await new CramFile(src(e.data.bam), src(e.data.bai)).open(); // CramFile: defined in cram.js, part of the same script
          f.refProvider = (refId, a, b) => new Promise(res => { const rid = ++refSeq; refWait.set(rid, res); self.postMessage({ op: 'ref', key, rid, refId, a, b }); });
        } else f = await new BamFile(src(e.data.bam), src(e.data.bai)).open();
        files.set(key, f); self.postMessage({ id, ok: true });
      }
      else if (op === 'close') files.delete(key);
      else if (op === 'fetch') {
        const bam = files.get(key), r = await bam.fetch(...e.data.args);
        const blocks = [...r.cov.blocks], cov = { start: r.cov.start, end: r.cov.end, depthOnly: !!r.cov.depthOnly, blocks };
        self.postMessage({ id, ok: true, res: { reads: r.reads, total: r.total, sampled: r.sampled, cov }, maxSpan: bam.maxSpan }, blocks.map(b => b[1].buffer));
      }
    } catch (err) { self.postMessage({ id, ok: false, err: { code: err.code || '', vars: err.vars || {}, message: String(err.message || err) } }); }
  };
}

if (typeof module !== 'undefined') module.exports = { nameHash, Downsampler, isBgzfHeader, BAM_LIMITS, bamRecordOk, addCoverage, CGV_ERR_EN, inflateRaw, inflateStream, inflateBgzf, parseBgzfBlocks, BGZF_CHECK, crc32, tryParseBamHeader, parseBamRecord, COV_BLOCK, NIB_IDX, DepthDiff, passFlags, clipLens, hardClipLens, readExtent, cigarStrClips, pairUnits, pairable, cigarStrRefLen, SparseCoverage, SparseRef, mergeIntervals, queryAt, ReadLayout, readKey, rowsInUse, matchContigs, contigKeys, contigKey, FaiBuilder, indexFasta, FileSource, UrlSource, CgvError, BamFile, BaiIndex, IndexedFasta, readBgzf, recTags, fillRefFromMD, computeFeatures, computeCoverage, baseAtRefPos, cigarString, packRows, passFilters, pairOrientation, insertSizeBounds, reg2bins };
