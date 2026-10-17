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
  async read(start, len) {
    const r = await this._fetch('bytes=' + start + '-' + (start + len - 1));
    const buf = new Uint8Array(await r.arrayBuffer());
    this.bytes += buf.length;
    if (r.status === 200) { // server ignored Range and sent the whole file
      this.noRange = true;
      return buf.subarray(start, Math.min(start + len, buf.length));
    }
    return buf;
  }
  async readAll() { const r = await this._fetch(null); const b = new Uint8Array(await r.arrayBuffer()); this.bytes += b.length; return b; }
}

// Errors carry a code and values; the UI turns them into text in the chosen language.
const CGV_ERR_EN = {
  network: 'Cannot reach {name}. Check the link and that the server allows CORS and Range requests.',
  auth: 'Access denied ({status}): {name}. The link has expired or a token is needed.',
  notfound: 'File not found (404): {name}', http: 'The server answered {status} for {name}',
  bgzf: 'The file is damaged or not BGZF (offset {offset})', bsize: 'A BGZF block has no BSIZE field',
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
async function inflateRaw(u8) {
  const ds = new DecompressionStream('deflate-raw');
  const w = ds.writable.getWriter();
  w.write(u8).catch(() => {}); w.close().catch(() => {});
  return new Uint8Array(await new Response(ds.readable).arrayBuffer());
}

function parseBgzfBlocks(buf, baseCoff) {
  const blocks = []; let p = 0;
  while (p + 18 <= buf.length) {
    if (buf[p] !== 31 || buf[p + 1] !== 139) throw new CgvError('bgzf', { offset: baseCoff + p });
    const xlen = buf[p + 10] | (buf[p + 11] << 8);
    let q = p + 12, bsize = -1; const xend = p + 12 + xlen;
    while (q + 4 <= xend) {
      const slen = buf[q + 2] | (buf[q + 3] << 8);
      if (buf[q] === 66 && buf[q + 1] === 67) bsize = buf[q + 4] | (buf[q + 5] << 8);
      q += 4 + slen;
    }
    if (bsize < 0) throw new CgvError('bsize');
    const total = bsize + 1;
    if (p + total > buf.length) break;
    blocks.push({ coff: baseCoff + p, cstart: p + 12 + xlen, clen: total - xlen - 20 });
    p += total;
  }
  return blocks;
}

// Reads and inflates the blocks from cStart through the block starting at cLast.
async function readBgzf(source, cStart, cLast) {
  const raw = await source.read(cStart, cLast - cStart + 65536);
  const blocks = parseBgzfBlocks(raw, cStart).filter(b => b.coff <= cLast);
  const parts = await Promise.all(blocks.map(b => inflateRaw(raw.subarray(b.cstart, b.cstart + b.clen))));
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
    const nref = this.dv.getInt32(4, true);
    this.refs = []; let p = 8; const dv = this.dv;
    for (let r = 0; r < nref; r++) {
      const start = p; const nbin = dv.getInt32(p, true); p += 4;
      for (let b = 0; b < nbin; b++) { const nchunk = dv.getInt32(p + 4, true); p += 8 + nchunk * 16; }
      const nintv = dv.getInt32(p, true); p += 4 + nintv * 8;
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

class BamFile {
  constructor(bamSource, baiSource) { this.src = bamSource; this.baiSrc = baiSource; this.maxSpan = 0; }
  async open() {
    const [idx] = await Promise.all([this.baiSrc.readAll(), this._readHeader()]);
    this.index = new BaiIndex(idx);
    return this;
  }
  async _readHeader() {
    let size = 1 << 16;
    for (let attempt = 0; attempt < 12; attempt++) {
      const raw = await this.src.read(0, size);
      const blocks = parseBgzfBlocks(raw, 0);
      if (!blocks.length) throw new CgvError('header');
      const parts = await Promise.all(blocks.map(b => inflateRaw(raw.subarray(b.cstart, b.cstart + b.clen))));
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
  // the filters. Returns { reads, cov (6 counters per bp over [beg, end)), total, sampled: [[s, e]...] }.
  async fetch(refId, beg, end, opts) {
    opts = opts || {};
    const flt = opts.filters || null;
    const ds = opts.downsample && opts.downsample.max > 0 && opts.downsample.window > 0 ? opts.downsample : null;
    const covLen = Math.max(0, end - beg), cov = new Uint32Array(covLen * 6);
    const out = { reads: [], cov, total: 0, sampled: [] };
    let win = -1, n = 0, kept = [], rng = 1;
    const flush = () => {
      if (kept.length) { kept.sort((a, b) => a.pos - b.pos); for (const r of kept) out.reads.push(r); }
      if (ds && n > ds.max) out.sampled.push([win * ds.window, (win + 1) * ds.window]);
      kept = []; n = 0;
    };
    const visit = (buf, dv, p, pos) => {
      const s = p + 4, flag = dv.getUint16(s + 14, true), mapq = buf[s + 9];
      if (flag & 4) return;
      const lrn = buf[s + 8], ncig = dv.getUint16(s + 12, true), cig = s + 32 + lrn;
      let refLen = 0;
      for (let i = 0; i < ncig; i++) { const c = dv.getUint32(cig + 4 * i, true), op = c & 15; if (op === 0 || op === 2 || op === 3 || op === 7 || op === 8) refLen += c >>> 4; }
      if (refLen > this.maxSpan) this.maxSpan = refLen;
      if (pos + Math.max(refLen, 1) <= beg) return;
      if (flt && !passFlags(flag, mapq, flt)) return;
      out.total++;
      addCoverage(buf, dv, p, cov, beg, covLen);
      if (!ds) { out.reads.push(parseBamRecord(buf, dv, p)); return; }
      const w = Math.floor(pos / ds.window);
      if (w !== win) { flush(); win = w; rng = (Math.imul(w + 1, 2654435761) >>> 0) || 1; }
      n++;
      if (kept.length < ds.max) kept.push(parseBamRecord(buf, dv, p));
      else { rng ^= rng << 13; rng >>>= 0; rng ^= rng >>> 17; rng ^= rng << 5; rng >>>= 0; const j = rng % n; if (j < ds.max) kept[j] = parseBamRecord(buf, dv, p); }
    };
    // chunks are disjoint and in file order, so reads arrive sorted by position
    for (const c of this.index.chunks(refId, beg, end)) await this._scanChunk(c, refId, beg, end, visit);
    flush();
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
      const parts = await Promise.all(use.map(b => inflateRaw(raw.subarray(b.cstart, b.cstart + b.clen))));
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
          if (bs < 32 || p + 4 + bs > buf.length) break;
          const inBlock = p - cl; // >= 0: the record starts in this block, at uoff inBlock + shift
          if (bcoff === cEnd && inBlock >= 0 && inBlock + shift >= uEnd) return;
          const refID = dv.getInt32(p + 4, true), pos = dv.getInt32(p + 8, true);
          if (!cached && inBlock >= 0) { this.blocks.set(bcoff, { refID, pos, uoff: inBlock + shift }); cached = true; }
          if (refID !== refId) { if (refID > refId || refID < 0) return; p += 4 + bs; continue; }
          if (pos >= end) return;
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
function matchContigs(bamRefs, faRefs) {
  const byKey = new Map();
  for (const f of faRefs) for (const k of contigKeys(f.name)) if (!byKey.has(k)) byKey.set(k, f.name);
  const matched = [], used = new Set(), left = [];
  for (const b of bamRefs) {
    let hit = null;
    for (const k of contigKeys(b.name)) if (byKey.has(k)) { hit = byKey.get(k); break; }
    if (hit) { matched.push([b.name, hit, 'name']); used.add(hit); } else left.push(b);
  }
  const count = (arr, len) => arr.filter(x => x.len === len).length;
  const faLeft = faRefs.filter(f => !used.has(f.name));
  for (const b of left) {
    if (count(left, b.len) !== 1) continue;
    const cand = faLeft.filter(f => f.len === b.len);
    if (cand.length === 1) { matched.push([b.name, cand[0].name, 'length']); used.add(cand[0].name); }
  }
  const m = new Set(matched.map(x => x[0]));
  return { matched, unmatched: bamRefs.filter(b => !m.has(b.name)).map(b => b.name) };
}

function passFlags(flag, mapq, flt) {
  if (mapq < flt.minMapq) return false;
  if (flt.hideDuplicates && flag & 1024) return false;
  if (flt.hideSecondary && flag & 256) return false;
  if (flt.hideSupplementary && flag & 2048) return false;
  if (flt.hideQcFail && flag & 512) return false;
  return true;
}

// Adds one record to coverage straight from the BAM buffer, without building an object.
const NIB_IDX = new Uint8Array([4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4]);
function addCoverage(buf, dv, p, cov, cs, clen) {
  const s = p + 4, lrn = buf[s + 8], ncig = dv.getUint16(s + 12, true), lseq = dv.getInt32(s + 16, true);
  const cig = s + 32 + lrn, seq = cig + 4 * ncig;
  let r = dv.getInt32(s + 4, true), q = 0;
  for (let i = 0; i < ncig; i++) {
    const c = dv.getUint32(cig + 4 * i, true), op = c & 15, n = c >>> 4;
    if (op === 0 || op === 7 || op === 8) {
      for (let k = 0; k < n; k++) {
        const j = r + k - cs;
        if (j >= 0 && j < clen) { const qi = q + k; const code = lseq ? (qi & 1 ? buf[seq + (qi >> 1)] & 15 : buf[seq + (qi >> 1)] >> 4) : 15; cov[j * 6 + NIB_IDX[code]]++; }
      }
      r += n; q += n;
    } else if (op === 2) { for (let k = 0; k < n; k++) { const j = r + k - cs; if (j >= 0 && j < clen) cov[j * 6 + 5]++; } r += n; }
    else if (op === 3) r += n;
    else if (op === 1 || op === 4) q += n;
  }
}

function tryParseBamHeader(d) {
  if (d.length < 12) return null;
  if (d[0] !== 66 || d[1] !== 65 || d[2] !== 77 || d[3] !== 1) throw new CgvError('notBam');
  const dv = new DataView(d.buffer, d.byteOffset, d.byteLength);
  const lText = dv.getInt32(4, true); let p = 8 + lText;
  if (p + 4 > d.length) return null;
  const text = new TextDecoder().decode(d.subarray(8, 8 + lText)).replace(/\0+$/, '');
  const nRef = dv.getInt32(p, true); p += 4;
  const refs = [], refIndex = new Map();
  for (let i = 0; i < nRef; i++) {
    if (p + 4 > d.length) return null;
    const ln = dv.getInt32(p, true); p += 4;
    if (p + ln + 4 > d.length) return null;
    const name = latin1(d, p, ln - 1); p += ln;
    const len = dv.getInt32(p, true); p += 4;
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
  while (p + 3 <= t.length) {
    const key = String.fromCharCode(t[p], t[p + 1]), type = String.fromCharCode(t[p + 2]); p += 3;
    let v;
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
        const sz = { c: 1, C: 1, s: 2, S: 2, i: 4, I: 4, f: 4 }[sub]; if (!sz) return (rec._tags = out);
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
      for (let i = 0; i < n; i++) {
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
      const parts = await Promise.all(blocks.map(bl => inflateRaw(raw.subarray(bl.cstart, bl.cstart + bl.clen))));
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
    for (;;) { const { value, done: end } = await r.read(); if (end) break; b.push(value); chunks.push(value); done += value.length; report(done); }
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
    for (const line of t.split('\n')) {
      const f = line.split('\t'); if (f.length < 5) continue;
      this.idx.set(f[0], { len: +f[1], off: +f[2], lb: +f[3], lw: +f[4] });
    }
    if (!this.idx.size) throw new CgvError('fai');
    if (this.gziSrc) {
      const g = await this.gziSrc.readAll(), dv = new DataView(g.buffer, g.byteOffset, g.byteLength);
      if (g.length < 8) throw new CgvError('gzi');
      const n = dv.getUint32(0, true) + dv.getUint32(4, true) * 4294967296;
      if (g.length < 8 + n * 16) throw new CgvError('gzi');
      this.gzi = [[0, 0]];
      for (let i = 0; i < n; i++) { const o = 8 + i * 16; this.gzi.push([dv.getUint32(o, true) + dv.getUint32(o + 4, true) * 4294967296, dv.getUint32(o + 8, true) + dv.getUint32(o + 12, true) * 4294967296]); }
    }
    return this;
  }
  entry(name) {
    if (this.alias && this.alias.has(name)) return this.alias.get(name);
    if (this.idx.has(name)) return this.idx.get(name);
    if (!this._keys) { this._keys = new Map(); for (const [n, e] of this.idx) { const k = contigKey(n); if (!this._keys.has(k)) this._keys.set(k, e); } }
    return this._keys.get(contigKey(name)) || null;
  }
  has(name) { return !!this.entry(name); }
  // Pairs the BAM's contigs with this reference's sequences, so lookups by BAM name work.
  // Returns { matched: [[bamName, fastaName, how]], unmatched: [bamName] }.
  linkTo(refs) {
    const r = matchContigs(refs, [...this.idx].map(([name, e]) => ({ name, len: e.len })));
    this.alias = new Map(r.matched.map(([b, f]) => [b, this.idx.get(f)]));
    return r;
  }
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
      if (op === 0 || op === 7 || op === 8) {
        for (let i = 0; i < n; i++) { const j = r + i - winStart; if (j >= 0 && j < winLen) cov[j * 6 + BASE_IDX[readBase(rec, q + i)]]++; }
        r += n; q += n;
      } else if (op === 2) { for (let i = 0; i < n; i++) { const j = r + i - winStart; if (j >= 0 && j < winLen) cov[j * 6 + 5]++; } r += n; }
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

// Dense row packing with a min-heap on row end.
function packRows(reads, padBp) {
  const rows = []; const heap = []; // [end, row]
  const less = (a, b) => a[0] < b[0] || (a[0] === b[0] && a[1] < b[1]);
  const push = x => { heap.push(x); let i = heap.length - 1; while (i) { const p = (i - 1) >> 1; if (less(heap[i], heap[p])) { [heap[i], heap[p]] = [heap[p], heap[i]]; i = p; } else break; } };
  const pop = () => { const top = heap[0], last = heap.pop(); if (heap.length) { heap[0] = last; let i = 0; for (;;) { const l = 2 * i + 1, r = l + 1; let m = i; if (l < heap.length && less(heap[l], heap[m])) m = l; if (r < heap.length && less(heap[r], heap[m])) m = r; if (m === i) break; [heap[i], heap[m]] = [heap[m], heap[i]]; i = m; } } return top; };
  for (const rec of reads) {
    if (heap.length && heap[0][0] + padBp <= rec.pos) { const [, row] = pop(); rows[row].push(rec); push([rec.end, row]); }
    else { rows.push([rec]); push([rec.end, rows.length - 1]); }
  }
  return rows;
}

if (typeof module !== 'undefined') module.exports = { matchContigs, contigKeys, contigKey, FaiBuilder, indexFasta, FileSource, UrlSource, CgvError, BamFile, BaiIndex, IndexedFasta, readBgzf, recTags, fillRefFromMD, computeFeatures, computeCoverage, baseAtRefPos, cigarString, packRows, passFilters, pairOrientation, insertSizeBounds, reg2bins };
