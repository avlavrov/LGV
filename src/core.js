// CGV core: источники данных, BGZF, BAI, BAM, FASTA, модель (покрытие, раскладка).
// Без внешних зависимостей. Распаковка — встроенный DecompressionStream('deflate-raw').
'use strict';

// ======================= Источники данных =======================
// Единый интерфейс: read(start, len) -> Uint8Array, readAll() -> Uint8Array.
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
    catch (e) { throw new CgvError('Нет доступа к ' + this.name + '. Проверьте ссылку и что сервер разрешает CORS и Range-запросы.', 'network'); }
    if (r.status === 401 || r.status === 403) throw new CgvError('Доступ запрещён (' + r.status + '): ' + this.name + '. Ссылка устарела или нужен токен.', 'auth');
    if (r.status === 404) throw new CgvError('Файл не найден (404): ' + this.name, 'notfound');
    if (!r.ok) throw new CgvError('Сервер ответил ' + r.status + ' на ' + this.name, 'http');
    return r;
  }
  async read(start, len) {
    const r = await this._fetch('bytes=' + start + '-' + (start + len - 1));
    const buf = new Uint8Array(await r.arrayBuffer());
    this.bytes += buf.length;
    if (r.status === 200) { // сервер проигнорировал Range — вернул весь файл
      this.noRange = true;
      return buf.subarray(start, Math.min(start + len, buf.length));
    }
    return buf;
  }
  async readAll() { const r = await this._fetch(null); const b = new Uint8Array(await r.arrayBuffer()); this.bytes += b.length; return b; }
}

class CgvError extends Error { constructor(msg, kind) { super(msg); this.kind = kind; } }

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
    if (buf[p] !== 31 || buf[p + 1] !== 139) throw new CgvError('Файл повреждён или это не BGZF (смещение ' + (baseCoff + p) + ')', 'format');
    const xlen = buf[p + 10] | (buf[p + 11] << 8);
    let q = p + 12, bsize = -1; const xend = p + 12 + xlen;
    while (q + 4 <= xend) {
      const slen = buf[q + 2] | (buf[q + 3] << 8);
      if (buf[q] === 66 && buf[q + 1] === 67) bsize = buf[q + 4] | (buf[q + 5] << 8);
      q += 4 + slen;
    }
    if (bsize < 0) throw new CgvError('В блоке BGZF нет поля BSIZE', 'format');
    const total = bsize + 1;
    if (p + total > buf.length) break;
    blocks.push({ coff: baseCoff + p, cstart: p + 12 + xlen, clen: total - xlen - 20 });
    p += total;
  }
  return blocks;
}

// Читает и распаковывает блоки с cStart до блока cLast включительно.
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
// Виртуальное смещение BGZF как число: coffset * 65536 + uoffset (точно до ~137 ГБ файла).
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
    if (u8.length < 8 || u8[0] !== 66 || u8[1] !== 65 || u8[2] !== 73 || u8[3] !== 1) throw new CgvError('Это не индекс BAI (.bai)', 'format');
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
      // объединяем, если разрыв меньше 64 КБ сжатых данных: один запрос дешевле двух
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
  constructor(bamSource, baiSource) { this.src = bamSource; this.baiSrc = baiSource; }
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
      if (!blocks.length) throw new CgvError('Не удалось прочитать заголовок BAM', 'format');
      const parts = await Promise.all(blocks.map(b => inflateRaw(raw.subarray(b.cstart, b.cstart + b.clen))));
      let total = 0; parts.forEach(p => total += p.length);
      const d = new Uint8Array(total); let o = 0; for (const p of parts) { d.set(p, o); o += p.length; }
      const h = tryParseBamHeader(d);
      if (h) { Object.assign(this, h); return; }
      if (raw.length < size) throw new CgvError('Заголовок BAM обрывается — файл неполный', 'format');
      size *= 2;
    }
    throw new CgvError('Слишком большой заголовок BAM', 'format');
  }
  refId(name) {
    if (this.refIndex.has(name)) return this.refIndex.get(name);
    const alt = name.startsWith('chr') ? name.slice(3) : 'chr' + name;
    if (this.refIndex.has(alt)) return this.refIndex.get(alt);
    if (name === 'chrM' || name === 'M') { for (const n of ['MT', 'chrMT']) if (this.refIndex.has(n)) return this.refIndex.get(n); }
    return -1;
  }
  // Риды, перекрывающие [beg, end) на refId. onChunk — прогресс.
  async fetch(refId, beg, end) {
    const chunks = this.index.chunks(refId, beg, end);
    const results = await Promise.all(chunks.map(c => this._readChunk(c, refId, beg, end)));
    const reads = [];
    for (const r of results) for (const x of r) reads.push(x);
    return reads;
  }
  async _readChunk(c, refId, beg, end) {
    const cBeg = Math.floor(c[0] / 65536), uBeg = c[0] % 65536, cEnd = Math.floor(c[1] / 65536), uEnd = c[1] % 65536;
    const { data, offs } = await readBgzf(this.src, cBeg, cEnd);
    const dv = new DataView(data.buffer, data.byteOffset, data.byteLength);
    let p = (offs.get(cBeg) || 0) + uBeg;
    const stop = offs.has(cEnd) ? offs.get(cEnd) + uEnd : data.length;
    const out = [];
    while (p + 4 <= data.length && p < stop) {
      const bs = dv.getInt32(p, true);
      if (bs < 32 || p + 4 + bs > data.length) break;
      const refID = dv.getInt32(p + 4, true), pos = dv.getInt32(p + 8, true);
      if (refID !== refId) { if (refID > refId || refID < 0) break; p += 4 + bs; continue; }
      if (pos >= end) break;
      const rec = parseBamRecord(data, dv, p);
      p += 4 + bs;
      if (rec.flag & 4) continue;
      if (rec.end <= beg) continue;
      out.push(rec);
    }
    return out;
  }
}

function tryParseBamHeader(d) {
  if (d.length < 12) return null;
  if (d[0] !== 66 || d[1] !== 65 || d[2] !== 77 || d[3] !== 1) throw new CgvError('Это не BAM-файл (нет сигнатуры BAM)', 'format');
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

// Восстанавливает базы референса по тегу MD: заполняет refArr (индексы от winStart).
function fillRefFromMD(rec, refArr, winStart) {
  const md = recTags(rec).MD; if (typeof md !== 'string') return;
  // последовательность выровненных позиций: [refPos, queryPos | -1 для D]
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

// Признаки рида для отрисовки: несовпадения, вставки, делеции, пропуски, soft clip.
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

// База рида в позиции референса (код символа) или 0 — нет покрытия, 45 ('-') — делеция.
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
class IndexedFasta {
  constructor(faSource, faiSource) { this.src = faSource; this.faiSrc = faiSource; }
  async open() {
    const t = new TextDecoder().decode(await this.faiSrc.readAll());
    this.idx = new Map();
    for (const line of t.split('\n')) {
      const f = line.split('\t'); if (f.length < 5) continue;
      this.idx.set(f[0], { len: +f[1], off: +f[2], lb: +f[3], lw: +f[4] });
    }
    if (!this.idx.size) throw new CgvError('Пустой или неверный индекс FASTA (.fai)', 'format');
    return this;
  }
  entry(name) { return this.idx.get(name) || this.idx.get(name.startsWith('chr') ? name.slice(3) : 'chr' + name); }
  // Возвращает Uint8Array кодов символов (верхний регистр) для [beg, end).
  async fetch(name, beg, end) {
    const e = this.entry(name); if (!e) return null;
    beg = Math.max(0, beg); end = Math.min(e.len, end); if (end <= beg) return new Uint8Array(0);
    const o0 = e.off + Math.floor(beg / e.lb) * e.lw + (beg % e.lb);
    const o1 = e.off + Math.floor(end / e.lb) * e.lw + (end % e.lb);
    const raw = await this.src.read(o0, o1 - o0);
    const out = new Uint8Array(end - beg); let j = 0;
    for (let i = 0; i < raw.length && j < out.length; i++) { const c = raw[i]; if (c === 10 || c === 13) continue; out[j++] = c & ~32; }
    return out;
  }
}

// ======================= Модель =======================
const F_PAIRED = 1, F_PROPER = 2, F_REV = 16, F_MREV = 32, F_R1 = 64, F_R2 = 128, F_SECONDARY = 256, F_QCFAIL = 512, F_DUP = 1024, F_SUPP = 2048;

function passFilters(rec, flt) {
  if (rec.mapq < flt.minMapq) return false;
  if (flt.hideDuplicates && rec.flag & F_DUP) return false;
  if (flt.hideSecondary && rec.flag & F_SECONDARY) return false;
  if (flt.hideSupplementary && rec.flag & F_SUPP) return false;
  if (flt.hideQcFail && rec.flag & F_QCFAIL) return false;
  return true;
}

// Покрытие: 6 счётчиков на позицию (A C G T N del).
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

// Ориентация пары для библиотеки FR: 'FR' норма, 'RF', 'FF', 'RR', 'other' (мейт на другой хромосоме), 'unpaired', 'mateUnmapped'
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

// Порог «необычного» размера вставки по загруженным данным (как в IGV: перцентили).
function insertSizeBounds(reads) {
  const v = [];
  for (const r of reads) if ((r.flag & F_PROPER) && r.tlen > 0) v.push(r.tlen);
  if (v.length < 20) return { lo: 50, hi: 1000 };
  v.sort((a, b) => a - b);
  return { lo: v[Math.floor(v.length * 0.005)], hi: v[Math.floor(v.length * 0.995)] };
}

// Плотная раскладка по рядам: куча (min по концу ряда).
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

if (typeof module !== 'undefined') module.exports = { FileSource, UrlSource, CgvError, BamFile, BaiIndex, IndexedFasta, readBgzf, recTags, fillRefFromMD, computeFeatures, computeCoverage, baseAtRefPos, cigarString, packRows, passFilters, pairOrientation, insertSizeBounds, reg2bins };
