// CGV: bzip2 and xz/LZMA2 decompression for CRAM blocks (methods 2 and 3). Written for CGV from the
// format descriptions (bzip2: Seward's format as documented by the bzip2 sources and Joe Tsai's
// "bzip2 format specification"; xz: "The .xz File Format" 1.1.0; LZMA: Igor Pavlov's LzmaSpec).
// No third-party code.
//
// Safety rules, as in the rest of CGV: every read is bounds-checked; output never grows past the size
// the CRAM block states; every length, table and distance taken from the data is checked before
// use; checksums (bzip2 block and stream CRCs, xz CRC32/CRC64) are verified; any damage ends in a
// CgvError, never a hang, a crash or wrong output passed on.
'use strict';
if (typeof CgvError === 'undefined' && typeof require !== 'undefined') Object.assign(globalThis, require('./core.js'));

const ZBad = what => new CgvError('cramBad', { what });

// ======================= bzip2 =======================
// CRC32 as bzip2 uses it: polynomial 0x04C11DB7, most significant bit first
const BZ_CRC = (() => { const t = new Uint32Array(256); for (let i = 0; i < 256; i++) { let c = i << 24; for (let k = 0; k < 8; k++) c = c & 0x80000000 ? (c << 1) ^ 0x04C11DB7 : c << 1; t[i] = c >>> 0; } return t; })();
class BzBits {
  constructor(u8) { this.u8 = u8; this.pos = 0; this.buf = 0; this.n = 0; }
  bits(k) { // k ≤ 24
    while (this.n < k) { if (this.pos >= this.u8.length) throw ZBad('a bzip2 block is cut short'); this.buf = ((this.buf << 8) | this.u8[this.pos++]) >>> 0; this.n += 8; }
    this.n -= k; return (this.buf >>> this.n) & ((1 << k) - 1);
  }
  bit() { return this.bits(1); }
}
function bzip2Decode(src, expectSize) {
  const out = new Uint8Array(expectSize);
  let o = 0;
  const r = new BzBits(src);
  while (r.pos < src.length) {
    // one stream: "BZh" and the block size level
    if (r.bits(8) !== 0x42 || r.bits(8) !== 0x5A || r.bits(8) !== 0x68) throw ZBad('not a bzip2 stream');
    const level = r.bits(8) - 0x30;
    if (level < 1 || level > 9) throw ZBad('an impossible bzip2 block size');
    const maxBlock = level * 100000;
    let combined = 0;
    for (;;) {
      const m1 = r.bits(24), m2 = r.bits(24);
      if (m1 === 0x177245 && m2 === 0x385090) { // end of stream
        const crc = (r.bits(16) * 65536 + r.bits(16)) >>> 0;
        if (crc !== combined) throw ZBad('bzip2 stream checksum mismatch');
        r.n -= r.n % 8; // streams start on a byte boundary
        break;
      }
      if (m1 !== 0x314159 || m2 !== 0x265359) throw ZBad('a damaged bzip2 block header');
      const blockCrc = (r.bits(16) * 65536 + r.bits(16)) >>> 0;
      if (r.bit()) throw ZBad('randomised bzip2 blocks (bzip2 before 0.9.5) are not read');
      const origPtr = r.bits(24);
      // symbols in use
      const used = r.bits(16), seqToUnseq = [];
      for (let i = 0; i < 16; i++) if (used & (0x8000 >> i)) { const g = r.bits(16); for (let j = 0; j < 16; j++) if (g & (0x8000 >> j)) seqToUnseq.push(i * 16 + j); }
      if (!seqToUnseq.length) throw ZBad('a bzip2 block uses no symbols');
      const alpha = seqToUnseq.length + 2;
      const nGroups = r.bits(3);
      if (nGroups < 2 || nGroups > 6) throw ZBad('an impossible number of bzip2 tables');
      const nSel = r.bits(15);
      if (nSel < 1) throw ZBad('a bzip2 block without selectors');
      const mtfG = [0, 1, 2, 3, 4, 5].slice(0, nGroups), selectors = new Uint8Array(Math.min(nSel, 18002));
      for (let i = 0; i < nSel; i++) {
        let j = 0; while (r.bit()) { if (++j >= nGroups) throw ZBad('a bzip2 selector out of range'); }
        const v = mtfG[j]; mtfG.splice(j, 1); mtfG.unshift(v);
        if (i < selectors.length) selectors[i] = v; // more than 18002 selectors: the rest are unused (as bzip2 1.0.8 does)
      }
      // Huffman tables, as delta-coded lengths
      const tables = [];
      for (let t = 0; t < nGroups; t++) {
        const lens = new Uint8Array(alpha);
        let len = r.bits(5);
        for (let i = 0; i < alpha; i++) {
          for (;;) { if (len < 1 || len > 20) throw ZBad('an impossible bzip2 code length'); if (!r.bit()) break; len += r.bit() ? -1 : 1; }
          lens[i] = len;
        }
        tables.push(bzTable(lens));
      }
      // symbols → MTF / run-length → bytes of the BWT block
      const tt = new Uint32Array(maxBlock);
      const mtf = new Uint8Array(256); for (let i = 0; i < 256; i++) mtf[i] = i;
      const cnt = new Uint32Array(256);
      const EOB = alpha - 1;
      let n = 0, groupPos = 0, sel = 0, tab = null, run = 0, runBit = 1;
      for (;;) {
        if (groupPos === 0) { if (sel >= selectors.length) throw ZBad('bzip2 data runs past its selectors'); tab = tables[selectors[sel++]]; groupPos = 50; }
        groupPos--;
        const sym = bzSym(r, tab);
        if (sym <= 1) { // RUNA, RUNB: a run of the front symbol
          run += (sym + 1) * runBit; runBit *= 2;
          if (run > maxBlock) throw ZBad('a bzip2 run past the block size');
          continue;
        }
        if (run) {
          if (n + run > maxBlock) throw ZBad('a bzip2 block larger than its size level');
          const b = seqToUnseq[mtf[0]]; cnt[b] += run; while (run--) tt[n++] = b; run = 0; runBit = 1;
        }
        if (sym === EOB) break;
        const k = sym - 1;
        if (k >= seqToUnseq.length) throw ZBad('a bzip2 symbol out of range');
        const v = mtf[k]; mtf.copyWithin(1, 0, k); mtf[0] = v;
        if (n >= maxBlock) throw ZBad('a bzip2 block larger than its size level');
        const b = seqToUnseq[v]; cnt[b]++; tt[n++] = b;
      }
      if (origPtr >= n) throw ZBad('a bzip2 block start out of range');
      // inverse BWT: tt[i] keeps the byte in its low 8 bits and the next index above them
      const cf = new Uint32Array(256); for (let i = 0, s = 0; i < 256; i++) { cf[i] = s; s += cnt[i]; }
      for (let i = 0; i < n; i++) { const b = tt[i] & 255; tt[cf[b]++] |= i << 8; }
      // undo the initial run-length coding (4 equal bytes and a count) while checking the CRC
      let crc = 0xFFFFFFFF, p = tt[origPtr] >>> 8, last = -1, same = 0;
      for (let i = 0; i < n; i++) {
        const e = tt[p]; const b = e & 255; p = e >>> 8;
        if (same === 4) {
          if (o + b > out.length) throw ZBad('bzip2 data larger than its block');
          for (let k = 0; k < b; k++) { out[o++] = last; crc = ((crc << 8) ^ BZ_CRC[((crc >>> 24) ^ last) & 255]) >>> 0; }
          same = 0; last = -1; continue;
        }
        if (b === last) same++; else { same = 1; last = b; }
        if (o >= out.length) throw ZBad('bzip2 data larger than its block');
        out[o++] = b; crc = ((crc << 8) ^ BZ_CRC[((crc >>> 24) ^ b) & 255]) >>> 0;
      }
      crc = (~crc) >>> 0;
      if (crc !== blockCrc) throw ZBad('bzip2 block checksum mismatch');
      combined = (((combined << 1) | (combined >>> 31)) ^ crc) >>> 0;
    }
    if (r.pos >= src.length) break;
  }
  if (o !== expectSize) throw ZBad('bzip2 data shorter than its block');
  return out;
}
// Canonical codes of one bzip2 table: limit / base / permutation per length.
function bzTable(lens) {
  let min = 32, max = 0; for (const l of lens) { if (l < min) min = l; if (l > max) max = l; }
  const perm = [], limit = new Int32Array(22).fill(-1), base = new Int32Array(22);
  for (let l = min; l <= max; l++) for (let i = 0; i < lens.length; i++) if (lens[i] === l) perm.push(i);
  let code = 0, idx = 0;
  for (let l = min; l <= max; l++) {
    let c = 0; for (const x of lens) if (x === l) c++;
    base[l] = idx - code; code += c; idx += c; limit[l] = code - 1; code <<= 1;
  }
  return { min, max, limit, base, perm };
}
function bzSym(r, t) {
  let l = t.min, v = r.bits(l);
  for (;;) {
    if (v <= t.limit[l]) { const i = t.base[l] + v; if (i < 0 || i >= t.perm.length) throw ZBad('a bzip2 code not in its table'); return t.perm[i]; }
    if (++l > t.max) throw ZBad('a bzip2 code not in its table');
    v = (v << 1) | r.bit();
  }
}

// ======================= xz / LZMA2 =======================
const XZ_CRC64 = (() => { // CRC-64/XZ (ECMA-182, reflected), as two 32-bit halves
  const lo = new Uint32Array(256), hi = new Uint32Array(256);
  const PL = 0xD7870F42, PH = 0xC96C5795;
  for (let i = 0; i < 256; i++) {
    let l = i, h = 0;
    for (let k = 0; k < 8; k++) { const b = l & 1; l = ((l >>> 1) | ((h & 1) << 31)) >>> 0; h >>>= 1; if (b) { l = (l ^ PL) >>> 0; h = (h ^ PH) >>> 0; } }
    lo[i] = l; hi[i] = h;
  }
  return { lo, hi };
})();
function xzCrc64(u8, a, b) {
  let l = 0xFFFFFFFF, h = 0xFFFFFFFF;
  for (let i = a; i < b; i++) { const k = (l ^ u8[i]) & 255; l = (((l >>> 8) | (h << 24)) ^ XZ_CRC64.lo[k]) >>> 0; h = ((h >>> 8) ^ XZ_CRC64.hi[k]) >>> 0; }
  return [(~l) >>> 0, (~h) >>> 0];
}
const XZ_CRC32 = (() => { const t = new Uint32Array(256); for (let n = 0; n < 256; n++) { let c = n; for (let k = 0; k < 8; k++) c = c & 1 ? 0xEDB88320 ^ (c >>> 1) : c >>> 1; t[n] = c >>> 0; } return t; })();
function xzCrc32(u8, a, b) { let c = 0xFFFFFFFF; for (let i = a; i < b; i++) c = XZ_CRC32[(c ^ u8[i]) & 255] ^ (c >>> 8); return (c ^ 0xFFFFFFFF) >>> 0; }
const le32 = (u, p) => (u[p] | (u[p + 1] << 8) | (u[p + 2] << 16) | (u[p + 3] << 24)) >>> 0;

function xzDecode(src, expectSize) {
  const out = new Uint8Array(expectSize);
  let p = 0, o = 0;
  const need = n => { if (p + n > src.length) throw ZBad('an xz stream is cut short'); };
  const varint = () => { let v = 0, m = 1; for (let i = 0; i < 9; i++) { need(1); const b = src[p++]; v += (b & 0x7F) * m; m *= 128; if (!(b & 0x80)) return v; } throw ZBad('an over-long number in an xz stream'); };
  while (p < src.length) {
    if (src[p] === 0 && o === expectSize) { p++; continue; } // stream padding
    // stream header
    need(12);
    const magic = [0xFD, 0x37, 0x7A, 0x58, 0x5A, 0x00];
    for (let i = 0; i < 6; i++) if (src[p + i] !== magic[i]) throw ZBad('not an xz stream');
    if (src[p + 6] !== 0 || (src[p + 7] & 0xF0)) throw ZBad('unknown xz stream flags');
    const check = src[p + 7] & 0x0F;
    if (xzCrc32(src, p + 6, p + 8) !== le32(src, p + 8)) throw ZBad('xz header checksum mismatch');
    const checkSize = check === 0 ? 0 : check === 1 ? 4 : check === 4 ? 8 : check === 10 ? 32 : -1;
    if (checkSize < 0) throw ZBad('an unknown xz check type ' + check);
    p += 12;
    // blocks, up to the index
    for (;;) {
      need(1);
      if (src[p] === 0) break; // the index follows
      const hs = (src[p] + 1) * 4, hStart = p;
      need(hs);
      if (xzCrc32(src, hStart, hStart + hs - 4) !== le32(src, hStart + hs - 4)) throw ZBad('xz block header checksum mismatch');
      p++;
      const flags = src[p++];
      if (flags & 0x3C) throw ZBad('unknown xz block flags');
      const nFilters = (flags & 3) + 1;
      let compSize = -1, uncompSize = -1;
      if (flags & 0x40) compSize = varint();
      if (flags & 0x80) uncompSize = varint();
      let dictSize = -1;
      for (let f = 0; f < nFilters; f++) {
        const id = varint(), psize = varint();
        need(psize);
        if (id !== 0x21 || f !== nFilters - 1) throw ZBad('xz filter ' + id + ' is not supported (LZMA2 only)');
        if (psize !== 1) throw ZBad('damaged LZMA2 filter properties');
        const d = src[p];
        if (d > 40) throw ZBad('an impossible LZMA2 dictionary size');
        dictSize = d === 40 ? 0xFFFFFFFF : (2 | (d & 1)) * 2 ** (11 + (d >> 1));
        p += psize;
      }
      p = hStart + hs; // header padding and its CRC
      const dataStart = p, oStart = o;
      const r = lzma2Decode(src, p, out, o);
      p = r.p; o = r.o;
      if (compSize >= 0 && p - dataStart !== compSize) throw ZBad('xz block size mismatch');
      if (uncompSize >= 0 && o - oStart !== uncompSize) throw ZBad('xz block output size mismatch');
      while ((p - dataStart) & 3) { need(1); if (src[p++] !== 0) throw ZBad('xz block padding is not zero'); }
      need(checkSize);
      if (check === 1 && xzCrc32(out, oStart, o) !== le32(src, p)) throw ZBad('xz block checksum mismatch');
      if (check === 4) { const [l, h] = xzCrc64(out, oStart, o); if (l !== le32(src, p) || h !== le32(src, p + 4)) throw ZBad('xz block checksum mismatch'); }
      // SHA-256 (check 10) is not verified; CRAM's own block CRC32 still covers the compressed data
      p += checkSize;
      void dictSize;
    }
    // index: record count, then (unpadded, uncompressed) sizes; padding and CRC32
    const iStart = p; p++;
    const nRec = varint();
    if (nRec > 1e6) throw ZBad('an impossible xz index');
    for (let i = 0; i < nRec; i++) { varint(); varint(); }
    while ((p - iStart) & 3) { need(1); if (src[p++] !== 0) throw ZBad('xz index padding is not zero'); }
    need(4); if (xzCrc32(src, iStart, p) !== le32(src, p)) throw ZBad('xz index checksum mismatch'); p += 4;
    // stream footer
    need(12);
    if (src[p + 10] !== 0x59 || src[p + 11] !== 0x5A) throw ZBad('a damaged xz stream footer');
    if (xzCrc32(src, p + 4, p + 10) !== le32(src, p)) throw ZBad('xz footer checksum mismatch');
    p += 12;
  }
  if (o !== expectSize) throw ZBad('xz data shorter than its block');
  return out;
}

// LZMA2 chunks of one xz block, decoded into `out` from position o. Returns the new positions.
function lzma2Decode(src, p, out, o) {
  const dec = new LzmaDecoder(out);
  let needProps = true, needDictReset = true, first = true;
  dec.dictStart = o;
  for (;;) {
    if (p >= src.length) throw ZBad('an LZMA2 stream is cut short');
    const c = src[p++];
    if (c === 0) return { p, o };
    if (c === 1 || c === 2) { // uncompressed chunk (1: with a dictionary reset)
      if (p + 2 > src.length) throw ZBad('an LZMA2 chunk is cut short');
      const n = ((src[p] << 8) | src[p + 1]) + 1; p += 2;
      if (c === 1) { dec.dictStart = o; needDictReset = false; } else if (needDictReset) throw ZBad('LZMA2 data before its dictionary reset');
      if (p + n > src.length || o + n > out.length) throw ZBad('an LZMA2 chunk past the end of its block');
      out.set(src.subarray(p, p + n), o); p += n; o += n;
      first = false;
      continue;
    }
    if (c < 0x80) throw ZBad('an unknown LZMA2 chunk type');
    if (p + 4 > src.length) throw ZBad('an LZMA2 chunk is cut short');
    const unpacked = ((c & 0x1F) << 16) + ((src[p] << 8) | src[p + 1]) + 1, packed = ((src[p + 2] << 8) | src[p + 3]) + 1;
    p += 4;
    const reset = (c >> 5) & 3;
    if (reset === 3) { dec.dictStart = o; needDictReset = false; } else if (needDictReset) throw ZBad('LZMA2 data before its dictionary reset');
    if (reset >= 2) {
      if (p >= src.length) throw ZBad('an LZMA2 chunk is cut short');
      let d = src[p++];
      if (d >= 9 * 5 * 5) throw ZBad('impossible LZMA properties');
      const lc = d % 9; d = (d / 9) | 0; const lp = d % 5, pb = (d / 5) | 0;
      if (lc + lp > 4) throw ZBad('impossible LZMA properties');
      dec.setProps(lc, lp, pb); needProps = false;
    } else if (needProps) throw ZBad('LZMA2 data before its properties');
    if (reset >= 1) dec.resetState();
    else if (first) throw ZBad('LZMA2 data before its state reset');
    if (p + packed > src.length || o + unpacked > out.length) throw ZBad('an LZMA2 chunk past the end of its block');
    o = dec.decodeChunk(src, p, p + packed, o, o + unpacked);
    p += packed;
    first = false;
  }
}

// The LZMA decoder (after LzmaSpec): range coder, adaptive bit models, literal/match/rep decoding.
// The dictionary is the output array itself: matches may only reach back to dictStart.
const LZ_BITS = 11, LZ_ONE = 1 << LZ_BITS, LZ_MOVE = 5, LZ_TOP = 1 << 24;
class LzmaDecoder {
  constructor(out) { this.out = out; this.dictStart = 0; this.lc = 3; this.lp = 0; this.pb = 2; this.lit = null; }
  setProps(lc, lp, pb) { this.lc = lc; this.lp = lp; this.pb = pb; this.lit = new Uint16Array(0x300 << (lc + lp)); }
  resetState() {
    this.lit.fill(LZ_ONE >> 1);
    const P = n => new Uint16Array(n).fill(LZ_ONE >> 1);
    this.isMatch = P(12 << 4); this.isRep = P(12); this.isRepG0 = P(12); this.isRepG1 = P(12); this.isRepG2 = P(12); this.isRep0Long = P(12 << 4);
    this.posSlot = P(4 << 6); this.posDec = P(115); this.align = P(16);
    this.lenDec = this.newLen(); this.repLenDec = this.newLen();
    this.state = 0; this.reps = [0, 0, 0, 0];
  }
  newLen() { const P = n => new Uint16Array(n).fill(LZ_ONE >> 1); return { choice: P(2), low: P(16 << 3), mid: P(16 << 3), high: P(256) }; }
  // range coder over src[p .. end)
  rcInit(src, p, end) {
    if (end - p < 5 || src[p] !== 0) throw ZBad('a damaged LZMA range coder start');
    this.src = src; this.p = p + 1; this.end = end; this.range = 0xFFFFFFFF; this.code = 0;
    for (let i = 0; i < 4; i++) this.code = ((this.code << 8) | src[this.p++]) >>> 0;
    if (this.code === this.range) throw ZBad('a damaged LZMA range coder start');
  }
  norm() { if (this.range < LZ_TOP) { if (this.p >= this.end) throw ZBad('an LZMA chunk is cut short'); this.range = (this.range << 8) >>> 0; this.code = ((this.code << 8) | this.src[this.p++]) >>> 0; } }
  bit(probs, i) {
    const v = probs[i], bound = (this.range >>> LZ_BITS) * v;
    let b;
    if (this.code < bound) { this.range = bound; probs[i] = v + ((LZ_ONE - v) >>> LZ_MOVE); b = 0; }
    else { this.range -= bound; this.code -= bound; probs[i] = v - (v >>> LZ_MOVE); b = 1; }
    this.norm();
    return b;
  }
  direct(n) { let r = 0; for (let i = 0; i < n; i++) { this.range >>>= 1; let t = 0; if (this.code >= this.range) { this.code -= this.range; t = 1; } r = (r << 1) | t; this.norm(); } return r >>> 0; }
  tree(probs, off, bits) { let m = 1; for (let i = 0; i < bits; i++) m = (m << 1) | this.bit(probs, off + m); return m - (1 << bits); }
  revTree(probs, off, bits) { let m = 1, s = 0; for (let i = 0; i < bits; i++) { const b = this.bit(probs, off + m); m = (m << 1) | b; s |= b << i; } return s; }
  len(L, posState) {
    if (!this.bit(L.choice, 0)) return this.tree(L.low, posState << 3, 3);
    if (!this.bit(L.choice, 1)) return 8 + this.tree(L.mid, posState << 3, 3);
    return 16 + this.tree(L.high, 0, 8);
  }
  decodeChunk(src, p, end, o, oEnd) {
    this.rcInit(src, p, end);
    const out = this.out, pbMask = (1 << this.pb) - 1, lpMask = (1 << this.lp) - 1, lc = this.lc;
    let state = this.state; const reps = this.reps;
    while (o < oEnd) {
      const posState = o & pbMask;
      if (!this.bit(this.isMatch, (state << 4) + posState)) { // literal
        const prev = o > this.dictStart ? out[o - 1] : 0;
        const base = 0x300 * (((o & lpMask) << lc) + (prev >> (8 - lc)));
        let sym = 1;
        if (state >= 7) { // after a match: decode against the byte at rep0
          const mi = o - reps[0] - 1; if (mi < this.dictStart) throw ZBad('an LZMA match before the start of its data');
          let match = out[mi];
          do { const mb = (match >> 7) & 1; match <<= 1; const b = this.bit(this.lit, base + ((1 + mb) << 8) + sym); sym = (sym << 1) | b; if (mb !== b) break; } while (sym < 0x100);
        }
        while (sym < 0x100) sym = (sym << 1) | this.bit(this.lit, base + sym);
        out[o++] = sym & 255;
        state = state < 4 ? 0 : state < 10 ? state - 3 : state - 6;
        continue;
      }
      let len;
      if (this.bit(this.isRep, state)) {
        if (o === this.dictStart) throw ZBad('an LZMA repeat at the start of its data');
        if (!this.bit(this.isRepG0, state)) {
          if (!this.bit(this.isRep0Long, (state << 4) + posState)) { // one byte at rep0
            state = state < 7 ? 9 : 11;
            const mi = o - reps[0] - 1; if (mi < this.dictStart) throw ZBad('an LZMA match before the start of its data');
            out[o] = out[mi]; o++; continue;
          }
        } else {
          let d;
          if (!this.bit(this.isRepG1, state)) d = reps[1];
          else { if (!this.bit(this.isRepG2, state)) d = reps[2]; else { d = reps[3]; reps[3] = reps[2]; } reps[2] = reps[1]; }
          reps[1] = reps[0]; reps[0] = d;
        }
        len = this.len(this.repLenDec, posState);
        state = state < 7 ? 8 : 11;
      } else {
        reps[3] = reps[2]; reps[2] = reps[1]; reps[1] = reps[0];
        len = this.len(this.lenDec, posState);
        state = state < 7 ? 7 : 10;
        // distance
        const ls = len < 4 ? len : 3;
        const slot = this.tree(this.posSlot, ls << 6, 6);
        let dist;
        if (slot < 4) dist = slot;
        else {
          const nb = (slot >>> 1) - 1;
          dist = ((2 | (slot & 1)) << nb) >>> 0;
          if (slot < 14) dist += this.revTree(this.posDec, dist - slot, nb);
          else { dist += this.direct(nb - 4) * 16; dist += this.revTree(this.align, 0, 4); }
          dist >>>= 0;
        }
        if (dist === 0xFFFFFFFF) throw ZBad('an LZMA end marker inside an LZMA2 chunk');
        reps[0] = dist;
      }
      len += 2;
      const from = o - reps[0] - 1;
      if (from < this.dictStart) throw ZBad('an LZMA match before the start of its data');
      if (o + len > oEnd) throw ZBad('an LZMA match past the end of its chunk');
      for (let k = 0; k < len; k++) out[o + k] = out[from + k]; // may overlap: byte by byte on purpose
      o += len;
    }
    this.state = state;
    return o;
  }
}

if (typeof module !== 'undefined') module.exports = { bzip2Decode, xzDecode };
