// CRAM block codecs, adapted from htscodecs (JavaScript reference decoders by James Bonfield,
// Genome Research Ltd), https://github.com/samtools/htscodecs, commit 0a487b5, javascript/
// (rans.js, rans4x16.js, arith_gen.js, arith_sh.js, byte_model.js, fqzcomp.js, tok3.js, iostream.js).
// Only decoding is taken. Changes from the original are listed in SECURITY.md and marked
// "CGV:" below. The original license follows.
//
// ---------------------------------------------------------------------------------------------
// Copyright (c) 2019, 2020 Genome Research Ltd.
// Author(s): James Bonfield
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    1. Redistributions of source code must retain the above copyright notice,
//       this list of conditions and the following disclaimer.
//
//    2. Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//    3. Neither the names Genome Research Ltd and Wellcome Trust Sanger
//       Institute nor the names of its contributors may be used to endorse
//       or promote products derived from this software without specific
//       prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH
// LTD OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// ---------------------------------------------------------------------------------------------
'use strict';
if (typeof CgvError === 'undefined' && typeof require !== 'undefined') Object.assign(globalThis, require('./core.js'));

// CGV: errors from the codecs are CgvErrors, so a damaged block gives a clear message.
const HtsError = what => new CgvError('cramBad', { what });

// CGV: IOStream (iostream.js) without bounds checks returned `undefined` past the end of the data,
// and frequency tables are read in loops that end only on a 0 byte, so a short block could loop
// forever. Every read here is checked and throws instead. Writing methods are not needed.
class HtsStream {
  constructor(buf, pos) { this.buf = buf; this.pos = pos || 0; this.length = buf.length; }
  EOF() { return this.pos >= this.length; }
  need(n) { if (this.pos + n > this.length) throw HtsError('a compressed block is cut short'); }
  ReadByte() { this.need(1); return this.buf[this.pos++]; }
  ReadData(len) { this.need(len); const a = this.buf.subarray(this.pos, this.pos + len); this.pos += len; return a; }
  ReadUint16() { this.need(2); const b = this.buf, p = this.pos; this.pos += 2; return b[p] | (b[p + 1] << 8); }
  // CGV: unsigned (the original used Buffer.readInt32LE)
  ReadUint32() { this.need(4); const b = this.buf, p = this.pos; this.pos += 4; return (b[p] | (b[p + 1] << 8) | (b[p + 2] << 16) | (b[p + 3] << 24)) >>> 0; }
  ReadUint7() {
    let i = 0, c, n = 0;
    do { c = this.ReadByte(); i = i * 128 + (c & 0x7f); if (++n > 5) throw HtsError('an over-long number in a compressed block'); } while (c & 0x80);
    return i;
  }
  ReadITF8() {
    let i = this.ReadByte();
    if (i >= 0xf0) { this.need(4); const b = this.buf, p = this.pos; i = ((i & 0x0f) << 28) + (b[p] << 20) + (b[p + 1] << 12) + (b[p + 2] << 4) + (b[p + 3] >> 4); this.pos += 4; }
    else if (i >= 0xe0) { this.need(3); const b = this.buf, p = this.pos; i = ((i & 0x0f) << 24) + (b[p] << 16) + (b[p + 1] << 8) + b[p + 2]; this.pos += 3; }
    else if (i >= 0xc0) { this.need(2); const b = this.buf, p = this.pos; i = ((i & 0x1f) << 16) + (b[p] << 8) + b[p + 1]; this.pos += 2; }
    else if (i >= 0x80) { this.need(1); i = ((i & 0x3f) << 8) + this.buf[this.pos++]; }
    return i;
  }
}

// ======================= rANS 4x8 (rans.js) =======================
// CGV: a symbol lookup that refuses frequencies summing to more than 4096 and marks the codes
// above the total as invalid (the original stepped past the table and wrote symbol 256).
function RansBuildC2S(C) {
  if (C[256] > 0x1000) throw HtsError('rANS frequencies sum to more than 4096');
  const C2S = new Int16Array(0x1000).fill(-1);
  let s = 0;
  for (let f = 0; f < C[256]; f++) { while (f >= C[s + 1]) s++; C2S[f] = s; }
  return C2S;
}
function RansAdvanceStep(R, c, f) { return f * (R >>> 12) + (R & 0xfff) - c; }
function RansRenorm(src, R) { while (R < (1 << 23)) R = R * 256 + src.ReadByte(); return R; }

function ReadFrequencies0(src, F, C) {
  for (let i = 0; i < 256; i++) F[i] = 0;
  let sym = src.ReadByte(), last_sym = sym, rle = 0;
  do {
    const f = src.ReadITF8();
    if (f < 0 || f > 0x1000) throw HtsError('an impossible rANS frequency'); // CGV
    F[sym] = f;
    if (rle > 0) { rle--; sym++; if (sym > 255) throw HtsError('a rANS frequency table runs past 255'); } // CGV
    else { sym = src.ReadByte(); if (sym === last_sym + 1) rle = src.ReadByte(); }
    last_sym = sym;
  } while (sym !== 0);
  C[0] = 0;
  for (let i = 0; i <= 255; i++) C[i + 1] = C[i] + F[i];
}
function RansDecode0(src, nbytes) {
  const F = new Array(256), C = new Array(257);
  ReadFrequencies0(src, F, C);
  const C2S = RansBuildC2S(C);
  const R = new Array(4);
  for (let i = 0; i < 4; i++) R[i] = src.ReadUint32();
  const output = new Uint8Array(nbytes); // CGV: Uint8Array instead of Buffer
  for (let i = 0; i < nbytes; i++) {
    const i4 = i % 4, f = R[i4] & 0xfff, s = C2S[f];
    if (s < 0) throw HtsError('a rANS code outside its frequency table'); // CGV
    output[i] = s;
    R[i4] = RansRenorm(src, RansAdvanceStep(R[i4], C[s], F[s]));
  }
  return output;
}
function ReadFrequencies1(src, F, C) {
  for (let i = 0; i < 256; i++) { F[i] = new Array(256).fill(0); C[i] = new Array(257).fill(0); }
  let sym = src.ReadByte(), last_sym = sym, rle = 0;
  do {
    ReadFrequencies0(src, F[sym], C[sym]);
    if (rle > 0) { rle--; sym++; if (sym > 255) throw HtsError('a rANS context table runs past 255'); } // CGV
    else { sym = src.ReadByte(); if (sym === last_sym + 1) rle = src.ReadByte(); }
    last_sym = sym;
  } while (sym !== 0);
}
function RansDecode1(src, nbytes) {
  const F = new Array(256), C = new Array(256);
  ReadFrequencies1(src, F, C);
  // CGV: tables only for contexts that occur; an unused context with no frequencies is an error if hit
  const C2S = C.map(c => (c[256] ? RansBuildC2S(c) : null));
  const R = new Array(4), L = [0, 0, 0, 0];
  for (let j = 0; j < 4; j++) R[j] = src.ReadUint32();
  const output = new Uint8Array(nbytes);
  const nbytes4 = Math.floor(nbytes / 4);
  const step = (j, at) => {
    const t = C2S[L[j]]; if (!t) throw HtsError('a rANS context without frequencies');
    const s = t[R[j] & 0xfff]; if (s < 0) throw HtsError('a rANS code outside its frequency table');
    output[at] = s;
    R[j] = RansRenorm(src, RansAdvanceStep(R[j], C[L[j]][s], F[L[j]][s]));
    L[j] = s;
  };
  let i = 0;
  for (; i < nbytes4; i++) for (let j = 0; j < 4; j++) step(j, i + j * nbytes4);
  // the remainder, when the size is not a multiple of 4, uses state 3 alone
  for (i = 4 * i; i < nbytes; i++) step(3, i);
  return output;
}
// Block decode entry point (rans.js decode). CGV: the stated sizes are checked against the block's.
function ransDecode4x8(src, expectSize) {
  const stream = new HtsStream(src);
  const order = stream.ReadByte(), n_in = stream.ReadUint32(), n_out = stream.ReadUint32();
  if (order > 1) throw HtsError('an unknown rANS order ' + order);
  if (n_in !== src.length - 9 || n_out !== expectSize) throw HtsError('rANS sizes do not match the block');
  return order === 0 ? RansDecode0(stream, n_out) : RansDecode1(stream, n_out);
}


// CGV: a cap on sizes taken from the data, checked before allocating
const HTS_MAX_OUT = 256 << 20;
const htsSize = (n, what) => { if (!(n >= 0) || n > HTS_MAX_OUT) throw HtsError('an impossible ' + what + ' size in a compressed block'); return n; };

// ======================= rANS 4x16 (rans4x16.js) =======================
function Rans16BuildC2S(C, bits) {
  const max = 1 << bits;
  if (C[256] > max) throw HtsError('rANS frequencies sum to more than ' + max); // CGV
  const C2S = new Int16Array(max).fill(-1);
  let s = 0;
  for (let f = 0; f < C[256]; f++) { while (f >= C[s + 1]) s++; C2S[f] = s; }
  return C2S;
}
function Rans16Advance(R, c, f, bits) { return f * (R >>> bits) + (R & ((1 << bits) - 1)) - c; }
function Rans16Renorm(src, R) { if (R < (1 << 15)) R = R * 65536 + src.ReadUint16(); return R; }

function Rans16DecodeRLEMeta(src, N, depth) {
  const u_meta_len = htsSize(src.ReadUint7(), 'RLE'), rle_len = htsSize(src.ReadUint7(), 'RLE');
  let rle_meta;
  if (u_meta_len & 1) rle_meta = src.ReadData((u_meta_len - 1) / 2);
  else { const comp_meta_len = src.ReadUint7(); rle_meta = Rans16Decode0(new HtsStream(src.ReadData(comp_meta_len)), u_meta_len / 2, N); }
  const meta = new HtsStream(rle_meta), L = new Uint8Array(256);
  let n = meta.ReadByte(); if (n === 0) n = 256;
  for (let i = 0; i < n; i++) L[meta.ReadByte()] = 1;
  return [L, meta, rle_len];
}
function Rans16DecodeRLE(buf, L, rle_meta, len) {
  const out = new Uint8Array(len);
  let j = 0;
  for (let i = 0; j < len; i++) {
    if (i >= buf.length) throw HtsError('RLE data run short'); // CGV
    const sym = buf[i];
    if (L[sym]) {
      const run = rle_meta.ReadUint7();
      if (run >= len - j) throw HtsError('an RLE run past the end of its block'); // CGV: the original could spin for 2^35 steps
      for (let r = 0; r <= run; r++) out[j++] = sym;
    } else out[j++] = sym;
  }
  return out;
}
function Rans16DecodePackMeta(src) {
  const nsym = src.ReadByte(), P = new Uint8Array(Math.max(nsym, 1));
  for (let i = 0; i < nsym; i++) P[i] = src.ReadByte();
  return [P, nsym, htsSize(src.ReadUint7(), 'packed')];
}
function htsDecodePack(data, P, nsym, len) {
  const out = new Uint8Array(len);
  let bits, mask;
  if (nsym <= 1) { out.fill(P[0]); return out; }
  else if (nsym <= 2) { bits = 1; mask = 1; } else if (nsym <= 4) { bits = 2; mask = 3; } else if (nsym <= 16) { bits = 4; mask = 15; }
  else throw HtsError('a packed block with ' + nsym + ' symbols'); // CGV: the original returned uninitialised memory
  const per = 8 / bits;
  if (data.length < Math.ceil(len / per)) throw HtsError('packed data run short');
  for (let i = 0, j = 0, v = 0; i < len; i++) { if (i % per === 0) v = data[j++]; out[i] = P[v & mask]; v >>= bits; }
  return out;
}
function Rans16DecodeStripe(src, len, depth) {
  const N = src.ReadByte();
  if (N < 1) throw HtsError('a striped block with no streams');
  for (let j = 0; j < N; j++) src.ReadUint7(); // compressed lengths, not needed to decode
  const T = [], ulen = [];
  for (let j = 0; j < N; j++) {
    ulen[j] = Math.floor(len / N) + ((len % N) > j ? 1 : 0);
    T[j] = Rans16DecodeStream(src, ulen[j], depth + 1);
    if (T[j].length < ulen[j]) throw HtsError('a stripe decodes short');
  }
  const out = new Uint8Array(len);
  for (let j = 0; j < N; j++) for (let i = 0; i < ulen[j]; i++) out[i * N + j] = T[j][i];
  return out;
}
function Rans16DecodeStream(stream, n_out, depth) {
  if (depth > 4) throw HtsError('compressed blocks nested too deep'); // CGV
  const format = stream.ReadByte();
  const order = format & 1, x32 = format & 4, stripe = format & 8, nosz = format & 16, cat = format & 32, rle = format & 64, pack = format & 128;
  const Nway = x32 ? 32 : 4;
  if (!nosz) n_out = htsSize(stream.ReadUint7(), 'output');
  if (stripe) return Rans16DecodeStripe(stream, n_out, depth);
  let P, nsym, pack_len, L, rle_meta, rle_len;
  if (pack) { pack_len = n_out; [P, nsym, n_out] = Rans16DecodePackMeta(stream); }
  if (rle) { rle_len = n_out; [L, rle_meta, n_out] = Rans16DecodeRLEMeta(stream, Nway, depth); }
  let buf;
  if (cat) buf = stream.ReadData(n_out);
  else if (order === 0) buf = Rans16Decode0(stream, n_out, Nway);
  else buf = Rans16Decode1(stream, n_out, Nway);
  if (rle) buf = Rans16DecodeRLE(buf, L, rle_meta, rle_len);
  if (pack) buf = htsDecodePack(buf, P, nsym, pack_len);
  return buf;
}
function Rans16ReadAlphabet(src) {
  const A = new Uint8Array(256);
  let rle = 0, sym = src.ReadByte(), last_sym = sym;
  do {
    A[sym] = 1;
    if (rle > 0) { rle--; sym++; if (sym > 255) throw HtsError('a rANS alphabet runs past 255'); } // CGV
    else { sym = src.ReadByte(); if (sym === last_sym + 1) rle = src.ReadByte(); }
    last_sym = sym;
  } while (sym !== 0);
  return A;
}
function Rans16Normalise(F, bits) {
  let tot = 0; for (let i = 0; i < 256; i++) tot += F[i];
  if (tot === 0 || tot === (1 << bits)) return;
  if (tot > (1 << bits)) throw HtsError('rANS frequencies sum to more than ' + (1 << bits)); // CGV
  let shift = 0; while (tot < (1 << bits)) { tot *= 2; shift++; }
  for (let i = 0; i < 256; i++) F[i] <<= shift;
}
function Rans16ReadFrequencies0(src, F, C) {
  F.fill(0);
  const A = Rans16ReadAlphabet(src);
  for (let i = 0; i < 256; i++) if (A[i]) F[i] = src.ReadUint7();
  Rans16Normalise(F, 12);
  C[0] = 0; for (let i = 0; i <= 255; i++) C[i + 1] = C[i] + F[i];
}
function Rans16Decode0(src, nbytes, N) {
  htsSize(nbytes, 'output');
  const F = new Array(256).fill(0), C = new Array(257);
  Rans16ReadFrequencies0(src, F, C);
  const C2S = Rans16BuildC2S(C, 12), R = new Array(N);
  for (let i = 0; i < N; i++) R[i] = src.ReadUint32();
  const output = new Uint8Array(nbytes);
  for (let i = 0; i < nbytes; i++) {
    const ix = i & (N - 1), s = C2S[R[ix] & 0xfff];
    if (s < 0) throw HtsError('a rANS code outside its frequency table'); // CGV
    output[i] = s;
    R[ix] = Rans16Renorm(src, Rans16Advance(R[ix], C[s], F[s], 12));
  }
  return output;
}
function Rans16ReadFrequencies1(src, F, C, shift) {
  for (let i = 0; i < 256; i++) { F[i] = new Array(256).fill(0); C[i] = new Array(257).fill(0); }
  const A = Rans16ReadAlphabet(src);
  for (let i = 0; i < 256; i++) {
    if (!A[i]) continue;
    let run = 0;
    for (let j = 0; j < 256; j++) {
      if (!A[j]) continue;
      if (run > 0) run--;
      else { F[i][j] = src.ReadUint7(); if (F[i][j] === 0) run = src.ReadByte(); }
    }
    Rans16Normalise(F[i], shift);
    C[i][0] = 0; for (let j = 0; j < 256; j++) C[i][j + 1] = C[i][j] + F[i][j];
  }
}
function Rans16Decode1(src, nbytes, N) {
  htsSize(nbytes, 'output');
  const comp = src.ReadByte(), shift = comp >> 4;
  if (shift < 1 || shift > 12) throw HtsError('an impossible rANS frequency precision'); // CGV
  let freq_src = src;
  if (comp & 1) {
    const ulen = htsSize(src.ReadUint7(), 'table'), clen = src.ReadUint7();
    freq_src = new HtsStream(Rans16Decode0(new HtsStream(src.ReadData(clen)), ulen, 4));
  }
  const F = new Array(256), C = new Array(256);
  Rans16ReadFrequencies1(freq_src, F, C, shift);
  const C2S = C.map(c => (c[256] ? Rans16BuildC2S(c, shift) : null)); // CGV: only contexts in use
  const R = new Array(N), L = new Array(N).fill(0);
  for (let j = 0; j < N; j++) R[j] = src.ReadUint32();
  const output = new Uint8Array(nbytes), mask = (1 << shift) - 1, nbytesx = Math.floor(nbytes / N);
  const step = (j, at) => {
    const t = C2S[L[j]]; if (!t) throw HtsError('a rANS context without frequencies');
    const s = t[R[j] & mask]; if (s < 0) throw HtsError('a rANS code outside its frequency table');
    output[at] = s;
    R[j] = Rans16Renorm(src, Rans16Advance(R[j], C[L[j]][s], F[L[j]][s], shift));
    L[j] = s;
  };
  let i = 0;
  for (; i < nbytesx; i++) for (let j = 0; j < N; j++) step(j, i + j * nbytesx);
  for (i = N * i; i < nbytes; i++) step(N - 1, i);
  return output;
}
function ransDecode4x16(src, expectSize) {
  const out = Rans16DecodeStream(new HtsStream(src), 0, 0);
  if (expectSize !== undefined && out.length !== expectSize) throw HtsError('rANS 4x16 output size differs from the block');
  return out;
}

// ======================= Adaptive arithmetic coder (arith_sh.js, byte_model.js) =======================
class RangeCoder {
  constructor() { this.low = 0; this.range = 0xffffffff; this.code = 0; }
  RangeStartDecode(src) { for (let i = 0; i < 5; i++) this.code = (this.code * 256 + src.ReadByte()) % 4294967296; }
  // CGV: >>> 0 in place of Math.floor; both values stay below 2^32 in a valid stream
  RangeGetFrequency(tot_freq) {
    this.range = (this.range / tot_freq) >>> 0;
    if (this.range === 0) throw HtsError('arithmetic coder range ran out'); // CGV: else the renorm loop below never ends
    return (this.code / this.range) >>> 0;
  }
  RangeDecode(src, sym_low, sym_freq) {
    this.code -= sym_low * this.range;
    this.range *= sym_freq;
    if (!(this.range > 0)) throw HtsError('arithmetic coder range ran out'); // CGV
    if (this.range < (1 << 24)) {
      // CGV: the byte reads of HtsStream.ReadByte, inlined; still bounds-checked
      const buf = src.buf; let pos = src.pos, range = this.range, code = this.code;
      while (range < (1 << 24)) { if (pos >= src.length) throw HtsError('a compressed block is cut short'); range *= 256; code = code * 256 + buf[pos++]; }
      src.pos = pos; this.range = range; this.code = code;
    }
  }
}
const BM_MAX_FREQ = (1 << 16) - 17, BM_STEP = 16;
class ByteModel {
  constructor(max_sym = 256) {
    this.total_freq = max_sym; this.max_sym = max_sym - 1;
    this.S = new Uint16Array(max_sym); this.F = new Uint32Array(max_sym);
    for (let i = 0; i <= this.max_sym; i++) { this.S[i] = i; this.F[i] = 1; }
  }
  ModelDecode(src, rc) {
    const freq = rc.RangeGetFrequency(this.total_freq), F = this.F, S = this.S, max = this.max_sym;
    let acc = 0, x = 0;
    while (x <= max && acc + F[x] <= freq) acc += F[x++];
    if (x > max) throw HtsError('an arithmetic code outside its model'); // CGV
    rc.RangeDecode(src, acc, F[x]);
    F[x] += BM_STEP; this.total_freq += BM_STEP;
    if (this.total_freq > BM_MAX_FREQ) this.ModelRenormalise();
    const sym = S[x];
    if (x > 0 && F[x] > F[x - 1]) {
      let t = F[x]; F[x] = F[x - 1]; F[x - 1] = t;
      t = S[x]; S[x] = S[x - 1]; S[x - 1] = t;
    }
    return sym;
  }
  ModelRenormalise() {
    this.total_freq = 0;
    for (let i = 0; i <= this.max_sym; i++) { this.F[i] -= Math.floor(this.F[i] / 2); this.total_freq += this.F[i]; }
  }
}

// ======================= Arithmetic codec (arith_gen.js) =======================
const ARITH_ORDER = 1, ARITH_EXT = 4, ARITH_STRIPE = 8, ARITH_NOSIZE = 16, ARITH_CAT = 32, ARITH_RLE = 64, ARITH_PACK = 128;
function arithDecodeStream(stream, n_out, depth) {
  if (depth > 4) throw HtsError('compressed blocks nested too deep'); // CGV
  const flags = stream.ReadByte();
  if (!(flags & ARITH_NOSIZE)) n_out = htsSize(stream.ReadUint7(), 'output');
  let e_len = n_out;
  const order = flags & ARITH_ORDER;
  if (flags & ARITH_STRIPE) {
    const N = stream.ReadByte(); if (N < 1) throw HtsError('a striped block with no streams');
    for (let j = 0; j < N; j++) stream.ReadUint7();
    const T = [], ulen = [];
    for (let j = 0; j < N; j++) { ulen[j] = Math.floor(n_out / N) + ((n_out % N) > j ? 1 : 0); T[j] = arithDecodeStream(stream, ulen[j], depth + 1); if (T[j].length < ulen[j]) throw HtsError('a stripe decodes short'); }
    const out = new Uint8Array(n_out);
    for (let j = 0; j < N; j++) for (let i = 0; i < ulen[j]; i++) out[i * N + j] = T[j][i];
    return out;
  }
  let P, nsym;
  if (flags & ARITH_PACK) { nsym = stream.ReadByte(); P = new Uint8Array(Math.max(nsym, 1)); for (let i = 0; i < nsym; i++) P[i] = stream.ReadByte(); e_len = htsSize(stream.ReadUint7(), 'packed'); }
  let data;
  if (flags & ARITH_CAT) data = stream.ReadData(e_len);
  else if (flags & ARITH_EXT) throw new CgvError('cramCodec', { codec: 'bzip2 inside arithmetic coding' }); // CGV: bzip2 comes in a later stage
  else if (flags & ARITH_RLE) data = arithDecodeRLE(stream, e_len, order);
  else data = arithDecodeOrder(stream, e_len, order);
  if (flags & ARITH_PACK) data = nsym > 16 ? data : htsDecodePack(data, P, nsym, n_out);
  return data;
}
function arithDecodeOrder(stream, n_out, order) {
  const output = new Uint8Array(n_out);
  let max_sym = stream.ReadByte(); if (max_sym === 0) max_sym = 256;
  const models = order ? new Array(max_sym) : [new ByteModel(max_sym)]; // CGV: order-1 models made on first use
  const rc = new RangeCoder(); rc.RangeStartDecode(stream);
  let last = 0;
  for (let i = 0; i < n_out; i++) {
    const m = order ? (models[last] || (models[last] = new ByteModel(max_sym))) : models[0];
    output[i] = m.ModelDecode(stream, rc);
    last = output[i];
    if (order && last >= max_sym) throw HtsError('an arithmetic symbol above its maximum');
  }
  return output;
}
function arithDecodeRLE(stream, n_out, order) {
  const output = new Uint8Array(n_out);
  let max_sym = stream.ReadByte(); if (max_sym === 0) max_sym = 256;
  const model_lit = order ? new Array(max_sym) : [new ByteModel(max_sym)];
  const model_run = new Array(258);
  const run = k => model_run[k] || (model_run[k] = new ByteModel(4));
  const rc = new RangeCoder(); rc.RangeStartDecode(stream);
  let last = 0, i = 0;
  while (i < n_out) {
    const ml = order ? (model_lit[last] || (model_lit[last] = new ByteModel(max_sym))) : model_lit[0];
    const sym = output[i] = ml.ModelDecode(stream, rc);
    if (order) { last = sym; if (last >= max_sym) throw HtsError('an arithmetic symbol above its maximum'); }
    let part = run(sym).ModelDecode(stream, rc), len = part, rctx = 256;
    while (part === 3) {
      part = run(rctx).ModelDecode(stream, rc); rctx = 257; len += part;
      if (len > n_out) throw HtsError('an RLE run past the end of its block'); // CGV: a run that costs no input could loop forever
    }
    if (i + len >= n_out + 1) throw HtsError('an RLE run past the end of its block');
    for (let j = 1; j <= len; j++) output[i + j] = sym;
    i += len + 1;
  }
  return output;
}
function arithDecode(src, expectSize) {
  const out = arithDecodeStream(new HtsStream(src), 0, 0);
  if (expectSize !== undefined && out.length !== expectSize) throw HtsError('arithmetic output size differs from the block');
  return out;
}

// ======================= FQZcomp qualities (fqzcomp.js) =======================
const FQZ_FLAG_DEDUP = 2, FQZ_FLAG_FLEN = 4, FQZ_FLAG_SEL = 8, FQZ_FLAG_QMAP = 16, FQZ_FLAG_PTAB = 32, FQZ_FLAG_DTAB = 64, FQZ_FLAG_QTAB = 128;
const FQZ_GFLAG_MULTI_PARAM = 1, FQZ_GFLAG_HAVE_STAB = 2, FQZ_GFLAG_DO_REV = 4;
function fqzReadArray(src, tab, size) {
  let j = 0, z = 0, last = -1;
  const R = [];
  while (z < size) {
    const run = src.ReadByte(); R[j++] = run; z += run;
    if (run === last) { let copy = src.ReadByte(); z += run * copy; while (copy--) R[j++] = run; }
    last = run;
    if (j > 4 * size + 16) throw HtsError('an FQZ table does not end'); // CGV
  }
  let i = 0; j = 0; z = 0;
  while (z < size) {
    let run_len = 0, part;
    do { if (j >= R.length) throw HtsError('an FQZ table runs short'); part = R[j++]; run_len += part; } while (part === 255); // CGV
    while (run_len-- && z < size) tab[z++] = i;
    i++;
  }
}
function fqzUpdateCtx(params, state, q) {
  let last = params.context;
  state.qctx = (state.qctx << params.qshift) + params.qtab[q];
  last += (state.qctx & ((1 << params.qbits) - 1)) << params.qloc;
  if (params.do_pos) last += params.ptab[Math.min(state.p, 1023)] << params.ploc;
  if (params.do_delta) { last += params.dtab[Math.min(state.delta, 255)] << params.dloc; state.delta += (state.prevq !== q) ? 1 : 0; state.prevq = q; }
  if (params.do_sel) last += state.s << params.sloc;
  state.p--;
  return last & 0xffff;
}
function fqzSingleParam(src) {
  const p = {};
  p.context = src.ReadUint16(); p.pflags = src.ReadByte();
  p.do_dedup = p.pflags & FQZ_FLAG_DEDUP; p.fixed_len = p.pflags & FQZ_FLAG_FLEN; p.do_sel = p.pflags & FQZ_FLAG_SEL;
  p.do_qmap = p.pflags & FQZ_FLAG_QMAP; p.do_pos = p.pflags & FQZ_FLAG_PTAB; p.do_delta = p.pflags & FQZ_FLAG_DTAB; p.do_qtab = p.pflags & FQZ_FLAG_QTAB;
  p.max_sym = src.ReadByte();
  let x = src.ReadByte(); p.qbits = x >> 4; p.qshift = x & 15;
  x = src.ReadByte(); p.qloc = x >> 4; p.sloc = x & 15;
  x = src.ReadByte(); p.ploc = x >> 4; p.dloc = x & 15;
  p.qmap = new Uint8Array(256);
  if (p.pflags & FQZ_FLAG_QMAP) for (let i = 0; i < p.max_sym; i++) p.qmap[i] = src.ReadByte();
  else for (let i = 0; i < 256; i++) p.qmap[i] = i;
  p.qtab = new Uint16Array(1024);
  if (p.qbits > 0 && (p.pflags & FQZ_FLAG_QTAB)) fqzReadArray(src, p.qtab, 256);
  else for (let i = 0; i < 256; i++) p.qtab[i] = i;
  p.ptab = new Uint16Array(1024);
  if (p.pflags & FQZ_FLAG_PTAB) fqzReadArray(src, p.ptab, 1024);
  p.dtab = new Uint16Array(256);
  if (p.pflags & FQZ_FLAG_DTAB) fqzReadArray(src, p.dtab, 256);
  return p;
}
function fqzParams(src) {
  const g = { max_sym: 0 };
  const vers = src.ReadByte();
  if (vers !== 5) throw HtsError('an unknown FQZcomp version ' + vers); // CGV: the original returned nothing
  const gflags = src.ReadByte();
  const nparam = (gflags & FQZ_GFLAG_MULTI_PARAM) ? src.ReadByte() : 1;
  if (nparam < 1) throw HtsError('FQZcomp without parameters');
  let max_sel = nparam > 1 ? nparam - 1 : 0; // CGV: the original read gflags.nparam, always undefined
  const stab = new Uint16Array(256);
  if (gflags & FQZ_GFLAG_HAVE_STAB) { max_sel = src.ReadByte(); fqzReadArray(src, stab, 256); }
  else { let i = 0; for (; i < nparam; i++) stab[i] = i; for (; i < 256; i++) stab[i] = nparam - 1; }
  g.do_rev = gflags & FQZ_GFLAG_DO_REV; g.stab = stab; g.max_sel = max_sel; g.params = [];
  for (let p = 0; p < nparam; p++) { g.params[p] = fqzSingleParam(src); if (g.max_sym < g.params[p].max_sym) g.max_sym = g.params[p].max_sym; }
  for (let s = 0; s <= max_sel; s++) if (stab[s] >= nparam) throw HtsError('an FQZcomp selector without parameters');
  return g;
}
function fqzDecode(src, expectSize) {
  const stream = new HtsStream(src);
  const n_out = htsSize(stream.ReadUint7(), 'output');
  const g = fqzParams(stream);
  const rev = [], q_lens = [];
  // CGV: the 65,536 quality models are made on first use, not all up front
  const qual = new Array(1 << 16).fill(null), qsyms = g.max_sym + 1;
  const mlen = [0, 1, 2, 3].map(() => new ByteModel(256)), mrev = new ByteModel(2), mdup = new ByteModel(2);
  const msel = g.max_sel > 0 ? new ByteModel(g.max_sel + 1) : null;
  const rc = new RangeCoder(); rc.RangeStartDecode(stream);
  const output = new Uint8Array(n_out);
  const st = { qctx: 0, prevq: 0, delta: 0, p: 0, s: 0, x: 0, len: 0, rec: 0 };
  let i = 0, params = null, last = 0;
  while (i < n_out) {
    if (st.p === 0) {
      // a new record (decode_fqz_new_record)
      st.s = msel ? msel.ModelDecode(stream, rc) : 0;
      st.x = g.stab[st.s];
      params = g.params[st.x];
      let len;
      if (params.fixed_len >= 0) {
        len = mlen[0].ModelDecode(stream, rc) | (mlen[1].ModelDecode(stream, rc) << 8) | (mlen[2].ModelDecode(stream, rc) << 16) | (mlen[3].ModelDecode(stream, rc) << 24);
        len >>>= 0;
        if (params.fixed_len > 0) params.fixed_len = -len;
      } else len = -params.fixed_len;
      if (len === 0 || len > n_out - i) throw HtsError('an FQZcomp record length outside its block'); // CGV
      st.len = len;
      if (g.do_rev) rev[st.rec] = mrev.ModelDecode(stream, rc);
      st.rec++;
      // CGV: a duplicate record copies the previous one. The original decoded the duplicate flag twice
      // and copied with an undefined length; this follows htscodecs' C decoder (fqzcomp_qual.c).
      if ((params.pflags & FQZ_FLAG_DEDUP) && mdup.ModelDecode(stream, rc)) {
        if (i < len) throw HtsError('an FQZcomp duplicate without a record before it');
        for (let x = 0; x < len; x++) output[i + x] = output[i + x - len];
        q_lens.push(len); i += len; continue;
      }
      q_lens.push(len);
      st.p = len; st.delta = 0; st.qctx = 0; st.prevq = 0;
      last = params.context;
    }
    let qmod = qual[last]; if (qmod === null) qmod = qual[last] = new ByteModel(qsyms);
    const Q = qmod.ModelDecode(stream, rc);
    output[i++] = params.qmap[Q];
    last = fqzUpdateCtx(params, st, Q);
  }
  if (g.do_rev) {
    for (let rec = 0, at = 0; at < n_out && rec < q_lens.length; at += q_lens[rec++]) {
      if (!rev[rec]) continue;
      for (let j = 0, k = q_lens[rec] - 1; j < k; j++, k--) { const t = output[at + j]; output[at + j] = output[at + k]; output[at + k] = t; }
    }
  }
  if (expectSize !== undefined && n_out !== expectSize) throw HtsError('FQZcomp output size differs from the block');
  return output;
}

// ======================= Read-name tokeniser (tok3.js) =======================
const TOK_TYPE = 0, TOK_STRING = 1, TOK_CHAR = 2, TOK_DIGITS0 = 3, TOK_DZLEN = 4, TOK_DUP = 5, TOK_DIFF = 6, TOK_DIGITS = 7, TOK_DELTA = 8, TOK_DELTA0 = 9, TOK_MATCH = 10, TOK_NOP = 11, TOK_END = 12;
function tokStreams(src, use_arith, nnames) {
  let t = -1;
  const B = [];
  while (!src.EOF()) {
    const ttype = src.ReadByte(), tok_new = ttype & 128, tok_dup = ttype & 64, type = ttype & 63;
    if (type > TOK_END) throw HtsError('an unknown name token type ' + type); // CGV
    if (tok_new) { t++; if (t > 255) throw HtsError('too many name tokens'); B[t] = new Array(13); }
    if (t < 0) throw HtsError('name token data before its first token'); // CGV
    if (type !== TOK_TYPE && tok_new) { const M = new Uint8Array(nnames).fill(TOK_MATCH); M[0] = type; B[t][TOK_TYPE] = new HtsStream(M); }
    if (tok_dup) {
      const dup_pos = src.ReadByte(), dup_type = src.ReadByte();
      if (!B[dup_pos] || dup_type > TOK_END || !B[dup_pos][dup_type]) throw HtsError('a name token copies a missing stream'); // CGV
      B[t][type] = new HtsStream(B[dup_pos][dup_type].buf);
    } else {
      const clen = src.ReadUint7(), data = src.ReadData(clen);
      B[t][type] = new HtsStream(use_arith ? arithDecode(data) : ransDecode4x16(data));
    }
  }
  return B;
}
function tokPad(val, len) { let s = String(val); while (s.length < len) s = '0' + s; return s; }
function tokName(B, N, T, n) {
  const tok = (t, k) => { const s = B[t] && B[t][k]; if (!s) throw HtsError('a name refers to a missing token stream'); return s; }; // CGV
  let type = tok(0, TOK_TYPE).ReadByte();
  const dist = tok(0, type).ReadUint32(), m = n - dist;
  if (m < 0 || m > n || (dist === 0 && type === TOK_DUP)) throw HtsError('a name refers to one that is not there'); // CGV
  if (type === TOK_DUP) { N[n] = N[m]; T[n] = T[m]; return N[n]; }
  let t = 1;
  N[n] = ''; T[n] = [];
  const prev = T[m] || [];
  do {
    type = tok(t, TOK_TYPE).ReadByte();
    let v;
    switch (type) {
      case TOK_CHAR: v = String.fromCharCode(tok(t, TOK_CHAR).ReadByte()); break;
      case TOK_STRING: { const s = tok(t, TOK_STRING); v = ''; for (let b = s.ReadByte(); b; b = s.ReadByte()) v += String.fromCharCode(b); break; }
      case TOK_DIGITS: v = tok(t, TOK_DIGITS).ReadUint32(); break;
      case TOK_DIGITS0: { const d = tok(t, TOK_DIGITS0).ReadUint32(), l = tok(t, TOK_DZLEN).ReadByte(); v = tokPad(d, l); break; }
      case TOK_DELTA: v = (Number(prev[t]) >>> 0) + tok(t, TOK_DELTA).ReadByte(); break;
      case TOK_DELTA0: { const d = (Number(prev[t]) >>> 0) + tok(t, TOK_DELTA0).ReadByte(); v = tokPad(d, String(prev[t] === undefined ? '' : prev[t]).length); break; }
      case TOK_MATCH: v = prev[t] === undefined ? '' : prev[t]; break;
      default: v = '';
    }
    T[n][t] = v;
    N[n] += v;
    if (++t > 256) throw HtsError('a name with too many tokens'); // CGV
  } while (type !== TOK_END);
  return N[n];
}
// Names separated by NUL bytes, as CRAM stores them in the RN block.
function tok3Decode(src, expectSize) {
  const s = new HtsStream(src);
  const ulen = htsSize(s.ReadUint32(), 'names'), nnames = s.ReadUint32(), use_arith = s.ReadByte();
  if (nnames > ulen || nnames > 1 << 24) throw HtsError('an impossible number of read names'); // CGV
  const B = tokStreams(s, use_arith, nnames), N = new Array(nnames), T = new Array(nnames);
  const out = new Uint8Array(ulen); let o = 0;
  for (let i = 0; i < nnames; i++) {
    const name = tokName(B, N, T, i);
    if (o + name.length + 1 > ulen) throw HtsError('read names longer than their block');
    for (let k = 0; k < name.length; k++) out[o++] = name.charCodeAt(k) & 255;
    out[o++] = 0;
  }
  if (o !== ulen || (expectSize !== undefined && ulen !== expectSize)) throw HtsError('read names do not fill their block');
  return out;
}

if (typeof module !== 'undefined') module.exports = { ransDecode4x8, ransDecode4x16, arithDecode, fqzDecode, tok3Decode, HtsStream };
