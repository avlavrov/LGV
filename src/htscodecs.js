// CRAM block codecs, adapted from htscodecs (JavaScript reference decoders by James Bonfield,
// Genome Research Ltd), https://github.com/samtools/htscodecs, commit 0a487b5, javascript/.
// Only decoding is taken. Changes from the original are listed in SECURITY.md and marked
// "CGV:" below. The original license follows.
//
// ---------------------------------------------------------------------------------------------
// Copyright (c) 2019-2020 Genome Research Ltd.
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

if (typeof module !== 'undefined') module.exports = { ransDecode4x8, HtsStream };
