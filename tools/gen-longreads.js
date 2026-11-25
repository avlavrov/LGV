// Writes a long-read test set for speed checks: one short circular-like contig (a plasmid) covered end
// to end by thousands of nanopore-like reads (substitutions, insertions, deletions every few dozen
// bases, base qualities). Not committed output; for profiling the drawing of long reads.
// Run: node tools/gen-longreads.js OUTDIR [reads] [contigLength]
'use strict';
const fs = require('fs'), path = require('path'), zlib = require('zlib');
const OUT = process.argv[2] || '.', N = +process.argv[3] || 5000, L = +process.argv[4] || 6891;
let seed = 12345;
const rnd = () => { seed = (seed + 0x6D2B79F5) | 0; let t = Math.imul(seed ^ (seed >>> 15), 1 | seed); t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296; };
const B = 'ACGT', ref = Array.from({ length: L }, () => B[Math.floor(rnd() * 4)]).join('');
const CODE = { A: 1, C: 2, G: 4, T: 8, N: 15 }, OPC = { M: 0, I: 1, D: 2, S: 4 };

function makeRead(i) {
  const start = Math.floor(rnd() * 200), end = L - Math.floor(rnd() * 200);
  const cig = [], seq = [], qual = []; let md = '', mdRun = 0, nm = 0, p = start;
  const push = (op, n) => { if (cig.length && cig[cig.length - 1][0] === op) cig[cig.length - 1][1] += n; else cig.push([op, n]); };
  if (rnd() < 0.5) { const n = 20 + Math.floor(rnd() * 80); for (let k = 0; k < n; k++) { seq.push(B[Math.floor(rnd() * 4)]); qual.push(10); } push('S', n); }
  while (p < end) {
    const r = rnd();
    if (r < 0.02) { // substitution
      let b; do { b = B[Math.floor(rnd() * 4)]; } while (b === ref[p]);
      seq.push(b); qual.push(8 + Math.floor(rnd() * 8)); push('M', 1); md += mdRun + ref[p]; mdRun = 0; nm++; p++;
    } else if (r < 0.035) { // insertion
      const n = 1 + Math.floor(rnd() * 3); for (let k = 0; k < n; k++) { seq.push(B[Math.floor(rnd() * 4)]); qual.push(9); } push('I', n); nm += n;
    } else if (r < 0.05) { // deletion
      const n = Math.min(end - p, 1 + Math.floor(rnd() * 3)); push('D', n); md += mdRun + '^' + ref.slice(p, p + n); mdRun = 0; nm += n; p += n;
    } else { seq.push(ref[p]); qual.push(15 + Math.floor(rnd() * 25)); push('M', 1); mdRun++; p++; }
  }
  md += mdRun;
  return { name: 'read' + i, pos: start, end: p, flag: rnd() < 0.5 ? 16 : 0, cig, seq: seq.join(''), qual, md, nm };
}
function reg2bin(beg, end) { --end; if (beg >> 14 === end >> 14) return ((1 << 15) - 1) / 7 + (beg >> 14); if (beg >> 17 === end >> 17) return ((1 << 12) - 1) / 7 + (beg >> 17); if (beg >> 20 === end >> 20) return ((1 << 9) - 1) / 7 + (beg >> 20); if (beg >> 23 === end >> 23) return ((1 << 6) - 1) / 7 + (beg >> 23); if (beg >> 26 === end >> 26) return ((1 << 3) - 1) / 7 + (beg >> 26); return 0; }
function encode(r) {
  const name = Buffer.from(r.name + '\0', 'latin1'), cig = Buffer.alloc(4 * r.cig.length);
  r.cig.forEach(([op, n], i) => cig.writeUInt32LE((n << 4) | OPC[op], 4 * i));
  const ls = r.seq.length, sq = Buffer.alloc((ls + 1) >> 1);
  for (let i = 0; i < ls; i++) sq[i >> 1] |= CODE[r.seq[i]] << (i & 1 ? 0 : 4);
  const tags = Buffer.concat([Buffer.from('NMi', 'latin1'), (() => { const b = Buffer.alloc(4); b.writeInt32LE(r.nm); return b; })(), Buffer.from('MDZ' + r.md + '\0', 'latin1')]);
  const core = Buffer.alloc(32);
  core.writeInt32LE(0, 0); core.writeInt32LE(r.pos, 4); core.writeUInt8(name.length, 8); core.writeUInt8(60, 9); core.writeUInt16LE(reg2bin(r.pos, r.end), 10);
  core.writeUInt16LE(r.cig.length, 12); core.writeUInt16LE(r.flag, 14); core.writeInt32LE(ls, 16); core.writeInt32LE(-1, 20); core.writeInt32LE(-1, 24); core.writeInt32LE(0, 28);
  const body = Buffer.concat([core, name, cig, sq, Buffer.from(r.qual), tags]), len = Buffer.alloc(4); len.writeInt32LE(body.length);
  return Buffer.concat([len, body]);
}
// BGZF writer that reports virtual offsets
const blocks = []; let filePos = 0, cur = [], curLen = 0;
function flushBlock() {
  if (!curLen) return;
  const data = Buffer.concat(cur), c = zlib.deflateRawSync(data, { level: 6 }), h = Buffer.from([31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, 66, 67, 2, 0, 0, 0]);
  h.writeUInt16LE(c.length + 25, 16);
  const t = Buffer.alloc(8); t.writeUInt32LE(zlib.crc32 ? zlib.crc32(data) : 0, 0); t.writeUInt32LE(data.length, 4);
  const blk = Buffer.concat([h, c, t]); blocks.push(blk); filePos += blk.length; cur = []; curLen = 0;
}
const voff = () => filePos * 65536 + curLen;
function write(b) { if (curLen + b.length > 60000) flushBlock(); cur.push(b); curLen += b.length; }
const text = `@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:plasmid\tLN:${L}\n@RG\tID:lr\tSM:LONGREADS\tPL:ONT\n`;
const hb = [Buffer.from('BAM\x01', 'latin1')], i4 = v => { const b = Buffer.alloc(4); b.writeInt32LE(v); return b; };
hb.push(i4(text.length), Buffer.from(text, 'latin1'), i4(1), i4('plasmid'.length + 1), Buffer.from('plasmid\0', 'latin1'), i4(L));
write(Buffer.concat(hb)); flushBlock();
const reads = Array.from({ length: N }, (_, i) => makeRead(i)).sort((a, b) => a.pos - b.pos);
const v0 = voff();
for (const r of reads) write(encode(r));
flushBlock(); const v1 = voff();
blocks.push(Buffer.from([31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, 66, 67, 2, 0, 27, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0]));
const u64 = v => { const b = Buffer.alloc(8); b.writeUInt32LE(v % 4294967296, 0); b.writeUInt32LE(Math.floor(v / 4294967296), 4); return b; };
const bins = new Map(); for (const r of reads) bins.set(reg2bin(r.pos, r.end), true);
// one chunk per bin covering all records (records of several bins interleave; fine for a test index)
const bai = [Buffer.from('BAI\x01', 'latin1'), i4(1), i4(bins.size)];
for (const b of bins.keys()) bai.push(i4(b), i4(1), u64(v0), u64(v1));
const nIntv = Math.ceil(L / 16384); bai.push(i4(nIntv)); for (let k = 0; k < nIntv; k++) bai.push(u64(v0));
fs.mkdirSync(OUT, { recursive: true });
fs.writeFileSync(path.join(OUT, 'longreads.bam'), Buffer.concat(blocks));
fs.writeFileSync(path.join(OUT, 'longreads.bam.bai'), Buffer.concat(bai));
let fa = '>plasmid\n'; for (let i = 0; i < L; i += 60) fa += ref.slice(i, i + 60) + '\n';
fs.writeFileSync(path.join(OUT, 'longreads.fa'), fa);
fs.writeFileSync(path.join(OUT, 'longreads.fa.fai'), `plasmid\t${L}\t9\t60\t61\n`);
console.log(`${N} reads, ${(Buffer.concat(blocks).length / 1048576).toFixed(1)} MB`);
