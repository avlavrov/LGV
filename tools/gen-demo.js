// Synthetic demo data for CGV: a made-up reference, read pairs with known variants,
// BAM (BGZF) + BAI + FASTA + FAI. Node built-in modules only.
// CGV_L=<length> sets the main contig length (default 40000).
'use strict';
const zlib = require('zlib');
const fs = require('fs');
const path = require('path');

const OUT = process.argv[2] || path.join(__dirname, '..', 'demo');
fs.mkdirSync(OUT, { recursive: true });

// ---------- deterministic random numbers ----------
let seed = 20261016;
function rnd() { seed |= 0; seed = (seed + 0x6D2B79F5) | 0; let t = Math.imul(seed ^ (seed >>> 15), 1 | seed); t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296; }
function rint(a, b) { return a + Math.floor(rnd() * (b - a + 1)); }
function gauss() { let u = 0, v = 0; while (u === 0) u = rnd(); while (v === 0) v = rnd(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v); }
const BASES = 'ACGT';
function randBase() { return BASES[Math.floor(rnd() * 4)]; }
function otherBase(b) { let x; do { x = randBase(); } while (x === b); return x; }
function revcomp(s) { const m = { A: 'T', C: 'G', G: 'C', T: 'A', N: 'N' }; let o = ''; for (let i = s.length - 1; i >= 0; i--) o += m[s[i]]; return o; }

// ---------- reference ----------
const L = +(process.env.CGV_L || 40000), L2 = 2000;
// CGV_OFFSET=<bp> places the data in the middle of a long contig padded with N (sparse, panel-like files)
const OFF = +(process.env.CGV_OFFSET || 0), LTOT = L + 2 * OFF;
const sh = (contig, p) => contig === 0 ? p + OFF : p;
let ref = ''; for (let i = 0; i < L; i++) ref += randBase();
// segmental duplication: 30000–30600 is copied to 33000–33600, so reads there get low MAPQ
ref = ref.slice(0, 33000) + ref.slice(30000, 30600) + ref.slice(33600);
let ref2 = ''; for (let i = 0; i < L2; i++) ref2 += randBase();

// ---------- variants (0-based coordinates) ----------
const V = {
  hetSnv: { pos: 12000 }, homSnv: { pos: 15000 },
  hetDel: { pos: 18000, len: 3 }, hetIns: { pos: 21000, seq: 'TGA' },
  bigDel: { pos: 25000, len: 3000 },
};
V.hetSnv.alt = otherBase(ref[V.hetSnv.pos]);
V.homSnv.alt = otherBase(ref[V.homSnv.pos]);

// haplotype = sequence + map from haplotype positions to reference (-1 = insertion)
function buildHap(h) {
  let seq = [], map = [];
  for (let i = 0; i < L; i++) {
    if (h === 0 && i >= V.bigDel.pos && i < V.bigDel.pos + V.bigDel.len) continue;
    if (h === 0 && i >= V.hetDel.pos && i < V.hetDel.pos + V.hetDel.len) continue;
    let b = ref[i];
    if (i === V.homSnv.pos) b = V.homSnv.alt;
    if (h === 1 && i === V.hetSnv.pos) b = V.hetSnv.alt;
    seq.push(b); map.push(i);
    if (h === 1 && i === V.hetIns.pos) for (const c of V.hetIns.seq) { seq.push(c); map.push(-1); }
  }
  return { seq: seq.join(''), map };
}
const haps = [buildHap(0), buildHap(1)];

// ---------- align a haplotype segment to the reference ----------
// returns {pos, cigar:[[op,len]], seq}; op: M I D S
function alignSegment(hap, a, len) {
  const ops = []; let pos = -1, prev = -1;
  const push = (op, n) => { if (ops.length && ops[ops.length - 1][0] === op) ops[ops.length - 1][1] += n; else ops.push([op, n]); };
  for (let k = a; k < a + len; k++) {
    const r = hap.map[k];
    if (r < 0) { push('I', 1); continue; }
    if (pos < 0) pos = r;
    else if (r > prev + 1) push('D', r - prev - 1);
    push('M', 1); prev = r;
  }
  // an insertion at either end becomes a soft clip
  if (ops[0][0] === 'I') ops[0][0] = 'S';
  if (ops[ops.length - 1][0] === 'I') ops[ops.length - 1][0] = 'S';
  // a large deletion: like a real aligner, soft-clip the shorter side
  for (let i = 0; i < ops.length; i++) {
    if (ops[i][0] === 'D' && ops[i][1] > 50) {
      const left = ops.slice(0, i).reduce((s, o) => s + (o[0] !== 'D' ? o[1] : 0), 0);
      const right = ops.slice(i + 1).reduce((s, o) => s + (o[0] !== 'D' ? o[1] : 0), 0);
      if (left >= right) { ops.splice(i, ops.length - i, ['S', right]); }
      else {
        // the read now starts on the right side
        let rp = pos; for (let j = 0; j < i; j++) if (ops[j][0] === 'M' || ops[j][0] === 'D') rp += ops[j][1];
        rp += ops[i][1]; pos = rp;
        ops.splice(0, i + 1, ['S', left]);
      }
      break;
    }
  }
  return { pos, cigar: ops, seq: hap.seq.slice(a, a + len) };
}

// ---------- read pairs ----------
const RL = 150;
const reads = [];
let pairId = 0;
// Binned qualities (2, 12, 23, 37), like current Illumina instruments: realistic and compress well.
function qualFor(isErr) { return isErr ? 12 : (rnd() < 0.08 ? 23 : 37); }
function makeRead(al, rev, contig, extra) {
  let s = al.seq.split(''); const q = [];
  for (let i = 0; i < s.length; i++) { const err = rnd() < 0.003; if (err) s[i] = otherBase(s[i]); q.push(qualFor(err)); }
  return Object.assign({ contig, pos: al.pos, cigar: al.cigar, seq: s.join(''), qual: q, rev }, extra);
}
function inRepeat(p, e) { return (p < 30600 && e > 30000) || (p < 33600 && e > 33000); }
const DEPTH = +(process.env.CGV_DEPTH || 30); // CGV_DEPTH=<x> for deep targeted-panel style data
const nPairs = Math.round(DEPTH * L / (2 * RL));
for (let n = 0; n < nPairs; n++) {
  const h = rnd() < 0.5 ? 0 : 1, hap = haps[h];
  let frag = Math.round(350 + 40 * gauss()); frag = Math.max(220, Math.min(600, frag));
  const a = rint(0, hap.seq.length - frag - 1);
  const r1 = alignSegment(hap, a, RL), r2 = alignSegment(hap, a + frag - RL, RL);
  const name = 'demo_' + (++pairId);
  const e1 = r1.pos + r1.cigar.reduce((s, o) => s + (o[0] === 'M' || o[0] === 'D' ? o[1] : 0), 0);
  const e2 = r2.pos + r2.cigar.reduce((s, o) => s + (o[0] === 'M' || o[0] === 'D' ? o[1] : 0), 0);
  let mq = 60;
  if (inRepeat(r1.pos, e1) || inRepeat(r2.pos, e2)) mq = [0, 0, 0, 1, 3][rint(0, 4)];
  else if (rnd() < 0.05) mq = rint(15, 45);
  const swap = rnd() < 0.5; // which read is first in pair
  const A = makeRead(r1, false, 0, { name, mapq: mq, hp: h + 1, first: !swap });
  const B = makeRead(r2, true, 0, { name, mapq: mq, hp: h + 1, first: swap });
  A.mate = B; B.mate = A; reads.push(A, B);
  if (rnd() < 0.03) { // PCR duplicate
    const nm = 'dup_' + pairId;
    const C = Object.assign({}, A, { name: nm, dup: true }), D = Object.assign({}, B, { name: nm, dup: true });
    C.mate = D; D.mate = C; reads.push(C, D);
  }
}
// chimeric pairs: reads near 22000 with mates on demo2
for (let n = 0; n < 24; n++) {
  const a = rint(21700, 22200);
  const r1 = alignSegment(haps[0], a, RL);
  const b = rint(900, 1100);
  const r2 = { pos: b, cigar: [['M', RL]], seq: revcomp(ref2.slice(b, b + RL)) };
  const name = 'chim_' + n;
  const A = makeRead(r1, false, 0, { name, mapq: 60, hp: 1, first: true });
  const B = makeRead(r2, true, 1, { name, mapq: 60, hp: 1, first: false });
  B.seq = ref2.slice(b, b + RL); // stored in reference orientation
  A.mate = B; B.mate = A; reads.push(A, B);
}

// ---------- BAM encoding ----------
const CODE = { '=': 0, A: 1, C: 2, M: 3, G: 4, R: 5, S: 6, V: 7, T: 8, W: 9, Y: 10, H: 11, K: 12, D: 13, B: 14, N: 15 };
const OPC = { M: 0, I: 1, D: 2, N: 3, S: 4, H: 5, P: 6, '=': 7, X: 8 };
const refs = [ref, ref2];
function endOf(r) { return r.pos + r.cigar.reduce((s, o) => s + (o[0] === 'M' || o[0] === 'D' || o[0] === 'N' ? o[1] : 0), 0); }
function mdTag(r) {
  const R = refs[r.contig]; let md = '', run = 0, rp = r.pos, qp = 0;
  for (const [op, n] of r.cigar) {
    if (op === 'M') { for (let i = 0; i < n; i++) { if (r.seq[qp + i] === R[rp + i]) run++; else { md += run + R[rp + i]; run = 0; } } rp += n; qp += n; }
    else if (op === 'I' || op === 'S') qp += n;
    else if (op === 'D') { md += run + '^' + R.slice(rp, rp + n); run = 0; rp += n; }
  }
  return md + run;
}
function reg2bin(beg, end) {
  --end;
  if (beg >> 14 === end >> 14) return ((1 << 15) - 1) / 7 + (beg >> 14);
  if (beg >> 17 === end >> 17) return ((1 << 12) - 1) / 7 + (beg >> 17);
  if (beg >> 20 === end >> 20) return ((1 << 9) - 1) / 7 + (beg >> 20);
  if (beg >> 23 === end >> 23) return ((1 << 6) - 1) / 7 + (beg >> 23);
  if (beg >> 26 === end >> 26) return ((1 << 3) - 1) / 7 + (beg >> 26);
  return 0;
}
function encode(r) {
  const end = endOf(r), m = r.mate;
  let flag = 1 | (r.rev ? 16 : 0) | (m.rev ? 32 : 0) | (r.first ? 64 : 128);
  if (r.dup) flag |= 1024;
  let tlen = 0;
  if (m.contig === r.contig) {
    const left = Math.min(r.pos, m.pos), right = Math.max(end, endOf(m));
    tlen = right - left; if (r.pos > m.pos || (r.pos === m.pos && r.rev)) tlen = -tlen;
    const fr = (r.pos <= m.pos ? !r.rev && m.rev : r.rev && !m.rev);
    if (fr && Math.abs(tlen) < 1000) flag |= 2;
  }
  const name = Buffer.from(r.name + '\0', 'latin1');
  const tags = [];
  const md = mdTag(r);
  tags.push(Buffer.concat([Buffer.from('MDZ', 'latin1'), Buffer.from(md + '\0', 'latin1')]));
  tags.push(Buffer.from([72, 80, 67, r.hp])); // HP:C
  tags.push(Buffer.concat([Buffer.from('RGZ', 'latin1'), Buffer.from('demo\0', 'latin1')]));
  const nm = Buffer.alloc(4); nm.write('NMC', 0, 'latin1'); nm[3] = (md.match(/[ACGT]/g) || []).length; tags.push(nm);
  const lseq = r.seq.length;
  const seqB = Buffer.alloc((lseq + 1) >> 1);
  for (let i = 0; i < lseq; i++) { const c = CODE[r.seq[i]]; if (i & 1) seqB[i >> 1] |= c; else seqB[i >> 1] = c << 4; }
  const qualB = Buffer.from(r.qual);
  const cig = Buffer.alloc(4 * r.cigar.length);
  r.cigar.forEach(([op, n], i) => cig.writeUInt32LE((n << 4) | OPC[op], 4 * i));
  const fix = Buffer.alloc(36);
  const body = Buffer.concat([name, cig, seqB, qualB, ...tags]);
  fix.writeInt32LE(32 + body.length, 0);
  fix.writeInt32LE(r.contig, 4); fix.writeInt32LE(sh(r.contig, r.pos), 8);
  fix[12] = name.length; fix[13] = r.mapq; fix.writeUInt16LE(reg2bin(sh(r.contig, r.pos), sh(r.contig, end)), 14);
  fix.writeUInt16LE(r.cigar.length, 16); fix.writeUInt16LE(flag, 18); fix.writeInt32LE(lseq, 20);
  fix.writeInt32LE(m.contig, 24); fix.writeInt32LE(sh(m.contig, m.pos), 28); fix.writeInt32LE(tlen, 32);
  return { buf: Buffer.concat([fix, body]), end };
}

// ---------- BGZF ----------
class Bgzf {
  constructor() { this.out = []; this.coff = 0; this.buf = []; this.len = 0; }
  voff() { return this.coff * 65536 + this.len; }
  write(b) { this.buf.push(b); this.len += b.length; }
  flush() {
    if (!this.len) return;
    const data = Buffer.concat(this.buf);
    const cdata = zlib.deflateRawSync(data, { level: 6 });
    const h = Buffer.from([31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, 66, 67, 2, 0, 0, 0]);
    h.writeUInt16LE(cdata.length + 25, 16);
    const t = Buffer.alloc(8); t.writeUInt32LE(zlib.crc32 ? zlib.crc32(data) : crc32(data), 0); t.writeUInt32LE(data.length, 4);
    const blk = Buffer.concat([h, cdata, t]);
    this.out.push(blk); this.coff += blk.length; this.buf = []; this.len = 0;
  }
  ensure(n) { if (this.len + n > 65280) this.flush(); }
  finish() { this.flush(); this.out.push(Buffer.from('1f8b08040000000000ff0600424302001b0003000000000000000000', 'hex')); return Buffer.concat(this.out); }
}
function crc32(buf) { let c, crc = 0xFFFFFFFF; for (let n = 0; n < buf.length; n++) { c = (crc ^ buf[n]) & 0xFF; for (let k = 0; k < 8; k++) c = c & 1 ? (c >>> 1) ^ 0xEDB88320 : c >>> 1; crc = (crc >>> 8) ^ c; } return (crc ^ 0xFFFFFFFF) >>> 0; }

reads.sort((a, b) => a.contig - b.contig || a.pos - b.pos);
const bg = new Bgzf();
const text = '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:demo\tLN:' + LTOT + '\n@SQ\tSN:demo2\tLN:' + L2 + '\n@RG\tID:demo\tSM:DEMO-SAMPLE\tPL:ILLUMINA\n@PG\tID:cgv-gen-demo\tPN:cgv-gen-demo\n';
const hparts = [Buffer.from('BAM\x01', 'latin1')];
let tb = Buffer.alloc(4); tb.writeInt32LE(text.length); hparts.push(tb, Buffer.from(text, 'latin1'));
tb = Buffer.alloc(4); tb.writeInt32LE(2); hparts.push(tb);
for (const [nm, ln] of [['demo', LTOT], ['demo2', L2]]) { const a = Buffer.alloc(4); a.writeInt32LE(nm.length + 1); const b = Buffer.alloc(4); b.writeInt32LE(ln); hparts.push(a, Buffer.from(nm + '\0', 'latin1'), b); }
bg.write(Buffer.concat(hparts)); bg.flush();

const idx = [new Map(), new Map()], lin = [[], []];
for (const r of reads) {
  const { buf, end } = encode(r);
  bg.ensure(buf.length);
  const v0 = bg.voff(); bg.write(buf); const v1 = bg.voff();
  const sp = sh(r.contig, r.pos), se = sh(r.contig, end);
  const bin = reg2bin(sp, se), m = idx[r.contig];
  if (!m.has(bin)) m.set(bin, []);
  const ch = m.get(bin); if (ch.length && ch[ch.length - 1][1] === v0) ch[ch.length - 1][1] = v1; else ch.push([v0, v1]);
  for (let w = sp >> 14; w <= (se - 1) >> 14; w++) if (lin[r.contig][w] === undefined) lin[r.contig][w] = v0;
}
const bam = bg.finish();

// ---------- BAI ----------
function u64(v) { const b = Buffer.alloc(8); b.writeUInt32LE(v % 4294967296, 0); b.writeUInt32LE(Math.floor(v / 4294967296), 4); return b; }
const bai = [Buffer.from('BAI\x01', 'latin1')];
let b4 = Buffer.alloc(4); b4.writeInt32LE(2); bai.push(b4);
for (let c = 0; c < 2; c++) {
  const m = idx[c]; b4 = Buffer.alloc(4); b4.writeInt32LE(m.size); bai.push(b4);
  for (const [bin, ch] of m) { const h = Buffer.alloc(8); h.writeUInt32LE(bin, 0); h.writeInt32LE(ch.length, 4); bai.push(h); for (const [x, y] of ch) bai.push(u64(x), u64(y)); }
  const li = lin[c]; let last = 0; for (let i = 0; i < li.length; i++) { if (li[i] === undefined) li[i] = last; last = li[i]; }
  b4 = Buffer.alloc(4); b4.writeInt32LE(li.length); bai.push(b4); for (const v of li) bai.push(u64(v));
}

// ---------- FASTA ----------
function fasta(name, s) { let o = '>' + name + '\n'; for (let i = 0; i < s.length; i += 60) o += s.slice(i, i + 60) + '\n'; return o; }
const pad = 'N'.repeat(OFF);
const fa1 = fasta('demo', OFF ? pad + ref + pad : ref), fa2 = fasta('demo2', ref2);
const fai = `demo\t${LTOT}\t6\t60\t61\ndemo2\t${L2}\t${fa1.length + 7}\t60\t61\n`;

fs.writeFileSync(path.join(OUT, 'demo.bam'), bam);
fs.writeFileSync(path.join(OUT, 'demo.bam.bai'), Buffer.concat(bai));
fs.writeFileSync(path.join(OUT, 'demo.fa'), fa1 + fa2);
fs.writeFileSync(path.join(OUT, 'demo.fa.fai'), fai);
fs.writeFileSync(path.join(OUT, 'truth.json'), JSON.stringify({
  reads: reads.length,
  variants: { hetSnv: [V.hetSnv.pos, ref[V.hetSnv.pos], V.hetSnv.alt], homSnv: [V.homSnv.pos, ref[V.homSnv.pos], V.homSnv.alt], hetDel: V.hetDel, hetIns: V.hetIns, bigDel: V.bigDel },
  countIn: (c, a, b) => 0,
  overlap_12000_12001: reads.filter(r => r.contig === 0 && r.pos <= 12000 && endOf(r) > 12000).length,
  overlap_demo2_1000_1200: reads.filter(r => r.contig === 1 && r.pos < 1200 && endOf(r) > 1000).length,
}, null, 1));
console.log('reads', reads.length, 'bam', bam.length, 'bytes');
