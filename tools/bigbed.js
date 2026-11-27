// Writes bigBed files (UCSC bbiFile.h / bPlusTree.h / cirTree.h, bigBed version 4) for LGV's gene
// tracks: BED records sorted by contig and start, compressed in blocks (zlib), an R tree over the
// blocks, the contig B+ tree, and optional extra B+ tree indexes on BED fields (e.g. the name), so a
// reader can fetch only what a view or a name needs. No zoom levels (gene tracks do not need
// coverage summaries). Checked against an independent reader (libBigWig) by tools/test-bigbed.js.
'use strict';
const zlib = require('zlib');

const BB_MAGIC = 0x8789F2EB, BPT_MAGIC = 0x78CA8C91, CIR_MAGIC = 0x2468ACE0;

class Out {
  constructor() { this.parts = []; this.len = 0; }
  push(b) { this.parts.push(b); this.len += b.length; }
  u8(v) { const b = Buffer.alloc(1); b.writeUInt8(v); this.push(b); }
  u16(v) { const b = Buffer.alloc(2); b.writeUInt16LE(v); this.push(b); }
  u32(v) { const b = Buffer.alloc(4); b.writeUInt32LE(v >>> 0); this.push(b); }
  u64(v) { const b = Buffer.alloc(8); b.writeBigUInt64LE(BigInt(v)); this.push(b); }
  f64(v) { const b = Buffer.alloc(8); b.writeDoubleLE(v); this.push(b); }
  fixed(s, n) { const b = Buffer.alloc(n); b.write(s, 'latin1'); this.push(b); }
  zeros(n) { this.push(Buffer.alloc(n)); }
  buf() { return Buffer.concat(this.parts); }
}

// B+ tree (bPlusTree.h) over sorted [key, value Buffer] items, written at `base`; returns a Buffer.
function bPlusTree(items, keySize, valSize, blockSize, base) {
  const o = new Out(), n = items.length;
  o.u32(BPT_MAGIC); o.u32(blockSize); o.u32(keySize); o.u32(valSize); o.u64(n); o.u64(0);
  if (!n) { o.u8(1); o.u8(0); o.u16(0); o.zeros(blockSize * (keySize + valSize)); return o.buf(); }
  // levels: leaves hold up to blockSize items; each level above has one item per node below
  let levels = 1; for (let k = blockSize; k < n; k *= blockSize) levels++;
  const nodeBytes = (leaf) => 4 + blockSize * (keySize + (leaf ? valSize : 8));
  // node counts per level, top (root) first
  const counts = []; for (let l = levels - 1; l >= 0; l--) counts.push(Math.ceil(n / blockSize ** (l + 1)));
  let off = base + 32; const levelOff = [];
  counts.forEach((c, i) => { levelOff.push(off); off += c * nodeBytes(i === levels - 1); });
  for (let li = 0; li < levels; li++) {
    const leaf = li === levels - 1, span = blockSize ** (levels - 1 - li), slot = blockSize * span;
    for (let start = 0; start < n; start += slot) {
      const kids = []; for (let i = start; i < Math.min(n, start + slot); i += span) kids.push(i);
      o.u8(leaf ? 1 : 0); o.u8(0); o.u16(kids.length);
      kids.forEach((i, j) => {
        o.fixed(items[i][0], keySize);
        if (leaf) o.push(items[i][1]);
        else o.u64(levelOff[li + 1] + (start / span + j) * nodeBytes(li + 1 === levels - 1));
      });
      o.zeros((blockSize - kids.length) * (keySize + (leaf ? valSize : 8)));
    }
  }
  return o.buf();
}

// R tree (cirTree.h) over leaf items {c0, s0, c1, s1, off, size}, written at `base`.
function cirTree(items, blockSize, itemsPerSlot, endFileOffset, base) {
  const o = new Out(), n = items.length;
  const b = n ? items.reduce((a, x) => ({ c0: Math.min(a.c0, x.c0), s0: a.c0 < x.c0 ? a.s0 : x.c0 < a.c0 ? x.s0 : Math.min(a.s0, x.s0), c1: Math.max(a.c1, x.c1), s1: a.c1 > x.c1 ? a.s1 : x.c1 > a.c1 ? x.s1 : Math.max(a.s1, x.s1) }), { c0: Infinity, s0: Infinity, c1: -1, s1: -1 }) : { c0: 0, s0: 0, c1: 0, s1: 0 };
  o.u32(CIR_MAGIC); o.u32(blockSize); o.u64(n); o.u32(b.c0); o.u32(b.s0); o.u32(b.c1); o.u32(b.s1); o.u64(endFileOffset); o.u32(itemsPerSlot); o.u32(0);
  // build levels bottom-up: level 0 = leaves (groups of items), above = groups of nodes
  const merge = g => ({ c0: g[0].c0, s0: g[0].s0, c1: g[g.length - 1].c1, s1: g.reduce((m, x) => (x.c1 > m.c1 || (x.c1 === m.c1 && x.s1 > m.s1)) ? x : m, g[0]).s1 });
  let level = []; for (let i = 0; i < n; i += blockSize) level.push(items.slice(i, i + blockSize));
  if (!level.length) level = [[]];
  const levels = [level];
  while (levels[levels.length - 1].length > 1) { const prev = levels[levels.length - 1], up = []; for (let i = 0; i < prev.length; i += blockSize) up.push(prev.slice(i, i + blockSize)); levels.push(up); }
  // write top-down; node size is fixed (blockSize slots) so offsets are computable
  const leafBytes = 4 + blockSize * 32, innerBytes = 4 + blockSize * 24;
  let off = base + 48; const levelOff = [];
  for (let l = levels.length - 1; l >= 0; l--) { levelOff[l] = off; off += levels[l].length * (l === 0 ? leafBytes : innerBytes); }
  // compute bounds of each node bottom-up
  levels[0].forEach(g => { g.__b = g.length ? merge(g) : { c0: 0, s0: 0, c1: 0, s1: 0 }; });
  for (let l = 1; l < levels.length; l++) levels[l].forEach(g => { g.__b = merge(g.map(x => x.__b)); });
  for (let l = levels.length - 1; l >= 0; l--) {
    let idx = 0;
    for (const g of levels[l]) {
      o.u8(l === 0 ? 1 : 0); o.u8(0); o.u16(g.length);
      for (const x of g) {
        if (l === 0) { o.u32(x.c0); o.u32(x.s0); o.u32(x.c1); o.u32(x.s1); o.u64(x.off); o.u64(x.size); }
        else { const bb = x.__b; o.u32(bb.c0); o.u32(bb.s0); o.u32(bb.c1); o.u32(bb.s1); o.u64(levelOff[l - 1] + idx * (l - 1 === 0 ? leafBytes : innerBytes)); idx++; }
      }
      o.zeros((blockSize - g.length) * (l === 0 ? 32 : 24));
    }
  }
  return o.buf();
}

// records: [{ chrom, start, end, rest }] (rest: the BED columns after the 3rd, tab-separated);
// chroms: [[name, size]]; opts: { autoSql, fieldCount, definedFieldCount, itemsPerSlot, blockSize,
// extraIndex: [{ field (0-based column), key: rec => string }] }
function writeBigBed(records, chroms, opts) {
  const itemsPerSlot = opts.itemsPerSlot || 512, blockSize = opts.blockSize || 256, level = opts.level ?? 9;
  const chromIds = new Map(chroms.map(([name], i) => [name, i]));
  const recs = records; // sorted in place (no copies: a full gene set is a few hundred thousand records)
  for (const r of recs) { r.cid = chromIds.get(r.chrom); if (r.cid === undefined) throw new Error('unknown contig ' + r.chrom); }
  recs.sort((a, b) => a.cid - b.cid || a.start - b.start || a.end - b.end);
  const out = new Out();
  // header placeholder (64), no zoom headers
  out.zeros(64);
  const autoSqlOffset = out.len; out.push(Buffer.from(opts.autoSql + '\0', 'latin1'));
  const totalSummaryOffset = out.len; out.zeros(40);
  const extensionOffset = out.len; out.zeros(64);
  const extra = opts.extraIndex || [];
  const extraListOffset = out.len; const extraListPos = [];
  for (let k = 0; k < extra.length; k++) { extraListPos.push(out.len); out.zeros(20); } // type, field count, offset, reserved; one field id
  // contig tree
  const chromTreeOffset = out.len;
  const keySize = Math.max(1, ...chroms.map(c => c[0].length));
  out.push(bPlusTree(chroms.map(([name, size], i) => { const v = Buffer.alloc(8); v.writeUInt32LE(i, 0); v.writeUInt32LE(size, 4); return [name, v]; }).sort((a, b) => a[0] < b[0] ? -1 : a[0] > b[0] ? 1 : 0), keySize, 8, Math.min(256, chroms.length || 1), chromTreeOffset));
  // data: record count, then blocks of up to itemsPerSlot records, one contig per block
  const fullDataOffset = out.len; out.u64(recs.length);
  const leaves = [], blockOf = new Array(recs.length); let maxUnc = 0, i = 0;
  let cov = 0, minV = Infinity, maxV = -Infinity, sum = 0, sq = 0;
  while (i < recs.length) {
    const j0 = i, cid = recs[i].cid, parts = [];
    while (i < recs.length && i - j0 < itemsPerSlot && recs[i].cid === cid) {
      const r = recs[i], b = Buffer.alloc(12); b.writeUInt32LE(cid, 0); b.writeUInt32LE(r.start, 4); b.writeUInt32LE(r.end, 8);
      parts.push(b, Buffer.from(r.rest + '\0', 'latin1')); blockOf[i] = leaves.length; i++;
    }
    const raw = Buffer.concat(parts); maxUnc = Math.max(maxUnc, raw.length);
    const z = zlib.deflateSync(raw, { level, memLevel: 9 });
    const g = recs.slice(j0, i);
    leaves.push({ c0: cid, s0: g[0].start, c1: cid, s1: Math.max(...g.map(r => r.end)), off: out.len, size: z.length });
    out.push(z);
  }
  // summary over covered bases (validCount), values all 1 here (BED has no values)
  for (const r of recs) { cov += r.end - r.start; }
  if (recs.length) { minV = 1; maxV = 1; sum = cov; sq = cov; } else { minV = maxV = 0; }
  const fullIndexOffset = out.len;
  out.push(cirTree(leaves, blockSize, itemsPerSlot, fullIndexOffset, fullIndexOffset));
  // extra indexes: B+ tree keyed by the field's value, value = (block offset, block size)
  const extraOffsets = [];
  for (const e of extra) {
    // one item per (key, block): a name in several blocks (the transcripts of a gene) is found in all
    const seen = new Set(), items = [];
    recs.forEach((r, k) => { const key = e.key(r); if (!key) return; const id = key + '\t' + blockOf[k]; if (seen.has(id)) return; seen.add(id); items.push([key, leaves[blockOf[k]]]); });
    items.sort((a, b) => a[0] < b[0] ? -1 : a[0] > b[0] ? 1 : a[1].off - b[1].off);
    const ks = Math.max(1, ...items.map(x => Buffer.byteLength(x[0], 'latin1')));
    extraOffsets.push(out.len);
    out.push(bPlusTree(items.map(([k, lf]) => { const v = Buffer.alloc(16); v.writeBigUInt64LE(BigInt(lf.off), 0); v.writeBigUInt64LE(BigInt(lf.size), 8); return [k, v]; }), ks, 16, 256, out.len));
  }
  out.u32(BB_MAGIC); // trailing magic, as UCSC writes
  const buf = out.buf();
  // header
  const h = buf; let p = 0;
  const w32 = v => { h.writeUInt32LE(v >>> 0, p); p += 4; }, w16 = v => { h.writeUInt16LE(v, p); p += 2; }, w64 = v => { h.writeBigUInt64LE(BigInt(v), p); p += 8; };
  w32(BB_MAGIC); w16(4); w16(0); w64(chromTreeOffset); w64(fullDataOffset); w64(fullIndexOffset);
  w16(opts.fieldCount); w16(opts.definedFieldCount); w64(autoSqlOffset); w64(totalSummaryOffset); w32(maxUnc); w64(extensionOffset);
  // total summary
  p = totalSummaryOffset; w64(cov); h.writeDoubleLE(minV, p); p += 8; h.writeDoubleLE(maxV, p); p += 8; h.writeDoubleLE(sum, p); p += 8; h.writeDoubleLE(sq, p); p += 8;
  // extension header: size, extra index count, extra index list offset
  p = extensionOffset; w16(64); w16(extra.length); w64(extra.length ? extraListOffset : 0);
  extra.forEach((e, k) => { p = extraListPos[k]; w16(0); w16(1); w64(extraOffsets[k]); w32(0); w16(e.field); w16(0); });
  return buf;
}
module.exports = { writeBigBed };
