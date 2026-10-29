// Fuzzing of every reader CGV has outside CRAM (CRAM is fuzzed in test-cram.js): BGZF, the BAM
// header and records (read through fetch, then through everything the app does with a record),
// BAI, FASTA with .fai and .gzi, FASTA indexing in the browser, the locus field and a pasted
// display profile. The rule checked: damaged input ends in a CgvError (a clear message) or in
// data, never in another exception, an endless loop or a huge allocation.
// Run: node tools/test-fuzz.js [rounds]   (deterministic: the same seed every run)
'use strict';
const fs = require('fs'), path = require('path'), zlib = require('zlib'), vm = require('vm');
const C = require('../src/core.js');
const D = path.join(__dirname, '..', 'demo');
const ROUNDS = +process.argv[2] || 300;
let fails = 0;
const ok = (cond, msg) => { console.log((cond ? 'OK   ' : 'FAIL ') + msg); if (!cond) fails++; };
let seed = 12345;
const rnd = n => { seed ^= seed << 13; seed >>>= 0; seed ^= seed >>> 17; seed ^= seed << 5; seed >>>= 0; return seed % n; };
const src = (u8, name) => new C.FileSource(new Blob([u8]), name || 'x');
const SLOW_MS = 5000; // any single case slower than this counts as a hang

// Runs f over many damaged inputs; tallies clear errors by code and fails on anything else.
async function campaign(title, n, make, run) {
  const codes = {}, bad = []; let fine = 0, slow = 0; const t0 = Date.now();
  for (let i = 0; i < n; i++) {
    const input = make(i), t = Date.now();
    try { await run(input); fine++; }
    catch (e) {
      if (e instanceof C.CgvError) codes[e.code] = (codes[e.code] || 0) + 1;
      else bad.push(`#${i}: ${e && e.constructor && e.constructor.name}: ${String(e && e.message).slice(0, 120)}`);
    }
    if (Date.now() - t > SLOW_MS) slow++;
  }
  ok(!bad.length && !slow, `${title}: ${n} cases, ${fine} read, clear errors ${JSON.stringify(codes)}, unexpected ${bad.length}, slow ${slow} (${Date.now() - t0} ms)`);
  for (const b of bad.slice(0, 5)) console.log('       ' + b);
}

// ---------- building BGZF with valid checksums, so damage gets past the CRC into the parsers ----------
function bgzfBlock(data) {
  const c = zlib.deflateRawSync(data, { level: 1 });
  const h = Buffer.from([31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, 66, 67, 2, 0, 0, 0]);
  h.writeUInt16LE(h.length + c.length + 8 - 1, 16);
  const tl = Buffer.alloc(8); tl.writeUInt32LE(zlib.crc32(data), 0); tl.writeUInt32LE(data.length, 4);
  return Buffer.concat([h, c, tl]);
}
// Original blocks: [{ coff, data }].
function bgzfSplit(buf) {
  return C.parseBgzfBlocks(buf, 0).map(b => ({ coff: b.coff, data: zlib.inflateRawSync(buf.subarray(b.cstart, b.cstart + b.clen)) }));
}
// Re-compresses (possibly changed) blocks; returns the file and old → new block offsets.
function bgzfJoin(blocks, oldEnd) {
  const parts = [], map = new Map(); let off = 0;
  for (const b of blocks) { map.set(b.coff, off); const x = bgzfBlock(b.data); parts.push(x); off += x.length; }
  map.set(oldEnd, off);
  return { file: Buffer.concat(parts), map };
}
// The BAI with every virtual offset moved to the re-compressed file.
function baiRemap(bai, map) {
  const out = Buffer.from(bai); let p = 8;
  const fix = o => { const v = out.readBigUInt64LE(o), co = Number(v >> 16n), u = v & 0xffffn; const n = map.has(co) ? map.get(co) : co; out.writeBigUInt64LE((BigInt(n) << 16n) | u, o); };
  const nref = out.readInt32LE(4);
  for (let r = 0; r < nref; r++) {
    const nbin = out.readInt32LE(p); p += 4;
    for (let b = 0; b < nbin; b++) {
      const bin = out.readUInt32LE(p), nch = out.readInt32LE(p + 4); p += 8;
      for (let c = 0; c < nch; c++) { if (bin !== 37450 || c === 0) { fix(p); fix(p + 8); } p += 16; }
    }
    const ni = out.readInt32LE(p); p += 4;
    for (let i = 0; i < ni; i++) { fix(p); p += 8; }
  }
  return out;
}

// Everything the app does with a record it got: none of it may throw on a record fetch accepted.
function useRecords(reads) {
  const ref = new Uint8Array(50000).fill(65);
  for (const r of reads) {
    C.recTags(r); C.cigarString(r); C.clipLens(r); C.hardClipLens(r);
    C.computeFeatures(r, p => ref[p] || 0);
    C.fillRefFromMD(r, ref, 0);
    C.baseAtRefPos(r, r.pos + 3); C.queryAt(r, r.pos + 3);
    C.pairOrientation(r);
  }
  C.computeCoverage(reads, 0, 50000);
  C.insertSizeBounds(reads);
  const units = C.pairUnits(reads, 10000, true);
  new C.ReadLayout(null, 2).update(units, () => '', null, u => u.key);
  new C.ReadLayout(null, 2).update(reads, r => r.flag & 16 ? 'r' : 'f');
}
async function readAllWays(bam, bai) {
  const f = await new C.BamFile(src(bam, 'f.bam'), src(bai, 'f.bam.bai')).open();
  for (let r = 0; r < Math.min(f.refs.length, 3); r++) {
    const a = await f.fetch(r, 0, 40000);
    useRecords(a.reads);
    await f.fetch(r, 9000, 9100, { downsample: { max: 5, window: 50 }, filters: { minMapq: 0, hideDuplicates: true, hideSecondary: true, hideSupplementary: false, hideQcFail: true } });
    await f.fetch(r, 0, 40000, { depthOnly: true });
  }
}

(async () => {
  const bamOrig = fs.readFileSync(path.join(D, 'demo.bam')), baiOrig = fs.readFileSync(path.join(D, 'demo.bam.bai'));
  const blocks = bgzfSplit(bamOrig);
  const joined = bgzfJoin(blocks, bamOrig.length);
  const bai0 = baiRemap(baiOrig, joined.map);
  // the re-compressed file must read exactly like the original before it is damaged
  {
    const a = await new C.BamFile(src(bamOrig), src(baiOrig)).open(), b = await new C.BamFile(src(joined.file), src(bai0)).open();
    const x = (await a.fetch(0, 0, 40000)).reads, y = (await b.fetch(0, 0, 40000)).reads;
    ok(x.length === y.length && x.length > 1000 && x.every((r, i) => r.name === y[i].name && r.pos === y[i].pos), `re-compressed demo BAM reads the same: ${x.length} reads`);
  }

  // 1. the raw file: bytes changed, cut, or shifted (the BGZF checksums catch most of it)
  await campaign('BAM file, raw damage (BGZF checks on)', ROUNDS, i => {
    const b = Buffer.from(bamOrig), k = i % 3;
    if (k === 0) { for (let j = 1 + rnd(8); j--;) b[rnd(b.length)] ^= 1 << rnd(8); return b; }
    if (k === 1) return b.subarray(0, rnd(b.length));
    const s = rnd(b.length), e = Math.min(b.length, s + rnd(4000)); return Buffer.concat([b.subarray(0, s), b.subarray(e)]);
  }, b => readAllWays(b, baiOrig));

  // 2. damage inside the decompressed data, re-compressed with valid checksums: reaches the parsers
  const recBlocks = blocks.map((b, i) => i).filter(i => blocks[i].data.length > 1000);
  const targeted = (data, i) => {
    // a field of some record set to an extreme value: block_size, l_read_name, n_cigar_op, l_seq, a CIGAR op
    const dv = new DataView(data.buffer, data.byteOffset, data.byteLength);
    let p = 0; const starts = [];
    while (p + 36 <= data.length) { const bs = dv.getInt32(p, true); if (bs < 32 || p + 4 + bs > data.length) break; starts.push(p); p += 4 + bs; }
    if (!starts.length) { data[rnd(data.length)] ^= 255; return; }
    const s = starts[rnd(starts.length)], v = [0, -1, 1, 31, 32, 0x7fffffff, 0x10000000, 65535, 255][rnd(9)];
    const field = [0, 12, 16, 20, 36][i % 5]; // block_size, l_read_name.., n_cigar_op.., l_seq, first CIGAR op
    if (field === 12) data[s + 12] = v & 255;
    else if (field === 16) dv.setUint16(s + 16, v & 0xffff, true);
    else if (s + field + 4 <= data.length) dv.setInt32(s + field, v | 0, true);
  };
  await campaign('BAM records, damage behind valid checksums', ROUNDS, i => {
    const copy = blocks.map(b => ({ coff: b.coff, data: Buffer.from(b.data) }));
    const blk = copy[recBlocks[rnd(recBlocks.length)]].data;
    if (i % 2) targeted(blk, i >> 1); else for (let j = 1 + rnd(20); j--;) blk[rnd(blk.length)] = rnd(256);
    const j = bgzfJoin(copy, bamOrig.length); return [j.file, baiRemap(baiOrig, j.map)];
  }, ([b, x]) => readAllWays(b, x));

  // 3. the header: magic, l_text, n_ref, names and lengths
  await campaign('BAM header, damage behind valid checksums', ROUNDS, i => {
    const copy = blocks.map(b => ({ coff: b.coff, data: Buffer.from(b.data) }));
    const h = copy[0].data, dv = new DataView(h.buffer, h.byteOffset, h.byteLength);
    const lText = dv.getInt32(4, true);
    if (i % 3 === 0) dv.setInt32(4, [-1, 0x7fffffff, lText + 1, lText - 1, 1 << 27][rnd(5)], true);
    else if (i % 3 === 1) dv.setInt32(8 + lText, [-5, 0x7fffffff, 3, 1 << 25][rnd(4)], true);
    else for (let j = 1 + rnd(6); j--;) h[rnd(Math.min(h.length, 8 + lText + 80))] = rnd(256);
    const j = bgzfJoin(copy, bamOrig.length); return [j.file, baiRemap(baiOrig, j.map)];
  }, ([b, x]) => readAllWays(b, x));

  // 4. the index
  await campaign('BAI index, damaged', ROUNDS, i => {
    const b = Buffer.from(baiOrig), k = i % 3;
    if (k === 0) { for (let j = 1 + rnd(10); j--;) b[8 + rnd(b.length - 8)] = rnd(256); return b; }
    if (k === 1) return b.subarray(0, rnd(b.length));
    const dv = new DataView(b.buffer, b.byteOffset, b.byteLength), at = 4 + 4 * rnd(Math.floor((b.length - 4) / 4));
    dv.setInt32(at, [-1, 0x7fffffff, 1 << 24, 0][rnd(4)], true); return b;
  }, b => readAllWays(bamOrig, b));

  // 5. BGZF on its own: odd XLEN, BSIZE, ISIZE values
  await campaign('BGZF block headers', ROUNDS, () => {
    const b = Buffer.from(bamOrig.subarray(0, 200000));
    const pos = C.parseBgzfBlocks(b, 0).map(x => x.coff), at = pos[rnd(pos.length)];
    const o = [10, 11, 16, 17][rnd(4)];
    if (rnd(2)) b[at + o] = rnd(256); else { const bl = C.parseBgzfBlocks(b, 0).find(x => x.coff === at); b.writeUInt32LE([0, 65537, 0xffffffff, 1][rnd(4)], bl.cstart + bl.clen + 4); }
    return b;
  }, async b => { for (const bl of C.parseBgzfBlocks(b, 0)) await C.inflateBgzf(b, bl); });

  // a deflate "bomb": a tiny block that inflates far past 64 KB must stop at the limit
  {
    const big = zlib.deflateRawSync(Buffer.alloc(50 << 20)), u8 = new Uint8Array(big);
    let code = '';
    try { await C.inflateBgzf(u8, { coff: 0, cstart: 0, clen: u8.length, crc: 0, isize: 65536 }); } catch (e) { code = e.code; }
    ok(code === 'bgzfCrc', `deflate bomb (${big.length} bytes → 50 MB) is stopped at the block size: ${code}`);
  }

  // 6. FASTA with .fai and .gzi
  const fa = fs.readFileSync(path.join(D, 'demo.fa')), fai = fs.readFileSync(path.join(D, 'demo.fa.fai'), 'utf8');
  const fetchSome = async f => { for (const n of f.idx.keys()) { await f.fetch(n, 0, 100); await f.fetch(n, 39000, 41000); await f.fetch(n, -5, 3); } };
  await campaign('FASTA index (.fai), damaged', ROUNDS, i => {
    const lines = fai.trim().split('\n').map(l => l.split('\t'));
    const L = lines[rnd(lines.length)], k = 1 + rnd(4);
    L[k] = ['-1', '0', '1e9', 'abc', '99999999999999999999', '', '3.5', String(+L[k] + 1)][rnd(8)];
    if (i % 7 === 0) L.length = 3;
    return lines.map(l => l.join('\t')).join('\n') + '\n';
  }, async t => { const f = await new C.IndexedFasta(src(fa, 'demo.fa'), src(Buffer.from(t), 'demo.fa.fai')).open(); await fetchSome(f); });

  // bgzip FASTA: blocks of 64 KB made here, with a .gzi as bgzip writes it
  const faBlocks = []; for (let o = 0; o < fa.length; o += 20000) faBlocks.push(fa.subarray(o, o + 20000));
  const fgz = Buffer.concat(faBlocks.map(bgzfBlock).concat([bgzfBlock(Buffer.alloc(0))]));
  const gzi = Buffer.alloc(8 + (faBlocks.length - 1) * 16); gzi.writeBigUInt64LE(BigInt(faBlocks.length - 1), 0);
  { let co = 0, uo = 0; faBlocks.forEach((b, i) => { if (i) { gzi.writeBigUInt64LE(BigInt(co), 8 + (i - 1) * 16); gzi.writeBigUInt64LE(BigInt(uo), 16 + (i - 1) * 16); } co += bgzfBlock(b).length; uo += b.length; }); }
  {
    const f = await new C.IndexedFasta(src(fgz, 'demo.fa.gz'), src(Buffer.from(fai), 'demo.fa.fai'), src(gzi, 'demo.fa.gz.gzi')).open();
    const plain = await new C.IndexedFasta(src(fa, 'demo.fa'), src(Buffer.from(fai), 'demo.fa.fai')).open();
    const a = await f.fetch('demo', 19990, 20110), b = await plain.fetch('demo', 19990, 20110);
    ok(Buffer.compare(Buffer.from(a), Buffer.from(b)) === 0, 'bgzip FASTA with .gzi reads like the plain one');
  }
  await campaign('FASTA .gzi, damaged', ROUNDS, i => {
    const b = Buffer.from(gzi);
    if (i % 3 === 0) b.writeBigUInt64LE(BigInt([0, 1, 1000, 2 ** 40][rnd(4)]), 0);
    else if (i % 3 === 1) return b.subarray(0, rnd(b.length));
    else for (let j = 1 + rnd(6); j--;) b[8 + rnd(b.length - 8)] = rnd(256);
    return b;
  }, async g => { const f = await new C.IndexedFasta(src(fgz, 'demo.fa.gz'), src(Buffer.from(fai), 'demo.fa.fai'), src(g, 'demo.fa.gz.gzi')).open(); await fetchSome(f); });
  await campaign('bgzip FASTA, damaged blocks', ROUNDS, () => {
    const b = Buffer.from(fgz); for (let j = 1 + rnd(5); j--;) b[rnd(b.length)] ^= 1 << rnd(8); return b;
  }, async b => { const f = await new C.IndexedFasta(src(b, 'demo.fa.gz'), src(Buffer.from(fai), 'demo.fa.fai'), src(gzi, 'demo.fa.gz.gzi')).open(); await fetchSome(f); });

  // 7. FASTA without an index: indexed in the browser (plain, bgzip and ordinary gzip)
  await campaign('FASTA indexed in the browser, damaged', ROUNDS, i => {
    const b = Buffer.from(fa.subarray(0, 30000 + rnd(30000)));
    for (let j = 1 + rnd(10); j--;) b[rnd(b.length)] = [10, 13, 62, 32, 0, 65][rnd(6)];
    const k = i % 3;
    return k === 0 ? b : k === 1 ? Buffer.concat([bgzfBlock(b.subarray(0, 20000)), bgzfBlock(b.subarray(20000)), bgzfBlock(Buffer.alloc(0))]) : zlib.gzipSync(b);
  }, async b => { const f = await new C.IndexedFasta(src(b, 'x.fa')).openUnindexed(b.length); await fetchSome(f); });

  // 8. the locus field and a pasted profile: the functions are taken from app.js as they are
  const app = fs.readFileSync(path.join(__dirname, '..', 'src', 'app.js'), 'utf8');
  const grab = (start, end) => { const a = app.indexOf(start), b = app.indexOf(end, a); if (a < 0 || b < 0) throw new Error('app.js changed: ' + start); return app.slice(a, b); };
  const Gn = require('../src/genes.js'); Gn.parseGenes(fs.readFileSync(path.join(__dirname, '..', 'data', 'mane.grch38.tsv'), 'utf8'));
  const ctx = { S: { contigs: [{ name: 'chr1', len: 248956422 }, { name: 'chrM', len: 16569 }, { name: 'chr17', len: 83257441 }] }, t: k => k, result: null,
    GENES: Gn.GENES, findGene: Gn.findGene, genesActive: () => true };
  ctx.contigId = n => ctx.S.contigs.findIndex(c => c.name === n || c.name === 'chr' + n);
  vm.createContext(ctx);
  vm.runInContext(grab('const DEFAULT_PROFILE', 'let P = DEFAULT_PROFILE;') + grab('function parseLocus', '// ======================= Input') +
    ';this.mergeProfile = mergeProfile; this.parseLocus = parseLocus; this.DEFAULT_PROFILE = DEFAULT_PROFILE;', ctx);
  {
    const loci = ['chr1:1-100', 'chr1:100', 'chrM:20000-30000', 'chr1:500-100', 'chr1:99999999999999999999-1', 'chr1:0-0', 'chr1:-5', 'chr1:1,000,000-2,000,000', '1:5', 'x'.repeat(5000), 'chr1:' + '9'.repeat(400), ':::', 'chr1:1-2-3', '', 'chrM:1.000-5', 'BRCA1', 'brca1', 'NM_007294.4', 'ENST00000357654', 'BRCA1:5', 'NOTAGENE', 'BRCA1 BRCA2'];
    const badLoc = [];
    for (const l of loci) {
      const r = ctx.parseLocus(l);
      if (r && !r.error) { const L = ctx.S.contigs[r.id].len; if (!(Number.isSafeInteger(r.s) && Number.isSafeInteger(r.e) && r.e > r.s && r.e <= L + 51 && r.s >= -50)) badLoc.push(l.slice(0, 40) + ' → ' + JSON.stringify(r)); }
    }
    const g = ctx.parseLocus('BRCA1');
    if (!g || g.error || ctx.S.contigs[g.id].name !== 'chr17' || g.s > 43044294 || g.e < 43125364) badLoc.push('BRCA1 is not found where MANE puts it: ' + JSON.stringify(g));
    const r1 = ctx.parseLocus('chr1:500-100'), r2 = ctx.parseLocus('chr1:100-500');
    if (JSON.stringify(r1) !== JSON.stringify(r2) || r2.s !== 99 || r2.e !== 500) badLoc.push('500-100 is not the same as 100-500: ' + JSON.stringify([r1, r2]));
    ok(!badLoc.length, `locus field: ${loci.length} odd inputs, every result inside its contig` + (badLoc.length ? ': ' + badLoc.join('; ') : ''));
    const evil = [
      { readHeight: 1e308, maxLoadMB: -1, maxReads: NaN, coverage: { height: 'x', alleleThreshold: 7 }, downsample: { maxReads: -3, window: 0, extra: 1 } },
      JSON.parse('{"__proto__": {"polluted": 1}, "constructor": {"prototype": {"polluted": 1}}, "filters": {"__proto__": {"polluted": 1}}}'),
      { colorBy: 'evil', language: '<img>', colorTag: '"><script>', mapqAlpha: { minAlpha: '1); background: url(x)' } },
      { mapqAlpha: null, filters: [1, 2], pairs: 'x', downsample: { enabled: 'yes' } },
    ];
    const D0 = ctx.DEFAULT_PROFILE, probs = [];
    for (const e of evil) {
      const p = ctx.mergeProfile(D0, e);
      const walk = (x, d, at) => { for (const k of Object.keys(x)) { if (!(k in d)) probs.push(at + k + ' kept'); else if (typeof x[k] !== typeof d[k]) probs.push(at + k + ' type'); else if (typeof x[k] === 'number' && !Number.isFinite(x[k])) probs.push(at + k + ' not finite'); else if (d[k] && typeof d[k] === 'object') walk(x[k], d[k], at + k + '.'); } };
      walk(p, D0, '');
      if (p.readHeight > 24 || p.maxLoadMB < 1 || p.coverage.alleleThreshold > 1 || p.downsample.window < 1 || p.colorBy === 'evil' || p.language !== 'en' || p.colorTag !== 'HP') probs.push('a value outside its range: ' + JSON.stringify(p).slice(0, 200));
    }
    ok(!probs.length && ({}).polluted === undefined && Object.prototype.polluted === undefined, `pasted profile: ${evil.length} hostile profiles, all values back in range, no prototype pollution` + (probs.length ? ': ' + probs.join('; ') : ''));
  }

  console.log(fails ? `\n${fails} check(s) failed` : '\nAll fuzz checks passed');
  process.exit(fails ? 1 : 0);
})().catch(e => { console.error(e); process.exit(1); });
