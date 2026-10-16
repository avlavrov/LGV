// Core checks on demo data: read counts, MD-rebuilt reference vs FASTA, known variants.
'use strict';
const fs = require('fs'), path = require('path');
const C = require('../src/core.js');
const D = path.join(__dirname, '..', 'demo');
const truth = JSON.parse(fs.readFileSync(path.join(D, 'truth.json')));
const blob = f => new Blob([fs.readFileSync(path.join(D, f))]);
let fails = 0;
const ok = (cond, msg) => { console.log((cond ? 'OK   ' : 'FAIL ') + msg); if (!cond) fails++; };

(async () => {
  const bam = await new C.BamFile(new C.FileSource(blob('demo.bam')), new C.FileSource(blob('demo.bam.bai'))).open();
  ok(bam.refs.length === 2 && bam.refs[0].name === 'demo' && bam.refs[0].len === 40000, 'header: ' + JSON.stringify(bam.refs) + ' samples ' + bam.samples);
  const fa = await new C.IndexedFasta(new C.FileSource(blob('demo.fa')), new C.FileSource(blob('demo.fa.fai'))).open();

  let t = Date.now();
  const all = await bam.fetch(0, 0, 40000);
  const all2 = await bam.fetch(1, 0, 2000);
  ok(all.length + all2.length === truth.reads, `all reads: ${all.length}+${all2.length} = ${truth.reads} (${Date.now() - t} ms)`);
  ok(new Set(all.map(r => r.name + (r.flag & 64))).size === all.length, 'no duplicates after merging chunks');

  const at = await bam.fetch(0, 12000, 12001);
  ok(at.length === truth.overlap_12000_12001, `overlapping 12000: ${at.length} = ${truth.overlap_12000_12001}`);
  const at2 = await bam.fetch(1, 1000, 1200);
  ok(at2.length === truth.overlap_demo2_1000_1200, `demo2:1000-1200: ${at2.length} = ${truth.overlap_demo2_1000_1200}`);

  // MD-rebuilt reference matches FASTA
  const ws = 0, refArr = new Uint8Array(40000);
  for (const r of all) C.fillRefFromMD(r, refArr, ws);
  const faSeq = await fa.fetch('demo', 0, 40000);
  let filled = 0, bad = 0; for (let i = 0; i < 40000; i++) if (refArr[i]) { filled++; if (refArr[i] !== faSeq[i]) bad++; }
  ok(bad === 0 && filled > 30000, `reference from MD: ${filled} filled, ${bad} differ from FASTA`);
  const fa2 = await fa.fetch('demo2', 995, 1005);
  ok(fa2.length === 10, 'FASTA, second contig: ' + String.fromCharCode(...fa2));

  // variants
  const [sp, sref, salt] = truth.variants.hetSnv;
  const cov = C.computeCoverage(all.filter(r => !(r.flag & 1024)), 0, 40000);
  const idx = { A: 0, C: 1, G: 2, T: 3 };
  const c12 = cov.slice(sp * 6, sp * 6 + 6);
  ok(c12[idx[salt]] > 5 && c12[idx[sref]] > 5, `heterozygous SNV ${sref}>${salt} @${sp}: ${Array.from(c12)}`);
  const [hp, , halt] = truth.variants.homSnv; const c15 = cov.slice(hp * 6, hp * 6 + 6);
  ok(c15[idx[halt]] / (c15[0] + c15[1] + c15[2] + c15[3]) > 0.9, `homozygous SNV @${hp}: ${Array.from(c15)}`);
  ok(cov[18001 * 6 + 5] > 5, '3 bp deletion: ' + cov[18001 * 6 + 5]);
  const ins = all.filter(r => C.cigarString(r).includes('3I')).length;
  ok(ins > 5, 'reads with 3 bp insertion: ' + ins);
  const refAt = p => faSeq[p];
  const mism = all.filter(r => r.pos <= sp && r.end > sp).flatMap(r => C.computeFeatures(r, refAt)).filter(f => f.t === 0 && f.pos === sp).length;
  ok(mism > 5, 'mismatches at SNV: ' + mism);
  const mid = cov[26500 * 6] + cov[26500 * 6 + 1] + cov[26500 * 6 + 2] + cov[26500 * 6 + 3];
  ok(mid > 5 && mid < 25, 'coverage inside large het deletion is ~half: ' + mid);
  const orient = {}; for (const r of all) { const o = C.pairOrientation(r); orient[o] = (orient[o] || 0) + 1; }
  console.log('pair orientations', orient, 'insert bounds', C.insertSizeBounds(all));
  ok(orient.other >= 24, 'chimeric pairs');
  const rows = C.packRows(all.slice(0, 3000), 2);
  let overlap = 0; for (const row of rows) for (let i = 1; i < row.length; i++) if (row[i].pos < row[i - 1].end + 2) overlap++;
  ok(overlap === 0, `layout: ${rows.length} rows, ${overlap} overlaps`);
  // bgzip-compressed FASTA with .gzi gives the same bases as plain FASTA
  {
    const zlib = require('zlib'), raw = fs.readFileSync(path.join(D, 'demo.fa'));
    const blocks = [], gzi = []; let coff = 0;
    for (let u = 0; u < raw.length; u += 1000) { // small blocks to cross block borders often
      const data = raw.subarray(u, u + 1000), c = zlib.deflateRawSync(data);
      const h = Buffer.from([31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, 66, 67, 2, 0, 0, 0]); h.writeUInt16LE(c.length + 25, 16);
      const tl = Buffer.alloc(8); tl.writeUInt32LE(zlib.crc32(data), 0); tl.writeUInt32LE(data.length, 4);
      if (u) gzi.push([coff, u]);
      const blk = Buffer.concat([h, c, tl]); blocks.push(blk); coff += blk.length;
    }
    const g = Buffer.alloc(8 + gzi.length * 16); g.writeUInt32LE(gzi.length, 0);
    gzi.forEach(([c, u], i) => { g.writeUInt32LE(c, 8 + i * 16); g.writeUInt32LE(u, 16 + i * 16); });
    const fz = await new C.IndexedFasta(new C.FileSource(new Blob([Buffer.concat(blocks)]), 'demo.fa.gz'), new C.FileSource(blob('demo.fa.fai')), new C.FileSource(new Blob([g]))).open();
    let same = true;
    for (const [b0, b1] of [[0, 10], [995, 2050], [12345, 13999], [39990, 40000]]) { const x = await fz.fetch('demo', b0, b1), y = await fa.fetch('demo', b0, b1); if (Buffer.compare(Buffer.from(x), Buffer.from(y))) same = false; }
    const z2 = await fz.fetch('demo2', 100, 1900), y2 = await fa.fetch('demo2', 100, 1900);
    ok(same && !Buffer.compare(Buffer.from(z2), Buffer.from(y2)), 'bgzip FASTA + .gzi matches plain FASTA');
    let threw = false; try { await new C.IndexedFasta(new C.FileSource(new Blob([Buffer.concat(blocks)]), 'demo.fa.gz'), new C.FileSource(blob('demo.fa.fai'))).open(); } catch (e) { threw = e.code === 'needGzi'; }
    ok(threw, '.fa.gz without .gzi gives a clear error');
  }
  console.log(fails ? `\n${fails} FAILED` : '\nAll checks passed');
  process.exit(fails ? 1 : 0);
})().catch(e => { console.error(e); process.exit(1); });
