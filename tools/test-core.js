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
  const all = (await bam.fetch(0, 0, 40000)).reads;
  const all2 = (await bam.fetch(1, 0, 2000)).reads;
  ok(all.length + all2.length === truth.reads, `all reads: ${all.length}+${all2.length} = ${truth.reads} (${Date.now() - t} ms)`);
  ok(new Set(all.map(r => r.name + (r.flag & 64))).size === all.length, 'no duplicates after merging chunks');

  const at = (await bam.fetch(0, 12000, 12001)).reads;
  ok(at.length === truth.overlap_12000_12001, `overlapping 12000: ${at.length} = ${truth.overlap_12000_12001}`);
  const at2 = (await bam.fetch(1, 1000, 1200)).reads;
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
  // FASTA without .fai: the index built in the browser matches the .fai file
  {
    const zlib = require('zlib'), raw = fs.readFileSync(path.join(D, 'demo.fa'));
    const want = fs.readFileSync(path.join(D, 'demo.fa.fai'), 'utf8').trim();
    const asFai = idx => [...idx].map(([n, e]) => [n, e.len, e.off, e.lb, e.lw].join('\t')).join('\n');
    const plain = await new C.IndexedFasta(new C.FileSource(new Blob([raw]), 'demo.fasta')).openUnindexed(raw.length);
    ok(asFai(plain.idx) === want, 'index built from plain FASTA equals .fai');
    const p1 = await plain.fetch('demo', 12345, 12400), p2 = await fa.fetch('demo', 12345, 12400);
    ok(!Buffer.compare(Buffer.from(p1), Buffer.from(p2)), 'bases from the built index match');
    // CRLF line endings
    const crlf = Buffer.from(raw.toString('latin1').replace(/\n/g, '\r\n'), 'latin1');
    const fc = await new C.IndexedFasta(new C.FileSource(new Blob([crlf]), 'crlf.fa')).openUnindexed(crlf.length);
    const c1 = await fc.fetch('demo2', 50, 1950), c2 = await fa.fetch('demo2', 50, 1950);
    ok(!Buffer.compare(Buffer.from(c1), Buffer.from(c2)), 'CRLF FASTA indexed and read correctly');
    // ordinary gzip (not bgzip)
    const gz = zlib.gzipSync(raw);
    const fg = await new C.IndexedFasta(new C.FileSource(new Blob([gz]), 'demo.fa.gz')).openUnindexed(gz.length);
    const g1 = await fg.fetch('demo', 30000, 30100), g2 = await fa.fetch('demo', 30000, 30100);
    ok(!Buffer.compare(Buffer.from(g1), Buffer.from(g2)), 'plain gzip FASTA indexed and read');
    // bgzip without .gzi: the .gzi table is rebuilt
    const blocks = []; for (let u = 0; u < raw.length; u += 900) {
      const d = raw.subarray(u, u + 900), c = zlib.deflateRawSync(d);
      const h = Buffer.from([31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, 66, 67, 2, 0, 0, 0]); h.writeUInt16LE(c.length + 25, 16);
      const tl = Buffer.alloc(8); tl.writeUInt32LE(zlib.crc32(d), 0); tl.writeUInt32LE(d.length, 4); blocks.push(h, c, tl);
    }
    blocks.push(Buffer.from('1f8b08040000000000ff0600424302001b0003000000000000000000', 'hex'));
    const bz = Buffer.concat(blocks);
    const fb = await new C.IndexedFasta(new C.FileSource(new Blob([bz]), 'demo.fa.gz')).openUnindexed(bz.length);
    let same = asFai(fb.idx) === want && !!fb.gzi;
    for (const [s0, s1] of [[0, 5], [880, 2100], [39000, 40000]]) { const x = await fb.fetch('demo', s0, s1), y = await fa.fetch('demo', s0, s1); if (Buffer.compare(Buffer.from(x), Buffer.from(y))) same = false; }
    ok(same, 'bgzip FASTA without .gzi: index and .gzi rebuilt, bases match');
    // lines of different lengths in the middle of a contig: clear error
    let code = null; try { await new C.IndexedFasta(new C.FileSource(new Blob(['>x\nACGT\nAC\nACGT\n']), 'bad.fa')).openUnindexed(17); } catch (e) { code = e.code; }
    ok(code === 'faIrregular', 'irregular line lengths give a clear error');
  }
  // streaming fetch: coverage counts every read, downsampling caps reads per window
  {
    const flt = { minMapq: 0, hideDuplicates: true, hideSecondary: true, hideSupplementary: false, hideQcFail: true };
    const full = await bam.fetch(0, 11000, 13000, { filters: flt });
    const ref = C.computeCoverage(full.reads, 11000, 2000);
    ok(full.cov.every((v, i) => v === ref[i]), 'coverage from the stream equals coverage from records');
    const ds = await bam.fetch(0, 11000, 13000, { filters: flt, downsample: { max: 5, window: 50 } });
    const perWin = {}; for (const r of ds.reads) { const w = Math.floor(r.pos / 50); perWin[w] = (perWin[w] || 0) + 1; }
    ok(Math.max(...Object.values(perWin)) <= 5 && ds.total === full.total && ds.sampled.length > 0, `downsampling: ${ds.reads.length} of ${ds.total} kept, ${ds.sampled.length} windows sampled, max ${Math.max(...Object.values(perWin))} per window`);
    ok(ds.cov.every((v, i) => v === full.cov[i]), 'coverage is the same with downsampling');
    ok(ds.reads.every((r, i) => !i || ds.reads[i - 1].pos <= r.pos), 'downsampled reads stay sorted');
    const again = await bam.fetch(0, 11000, 13000, { filters: flt });
    ok(again.reads.length === full.reads.length && again.total === full.total, 'second fetch through the block cache gives the same reads');
    ok(['demo', 'Demo', 'DEMO'].every(n => bam.refId(n) === 0) && C.contigKey('Chr7') === C.contigKey('7') && C.contigKey('chrM') === C.contigKey('MT'), 'contig names are case-insensitive, chr prefix and M/MT optional');
  }
  // matching BAM contigs to differently named reference sequences
  {
    const bamRefs = [{ name: 'LRG_1730', len: 5000 }, { name: 'chr7', len: 159345973 }, { name: 'NG_012345.1', len: 777 }, { name: 'Contig9', len: 4242 }];
    const fa = [{ name: 'gi|224589800|ref|NC_000007.13|', len: 1 }, { name: 'LRG_1730t1', len: 5000 }, { name: 'Chr7', len: 159345973 }, { name: 'NG_012345', len: 777 }, { name: 'other', len: 999 }];
    const m = C.matchContigs(bamRefs, fa), got = Object.fromEntries(m.matched.map(x => [x[0], x[1] + '/' + x[2]]));
    ok(got.chr7 === 'Chr7/name' && got['NG_012345.1'] === 'NG_012345/name' && got.LRG_1730 === 'LRG_1730t1/length' && m.unmatched.join() === 'Contig9',
      'contig matching: case, version, unique length; ' + JSON.stringify(got) + ' unmatched ' + m.unmatched);
    const ncbi = C.matchContigs([{ name: 'NC_000007.13', len: 10 }], [{ name: 'gi|224589800|ref|NC_000007.13|', len: 10 }]);
    ok(ncbi.matched.length === 1 && ncbi.matched[0][2] === 'name', 'NCBI-style FASTA headers match by accession');
    const amb = C.matchContigs([{ name: 'a', len: 5 }, { name: 'b', len: 5 }], [{ name: 'x', len: 5 }, { name: 'y', len: 5 }]);
    ok(amb.matched.length === 0, 'equal lengths on several contigs are not guessed');
  }
  console.log(fails ? `\n${fails} FAILED` : '\nAll checks passed');
  process.exit(fails ? 1 : 0);
})().catch(e => { console.error(e); process.exit(1); });
