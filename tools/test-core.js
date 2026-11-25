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
    let same = true; for (let j = 0; j < 2000; j++) for (let k = 0; k < 6; k++) if (full.cov.get(11000 + j, k) !== ref[j * 6 + k]) same = false;
    ok(same, 'coverage from the stream equals coverage from records');
    const ds = await bam.fetch(0, 11000, 13000, { filters: flt, downsample: { max: 5, window: 50 } });
    const perWin = {}; for (const r of ds.reads) { const w = Math.floor(r.pos / 50); perWin[w] = (perWin[w] || 0) + 1; }
    ok(Math.max(...Object.values(perWin)) <= 10 && ds.total === full.reads.filter(r => r.pos >= 11000).length && ds.sampled.length > 0, `downsampling: ${ds.reads.length} of ${ds.total} kept, ${ds.sampled.length} windows sampled, max ${Math.max(...Object.values(perWin))} per window`);
    let same2 = true; for (let j = 11000; j < 13000; j++) for (let k = 0; k < 6; k++) if (ds.cov.get(j, k) !== full.cov.get(j, k)) same2 = false;
    ok(same2, 'coverage is the same with downsampling');
    ok(ds.reads.every((r, i) => !i || ds.reads[i - 1].pos <= r.pos), 'downsampled reads stay sorted');
    {
      // pairs: a kept read's mate (in the same query) is kept far more often than by a random choice
      const D = { max: 5, window: 50 }, wide = await bam.fetch(0, 10000, 16000, { filters: flt, downsample: D }), all = await bam.fetch(0, 10000, 16000, { filters: flt });
      const keptN = new Map(); for (const r of wide.reads) keptN.set(r.name, (keptN.get(r.name) || 0) + 1);
      const allN = new Map(); for (const r of all.reads) if (r.pos >= 10000) allN.set(r.name, (allN.get(r.name) || 0) + 1);
      let both = 0, one = 0; for (const [nm, c] of keptN) if (allN.get(nm) === 2) { if (c === 2) both++; else one++; }
      const frac = wide.reads.length / [...allN.values()].reduce((a, b) => a + b, 0);
      ok(both / (both + one) > 0.6, `downsampling keeps pairs: of the pairs with a read kept, ${both} kept both, ${one} one (${(100 * both / (both + one)).toFixed(0)}%; random would give about ${(100 * frac).toFixed(0)}%)`);
      // a name asked for is kept even where its window is full
      const lost = all.reads.find(r => r.pos >= 10000 && !keptN.has(r.name));
      const forced = await bam.fetch(0, 10000, 16000, { filters: flt, downsample: D, keep: [lost.name] });
      ok(forced.reads.some(r => r.name === lost.name) && forced.reads.every((r, i) => !i || forced.reads[i - 1].pos <= r.pos), 'a read asked for by name is kept through downsampling, order kept');
      const again = await bam.fetch(0, 10000, 16000, { filters: flt, downsample: D });
      ok(again.reads.map(r => r.name + r.pos).join() === wide.reads.map(r => r.name + r.pos).join(), 'the same reads are kept on every load');
    }
    {
      const fullC = await bam.fetch(0, 5000, 30000, { filters: flt }), dep = await bam.fetch(0, 5000, 30000, { filters: flt, depthOnly: true });
      const A = new Uint32Array(6), sm = dep.cov.summary(), want = new Uint32Array(sm.max5.length);
      for (let pos = 5000; pos < 30000; pos++) { const a = fullC.cov.counts(pos, A), t = a[0] + a[1] + a[2] + a[3] + a[4], k = Math.floor(pos / sm.bin) - sm.b0; if (t > want[k]) want[k] = t; }
      let bad = 0; for (let k = 0; k < want.length; k++) if (want[k] !== sm.max5[k]) bad++;
      ok(bad === 0 && sm.bin === 16 && dep.cov.depthOnly && dep.total === fullC.total, 'depth-only coverage: highest depth per 16 bp equals the sum of base counts (' + bad + ' bins differ)');
      const full64 = fullC.cov.summary(); let bad64 = 0;
      for (let k = 0; k < full64.max5.length; k++) { const p0 = (full64.b0 + k) * 64; let x = 0; for (let p = Math.max(p0, 5000); p < Math.min(p0 + 64, 30000); p++) { const a = fullC.cov.counts(p, A); x = Math.max(x, a[0] + a[1] + a[2] + a[3] + a[4]); } if (x !== full64.max5[k]) bad64++; }
      ok(bad64 === 0, 'the 64 bp summary of detailed coverage matches the base counts');
      // read in two parts (as the page does with several workers) and merged: the same summary
      const p1 = await bam.fetch(0, 5000, 17008, { filters: flt, depthOnly: true, covBin: 16 }), p2 = await bam.fetch(0, 17008, 30000, { filters: flt, depthOnly: true, covBin: 16 });
      const mg = new C.CovSummary(5000, 30000, 16);
      for (const [p, c0] of [[p1, 5000], [p2, 17008]]) { const s1 = p.cov.summary(), off = Math.floor(c0 / 16) - mg.b0; s1.max5.forEach((v, j) => { if (v > mg.max5[off + j]) mg.max5[off + j] = v; }); }
      ok(mg.max5.every((v, k) => v === sm.max5[k]), 'a window read in two parts gives the same depth summary');
      // a 200 kb aligned stretch crosses 64 kb blocks with no other marks; a second read starts after it
      const dd = new C.DepthDiff(0, 400000, 64), cv = new C.SparseCoverage(0, 400000);
      dd.addRec({ pos: 1000, cigar: [(200000 << 4) | 0] }); dd.addRec({ pos: 150000, cigar: [(100 << 4) | 0] }); dd.addRec({ pos: 300000, cigar: [(50 << 4) | 0] }); dd.finish(cv);
      const ds = cv.summary();
      ok(ds.maxIn(1000, 2000) === 1 && ds.maxIn(130000, 140000) === 1 && ds.maxIn(150000, 150100) === 2 && ds.maxIn(201024, 299000) === 0 && ds.maxIn(300000, 300050) === 1 && ds.maxIn(0, 960) === 0, 'depth-only: a long read keeps its depth through empty 64 kb blocks');
      ok(dep.reads.length === 0, 'a depth-only (wide) window returns no reads');
      const ix = bam.indexDepth(0, 0, 40000);
      ok(ix && ix.win === 16384 && ix.bytes.length === 3 && ix.bytes.every(v => Number.isFinite(v) && v >= 0) && ix.bytes.reduce((x, y) => x + y, 0) > 0, 'depth from the index: bytes per 16 kb window (' + (ix ? [...ix.bytes].map(Math.round).join(', ') : 'none') + ')');
      const capped = await bam.fetch(0, 5000, 30000, { filters: flt, maxReads: 1000 });
      let same3 = true; for (let j = 5000; j < 30000; j += 7) for (let k = 0; k < 6; k++) if (capped.cov.get(j, k) !== fullC.cov.get(j, k)) same3 = false;
      ok(capped.tooMany && capped.reads.length === 0 && same3 && capped.total === fullC.total, `past maxReads (${fullC.reads.length} > 1000): no reads, the limit is reported, coverage complete`);
    }
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
    const b38 = C.matchContigs([{ name: 'chr1', len: 248956422 }, { name: 'chrM', len: 16569 }], [{ name: 'chr1', len: 249250621 }, { name: 'chrM', len: 16571 }]);
    ok(b38.matched.length === 0 && b38.conflicts.length === 2 && b38.conflicts[0][0] === 'chr1', 'same names with other lengths (another build) are conflicts, not matches');
  }
  // stable layout
  {
    const all = (await bam.fetch(0, 0, 40000)).reads;
    const win = (a, b) => all.filter(r => r.end > a && r.pos < b);
    const noOverlap = groups => groups.every(g => g.rows.every(row => row.every((r, i) => !i || row[i - 1].end + 2 <= r.pos)));
    const L = new C.ReadLayout('t', 2);
    const A = L.update(win(10000, 20000), () => '');
    const depthAt = p => all.filter(r => r.pos <= p && r.end + 2 > p).length;
    const rowOf = res => { const m = new Map(); for (const g of res.groups) g.rows.forEach((row, i) => row.forEach(r => m.set(C.readKey(r), i))); return m; };
    const ra = rowOf(A);
    ok(noOverlap(A.groups), 'layout: no overlaps in rows');
    ok(win(10000, 20000).every(r => ra.get(C.readKey(r)) < depthAt(r.pos)), 'layout: each read sits no lower than the depth at its start (rows filled from the top)');
    const B = L.update(win(15000, 25000), () => ''), rb = rowOf(B);
    const common = win(15000, 20000);
    ok(common.every(r => ra.get(C.readKey(r)) === rb.get(C.readKey(r))) && noOverlap(B.groups), `layout: ${common.length} reads shared by two windows keep their rows`);
    const B2 = L.update(win(15000, 25000), () => ''), rb2 = rowOf(B2);
    ok([...rb].every(([k, v]) => rb2.get(k) === v), 'layout: same reads again give the same rows');
    const back = L.update(win(5000, 15000), () => ''), rk = rowOf(back);
    const c2 = all.filter(r => r.pos < 15000 && r.end > 15000); // loaded in both windows
    ok(c2.length > 0 && c2.every(r => rk.get(C.readKey(r)) === rb.get(C.readKey(r))) && noOverlap(back.groups), `layout: panning back, ${c2.length} reads loaded in both windows keep their rows`);
    // sorting places reads covering the position at the top, in key order
    const S = new C.ReadLayout('s', 2), pos = 12000;
    const prio = { pos, key: r => [r.flag & 16 ? 1 : 0, r.mapq] };
    const res = S.update(win(11000, 13000), () => '', prio), rows = res.groups[0].rows;
    const cover = win(11000, 13000).filter(r => r.pos <= pos && r.end > pos);
    const top = rows.slice(0, cover.length).map(row => row.find(r => r.pos <= pos && r.end > pos));
    ok(top.every(Boolean) && top.every((r, i) => !i || prio.key(top[i - 1])[0] <= prio.key(r)[0]) && noOverlap(res.groups), `sort: ${cover.length} reads at ${pos} fill the top rows in key order`);
    // downsampling gives the same reads for a window whatever the query start
    const flt = { minMapq: 0, hideDuplicates: true, hideSecondary: true, hideSupplementary: false, hideQcFail: true }, ds = { max: 5, window: 50 };
    const q1 = await bam.fetch(0, 10000, 14000, { filters: flt, downsample: ds }), q2 = await bam.fetch(0, 11000, 16000, { filters: flt, downsample: ds });
    // (past an insert size from the later start: a mate before the query start cannot decide for its pair)
    const keys = q => new Set(q.reads.filter(r => r.pos >= 12000 && r.pos < 14000).map(C.readKey));
    const k1 = keys(q1), k2 = keys(q2);
    ok(k1.size === k2.size && [...k1].every(k => k2.has(k)), `downsampling: ${k1.size} kept reads in 12000-14000 are the same for both queries`);
  }
  // sparse coverage over a huge window allocates blocks only where reads are
  {
    const huge = await bam.fetch(0, 0, 40000000 > 40000 ? 40000 : 0);
    const wide = new C.SparseCoverage(0, 200000000);
    ok(wide.blocks.size === 0 && wide.get(150000000, 0) === 0, 'empty sparse coverage costs nothing');
    ok(huge.cov.blocks.size <= 1, `coverage blocks for the demo contig: ${huge.cov.blocks.size}`);
    const r = new C.SparseRef(); r.add(100, new Uint8Array([65, 67])); r.add(10, new Uint8Array([71]));
    ok(r.at(10) === 71 && r.at(101) === 67 && r.at(50) === 0, 'sparse reference lookup');
    const rec = huge.reads.find(x => C.cigarString(x).includes('3D'));
    const dpos = rec ? rec.pos + parseInt(C.cigarString(rec)) : -1;
    ok(rec && C.queryAt(rec, dpos).del && C.queryAt(rec, rec.pos).qi === 0, 'base under a click: query index or deletion');
  }
  // pair view units
  {
    const all = (await bam.fetch(0, 0, 40000)).reads;
    const units = C.pairUnits(all, 10000);
    const both = units.filter(u => u.reads.length === 2), single = units.filter(u => u.reads.length === 1);
    ok(both.length > 3000 && both.every(u => u.reads[0].name === u.reads[1].name && u.pos === Math.min(u.reads[0].pos, u.reads[1].pos) && u.end === Math.max(u.reads[0].end, u.reads[1].end)), `pairs: ${both.length} pairs joined, extent spans both mates`);
    ok(single.some(u => u.reads[0].mateRef !== u.reads[0].refID), 'mates on another contig stay single');
    ok(units.every((u, i) => !i || units[i - 1].pos <= u.pos), 'pair units sorted by start');
    // a half pair (mate not loaded) gets the same extent as the full pair when MC is present, so layout is stable
    const win = all.filter(r => r.end > 12000 && r.pos < 13000), half = C.pairUnits(win, 10000).filter(u => u.reads.length === 1 && u.mate);
    const full = new Map(units.map(u => [u.key, u]));
    const hasMC = win.some(r => C.recTags(r).MC);
    ok(half.length > 0 && (!hasMC || half.every(u => full.get(u.key).pos === u.pos && full.get(u.key).end === u.end)) || (half.length > 0 && !hasMC), `half pairs at the window edge: ${half.length} (MC tags ${hasMC ? 'present' : 'absent'})`);
    const L = new C.ReadLayout('p', 2), A = L.update(C.pairUnits(win, 10000), () => '', null, u => u.key);
    const B = L.update(C.pairUnits(all.filter(r => r.end > 11500 && r.pos < 13500), 10000), () => '', null, u => u.key);
    const rows = res => { const m = new Map(); res.groups[0].rows.forEach((row, i) => row.forEach(u => m.set(u.key, i))); return m; };
    const ra = rows(A), rb = rows(B); let moved = 0, n = 0; for (const [k, v] of ra) if (rb.has(k)) { n++; if (rb.get(k) !== v) moved++; }
    const noOv = B.groups[0].rows.every(row => row.every((u, i) => !i || row[i - 1].end + 2 <= u.pos));
    ok(noOv, `pair layout: no overlaps after the window grows; ${moved} of ${n} units moved (only grown pairs may move)`);
    ok(C.cigarStrRefLen('10S50M2I3D40M5N10M') === 108, 'reference length of an MC cigar');
  }
  console.log(fails ? `\n${fails} FAILED` : '\nAll checks passed');
  process.exit(fails ? 1 : 0);
})().catch(e => { console.error(e); process.exit(1); });
