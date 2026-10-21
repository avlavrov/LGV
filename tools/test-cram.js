// CRAM checks against the same data as BAM (test/cram/*.cram were made from demo/demo.bam with samtools,
// see tools/make-cram-tests.sh). Each stage of CRAM support adds its checks here.
'use strict';
const fs = require('fs'), path = require('path');
const C = require('../src/core.js'), R = require('../src/cram.js');
const T = path.join(__dirname, '..', 'test', 'cram'), D = path.join(__dirname, '..', 'demo');
const src = (dir, f) => new C.FileSource(new Blob([fs.readFileSync(path.join(dir, f))]), f);
let fails = 0;
const ok = (cond, msg) => { console.log((cond ? 'OK   ' : 'FAIL ') + msg); if (!cond) fails++; };
const throwsCode = async (f, code) => { try { await f(); return 'no error'; } catch (e) { return e.code === code ? true : 'got ' + (e.code || e.message); } };

(async () => {
  const bam = await new C.BamFile(src(D, 'demo.bam'), src(D, 'demo.bam.bai')).open();
  const files = fs.readdirSync(T).filter(f => f.endsWith('.cram'));
  ok(files.length >= 4, 'test CRAM files present: ' + files.join(', '));

  // ---- low-level pieces ----
  {
    const enc = v => { // ITF8 encoder for the test only
      v >>>= 0;
      if (v < 0x80) return [v];
      if (v < 0x4000) return [0x80 | (v >> 8), v & 0xFF];
      if (v < 0x200000) return [0xC0 | (v >> 16), (v >> 8) & 0xFF, v & 0xFF];
      if (v < 0x10000000) return [0xE0 | (v >>> 24), (v >> 16) & 0xFF, (v >> 8) & 0xFF, v & 0xFF];
      return [0xF0 | (v >>> 28), (v >>> 20) & 0xFF, (v >>> 12) & 0xFF, (v >>> 4) & 0xFF, v & 0x0F];
    };
    const vals = [0, 1, 127, 128, 16383, 16384, 2097151, 2097152, 268435455, 268435456, 2147483647, -1, -2, -4542278];
    const bad = vals.filter(v => { const r = new R.CramReader(Uint8Array.from(enc(v))); return r.itf8() !== v || r.pos !== r.end; });
    ok(!bad.length, 'ITF8 round trip on edge values' + (bad.length ? ' (bad: ' + bad + ')' : ''));
    ok(R.crc32(new TextEncoder().encode('123456789'), 0, 9) === 0xCBF43926, 'CRC32 check value');
    const cut = await throwsCode(() => new R.CramReader(Uint8Array.from([0xE0, 1])).itf8(), 'cramBad');
    ok(cut === true, 'a value cut short is a CgvError, not a crash (' + cut + ')');
  }

  const fa = await new C.IndexedFasta(src(D, 'demo.fa'), src(D, 'demo.fa.fai')).open();
  let unmappedInBam = 0; // demo data has none; kept for clarity of the count below
  // fields that must be equal; MD and NM are left out: samtools drops them when writing CRAM and
  // htslib regenerates them on reading only when asked
  const sameRead = (x, y) => {
    for (const k of ['name', 'flag', 'refID', 'pos', 'end', 'mapq', 'lseq', 'mateRef', 'matePos', 'tlen']) if (x[k] !== y[k]) return `${k} ${x[k]} vs ${y[k]}`;
    if (C.cigarString(x) !== C.cigarString(y)) return `CIGAR ${C.cigarString(x)} vs ${C.cigarString(y)}`;
    if (Buffer.compare(Buffer.from(x.seq), Buffer.from(y.seq))) return 'bases differ';
    if (Buffer.compare(Buffer.from(x.qual), Buffer.from(y.qual))) return 'qualities differ';
    const tx = Object.assign({}, C.recTags(x)), ty = Object.assign({}, C.recTags(y)); delete tx.MD; delete tx.NM; delete ty.MD; delete ty.NM;
    if (JSON.stringify(Object.entries(tx).sort()) !== JSON.stringify(Object.entries(ty).sort())) return `tags ${JSON.stringify(tx)} vs ${JSON.stringify(ty)}`;
    return null;
  };
  for (const f of files) {
    const cram = await new R.CramFile(src(T, f), src(T, f + '.crai')).open();
    const sameRefs = JSON.stringify(cram.refs) === JSON.stringify(bam.refs);
    ok(sameRefs, `${f}: CRAM ${cram.version}, contigs as in BAM (${cram.refs.map(r => r.name + ':' + r.len).join(', ')})`);
    ok(JSON.stringify(cram.samples) === JSON.stringify(bam.samples), `${f}: sample ${cram.samples.join(', ')} as in BAM`);
    ok(cram.refId('DEMO') === 0 && cram.refId('demo2') === 1, `${f}: contig names case-insensitive`);
    const list = await cram.scan();
    const recs = list.reduce((a, c) => a + c.header.nRecords, 0);
    const methods = [...new Set(list.flatMap(c => c.methods))];
    cram.decodable = methods.every(m => m === 'raw' || m === 'gzip');
    ok(list[list.length - 1].header.eof && recs === 8286, `${f}: ${list.length} containers to the end marker, all checksums correct, ${recs} records (BAM: 8286); block compression: ${methods.join(', ')}`);
    ok(cram.index.slices > 0 && cram.estimateBytes(0, 11000, 13000) > 0 && cram.estimateBytes(0, 0, 1) >= 0, `${f}: index has ${cram.index.slices} slices; size estimate for a window works`);
    // ---- records: every read equal to the BAM's (stage 2) ----
    if (!cram.decodable) {
      const why = await throwsCode(() => cram.fetch(0, 11000, 13000), 'cramCodec');
      ok(why === true, `${f}: blocks in codecs of later stages give a clear message`);
      continue;
    }
    cram.refProvider = async (id, a, b) => fa.fetch(cram.refs[id].name, a, b);
    let diffs = [], count = 0;
    for (let id = 0; id < bam.refs.length; id++) {
      const A = (await bam.fetch(id, 0, bam.refs[id].len)).reads, B = (await cram.fetch(id, 0, cram.refs[id].len)).reads;
      if (A.length !== B.length) diffs.push(`${bam.refs[id].name}: ${A.length} reads in BAM, ${B.length} in CRAM`);
      for (let i = 0; i < Math.min(A.length, B.length); i++) { const d = sameRead(A[i], B[i]); if (d) diffs.push(`${A[i].name}: ${d}`); }
      count += B.length;
    }
    ok(!diffs.length && count === 8286 - unmappedInBam, `${f}: all ${count} mapped reads equal to BAM (name, flags, position, CIGAR, bases, qualities, mate, TLEN, tags)` + (diffs.length ? ' — ' + diffs.length + ' differ: ' + diffs.slice(0, 4).join('; ') : ''));
    // the same window with downsampling and filters picks the same reads as BAM
    const opts = { filters: { minMapq: 10, hideDuplicates: true, hideSecondary: true, hideSupplementary: false, hideQcFail: true }, downsample: { max: 5, window: 50 } };
    const a = await bam.fetch(0, 11000, 13000, opts), b = await cram.fetch(0, 11000, 13000, opts);
    const keys = x => x.reads.map(r => r.name + ':' + r.flag + ':' + r.pos).join();
    let covSame = true; const ca = new Uint32Array(6), cb = new Uint32Array(6);
    for (let p = 11000; p < 13000 && covSame; p++) { a.cov.counts(p, ca); b.cov.counts(p, cb); covSame = ca.every((v, k) => v === cb[k]); }
    ok(keys(a) === keys(b) && a.total === b.total && covSame, `${f}: a filtered, downsampled window gives the same ${b.reads.length} reads, total ${b.total} and coverage as BAM`);
  }

  // ---- htslib's own test cases (test/cram/htslib: SAM files from htslib's test suite, made into BAM and
  // CRAM with samtools): clips, padding, unmapped mates, supplementary reads, odd tags, long reads ----
  {
    const H = path.join(T, 'htslib');
    const cases = fs.readdirSync(H).filter(f => f.endsWith('.cram')).map(f => f.slice(0, -5));
    let pass = 0; const bad = [];
    for (const n of cases) {
      try {
        const fa2 = await new C.IndexedFasta(src(H, n.split('_')[0] + '.fa'), src(H, n.split('_')[0] + '.fa.fai')).open();
        const b = await new C.BamFile(src(H, n + '.bam'), src(H, n + '.bam.bai')).open();
        const c = await new R.CramFile(src(H, n + '.cram'), src(H, n + '.cram.crai')).open();
        c.refProvider = async (id, a, e) => fa2.fetch(c.refs[id].name, a, e);
        const diffs = [];
        for (let id = 0; id < b.refs.length; id++) {
          const A = (await b.fetch(id, 0, 1 << 29)).reads, B = (await c.fetch(id, 0, 1 << 29)).reads; // whole contig, even reads past a too-short @SQ length
          if (A.length !== B.length) diffs.push(`${b.refs[id].name}: ${A.length} vs ${B.length} reads`);
          for (let i = 0; i < Math.min(A.length, B.length); i++) { const d = sameRead(A[i], B[i]); if (d) diffs.push(`${A[i].name}: ${d}`); }
        }
        if (diffs.length) bad.push(n + ': ' + diffs.slice(0, 2).join('; ')); else pass++;
      } catch (e) { bad.push(n + ': ' + (e.code ? e.message : e.stack.split('\n').slice(0, 2).join(' '))); }
    }
    ok(!bad.length, `htslib test cases: ${pass} of ${cases.length} CRAM files read exactly like their BAM` + (bad.length ? '\n     ' + bad.join('\n     ') : ''));
  }

  // ---- damaged files: a CgvError with a clear message, never a crash or a hang ----
  {
    const buf = fs.readFileSync(path.join(T, 'v31.cram'));
    const idx = src(T, 'v31.cram.crai');
    const open = u8 => new R.CramFile(new C.FileSource(new Blob([u8]), 'x.cram'), idx).open();
    const r1 = await throwsCode(() => open(Buffer.from('BAM\x01' + 'x'.repeat(40))), 'notCram');
    ok(r1 === true, 'a non-CRAM file is refused (' + r1 + ')');
    const v2 = Buffer.from(buf); v2[4] = 2; v2[5] = 1;
    const r2 = await throwsCode(() => open(v2), 'cramVersion');
    ok(r2 === true, 'CRAM 2.1 is refused with a clear message (' + r2 + ')');
    const flip = Buffer.from(buf); flip[40] ^= 0x55; // inside the header container's header or first block
    const r3 = await throwsCode(() => open(flip), 'cramCrc');
    ok(r3 === true, 'a changed byte is caught by the checksum (' + r3 + ')');
    const cut = buf.subarray(0, 60);
    const r4 = await throwsCode(() => open(cut), 'cramBad');
    ok(r4 === true, 'a cut-off file is refused (' + r4 + ')');
    const badIdx = await throwsCode(() => new R.CramFile(src(T, 'v31.cram'), new C.FileSource(new Blob(['not an index\n']), 'x.crai')).open(), 'notCrai');
    ok(badIdx === true, 'a broken .crai is refused (' + badIdx + ')');
    // many random corruptions: every one must end in a CgvError or open cleanly, within the time limit
    let crashes = 0, n = 0; const t0 = Date.now(); let seed = 12345;
    const rnd = () => (seed = (Math.imul(seed, 1103515245) + 12345) >>> 0) / 4294967296;
    for (; n < 400; n++) {
      const m = Buffer.from(buf.subarray(0, 4096 + Math.floor(rnd() * 8192)));
      for (let k = 0; k < 1 + Math.floor(rnd() * 4); k++) m[26 + Math.floor(rnd() * (m.length - 26))] = Math.floor(rnd() * 256);
      try { const c = await open(m); await c.scan(3).catch(e => { if (!e.code) throw e; }); }
      catch (e) { if (!e.code) { crashes++; if (crashes < 3) console.log('   crash:', e.message); } }
    }
    ok(crashes === 0, `fuzzing: ${n} corrupted files, ${crashes} ended in an unexpected error, ${Date.now() - t0} ms`);
    // the same for record decoding: damaged data containers of a readable CRAM, read through fetch
    {
      const g = fs.readFileSync(path.join(T, 'g30_small.cram')), gi = src(T, 'g30_small.cram.crai');
      const first = await new R.CramFile(new C.FileSource(new Blob([g]), 'g.cram'), gi).open();
      let crash2 = 0, okReads = 0; const t1 = Date.now(), codes = {};
      R.CRAM_CHECK.crc = false; // let damaged bytes through to the decoder itself
      for (let k = 0; k < 300; k++) {
        const m = Buffer.from(g);
        for (let j = 0; j < 1 + Math.floor(rnd() * 3); j++) m[first.firstContainer + Math.floor(rnd() * (m.length - first.firstContainer))] = Math.floor(rnd() * 256);
        try {
          const c = await new R.CramFile(new C.FileSource(new Blob([m]), 'g.cram'), gi).open();
          c.refProvider = async (id, a, b) => fa.fetch(c.refs[id].name, a, b);
          okReads += (await c.fetch(0, Math.floor(rnd() * 30000), 40000)).reads.length;
        } catch (e) { codes[e.code || 'unexpected'] = (codes[e.code || 'unexpected'] || 0) + 1; if (!e.code) { crash2++; if (crash2 < 3) console.log('   crash:', e.stack.split('\n').slice(0, 3).join(' | ')); } }
      }
      R.CRAM_CHECK.crc = true;
      ok(crash2 === 0, `fuzzing records (checksums off, so damage reaches the decoder): 300 damaged CRAM files read through fetch, ${crash2} ended in an unexpected error (clear errors: ${JSON.stringify(codes)}), ${Date.now() - t1} ms`);
    }
  }

  console.log(fails ? `\n${fails} check(s) failed` : '\nAll CRAM checks passed');
  process.exit(fails ? 1 : 0);
})().catch(e => { console.error(e); process.exit(1); });
