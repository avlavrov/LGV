// Checks tools/bigbed.js against an independent reader (libBigWig, built from its source:
// see the usage line). Writes the bundled MANE transcripts as a bigBed and compares every record the
// reader returns with what was written, field for field.
// Usage: node tools/test-bigbed.js PATH/TO/bbdump   (bbdump: libBigWig's bigBed record dumper)
'use strict';
const fs = require('fs'), path = require('path'), os = require('os'), { execFileSync } = require('child_process');
const { writeBigBed } = require('./bigbed.js');
const dump = process.argv[2]; if (!dump) { console.error('usage: node tools/test-bigbed.js bbdump'); process.exit(2); }
let fails = 0; const ok = (c, m) => { console.log((c ? 'OK   ' : 'FAIL ') + m); if (!c) fails++; };
const rows = fs.readFileSync(path.join(__dirname, '..', 'data', 'mane.grch38.tsv'), 'utf8').split('\n').filter(l => l && l[0] !== '#').map(l => l.split('\t'));
const sizes = {}; for (const f of rows) sizes[f[1]] = Math.max(sizes[f[1]] || 0, +f[3] + +f[4] + 1000);
const recs = rows.map(f => { const s = +f[3]; return { chrom: f[1], start: s, end: s + +f[4], rest: [f[9], 0, f[2], s + +f[5], s + +f[6], 0, f[7].split(',').length, f[8] + ',', f[7] + ',', f[0], f[10]].join('\t') }; });
for (const ips of [64, 512, 4096]) {
  const buf = writeBigBed(recs, Object.entries(sizes), { autoSql: 'table t "test" ( string chrom; uint chromStart; uint chromEnd; string name; )', fieldCount: 14, definedFieldCount: 12, itemsPerSlot: ips, extraIndex: [{ field: 3, key: r => r.rest.split('\t')[0] }, { field: 12, key: r => r.rest.split('\t')[9] }] });
  const f = path.join(os.tmpdir(), 'lgv-test.bb'); fs.writeFileSync(f, buf);
  const out = execFileSync(dump, [f], { maxBuffer: 1 << 30 }).toString().trim().split('\n');
  const want = recs.map(r => [r.chrom, r.start, r.end, r.rest].join('\t'));
  const sort = a => a.slice().sort();
  const A = sort(out), B = sort(want); let diff = A.length !== B.length;
  for (let i = 0; i < A.length && !diff; i++) if (A[i] !== B[i]) { diff = true; console.log(' first difference:\n  ' + A[i] + '\n  ' + B[i]); }
  ok(!diff, `bigBed of ${recs.length} transcripts, ${ips} per block: libBigWig reads back every record unchanged (${(buf.length / 1048576).toFixed(2)} MB)`);
}
// LGV's own reader (src/signal.js) on the same file: regions and the name indexes
(async () => {
  const { BigBedFile } = require('../src/signal.js');
  const buf = writeBigBed(recs, Object.entries(sizes), { autoSql: 'table t "test (one per gene)" ( string chrom; uint chromStart; uint chromEnd; string name; uint score; char[1] strand; uint thickStart; uint thickEnd; uint reserved; int blockCount; int[blockCount] blockSizes; int[blockCount] chromStarts; string geneName; string refseq; )', fieldCount: 14, definedFieldCount: 12, itemsPerSlot: 512, extraIndex: [{ field: 3, key: r => r.rest.split('\t')[0].replace(/\.\d+$/, '').toUpperCase() }, { field: 12, key: r => r.rest.split('\t')[9].toUpperCase() }] });
  const src = { name: 'test.bb', size: buf.length, read: async (p, n) => new Uint8Array(buf.buffer, buf.byteOffset + p, Math.max(0, Math.min(n, buf.length - p))) };
  const bb = await new BigBedFile(src).open();
  let bad = 0, n = 0;
  for (const [chr, a, b] of [['chr17', 43000000, 43200000], ['chr7', 117400000, 117700000], ['chrX', 31000000, 33400000], ['chr1', 0, 3000000]]) {
    const got = await bb.features(chr, a, b), want = recs.filter(r => r.chrom === chr && r.start < b && r.end > a);
    n += got.length;
    if (got.length !== want.length) bad++;
    const w = new Map(want.map(r => [r.rest.split('\t')[0], r]));
    for (const f of got) { const r = w.get(f.name); if (!r || r.start !== f.start || r.end !== f.end || f.exons.length / 2 !== +r.rest.split('\t')[6] || f.attrs.geneName !== r.rest.split('\t')[9]) bad++; }
  }
  ok(!bad && n > 50, `LGV reads regions of the bigBed: ${n} transcripts in 4 regions, all fields as written`);
  const g = await bb.findName('brca1'), t = await bb.findName('ENST00000357654'), none = await bb.findName('NOSUCHGENE');
  ok(g.length === 1 && g[0].attrs.geneName === 'BRCA1' && t.length === 1 && t[0].start === g[0].start && none.length === 0, 'name indexes find a gene (any case) and a transcript id without version');
  console.log(fails ? `\n${fails} FAILED` : '\nAll checks passed'); process.exit(fails ? 1 : 0);
})().catch(e => { console.error(e); process.exit(1); });
