// Checks freshly built gene bigBeds before they are published (used by the monthly data update):
//   node tools/check-genes.js DIR [data/versions.json]
// - LGV's reader finds control genes by name on the chromosome and near the place MANE has them
//   (BRCA1, CFTR, DMD, SMN1, TTN) and IGHM on chr14 (not in MANE);
// - counts are plausible and not more than 5% below the last published ones.
// Prints a short report (used in the update e-mail); exits 1 on any failure.
'use strict';
const fs = require('fs'), path = require('path');
const { BigBedFile } = require('../src/signal.js');
const DIR = process.argv[2], VERS = process.argv[3] || path.join(__dirname, '..', 'data', 'versions.json');
const prev = fs.existsSync(VERS) ? JSON.parse(fs.readFileSync(VERS, 'utf8')) : {};
let fails = 0; const say = (c, m) => { console.log((c ? 'OK   ' : 'FAIL ') + m); if (!c) fails++; };
const open = async f => { const buf = fs.readFileSync(path.join(DIR, f)); return { bb: await new BigBedFile({ name: f, size: buf.length, read: async (p, n) => new Uint8Array(buf.buffer, buf.byteOffset + p, Math.max(0, Math.min(n, buf.length - p))) }).open(), size: buf.length }; };
const mane = new Map(fs.readFileSync(path.join(__dirname, '..', 'data', 'mane.grch38.tsv'), 'utf8').split('\n').filter(l => l && l[0] !== '#').map(l => l.split('\t')).filter(f => f[11] === 'S').map(f => [f[0], { chr: f[1], s: +f[3], e: +f[3] + +f[4] }]));
(async () => {
  const counts = {};
  for (const [k, f, min] of [['canonical', 'ensembl.canonical.grch38.bb', 50000], ['all', 'ensembl.all.grch38.bb', 180000]]) {
    const { bb, size } = await open(f);
    let n = 0; for (const [name, c] of bb.chrs) { void name; const r = await bb.features(c.name, 0, c.size, true); n += r.length; }
    counts[k] = n;
    const was = prev.genes && prev.genes[k] || 0;
    say(n >= min && (!was || n >= was * 0.95), `${k}: ${n} transcripts (last published ${was || 'none'}), ${(size / 1048576).toFixed(1)} MB`);
    for (const g of ['BRCA1', 'CFTR', 'DMD', 'SMN1', 'TTN']) {
      const hits = await bb.findName(g), m = mane.get(g);
      // all of them on the MANE chromosome, within 1 Mb of the MANE transcript (other transcripts of a
      // gene may start or end outside it), and at least one overlapping it
      say(hits.length > 0 && hits.every(h => h.chr === m.chr && h.end > m.s - 1e6 && h.start < m.e + 1e6) && hits.some(h => h.start < m.e && h.end > m.s), `${k}: ${g} found on ${m.chr} at its MANE place (${hits.length} transcript${hits.length === 1 ? '' : 's'})`);
    }
    const ig = await bb.findName('IGHM');
    say(ig.length > 0 && ig.every(h => h.chr === 'chr14' && h.start > 105500000 && h.end < 106900000), `${k}: IGHM found on chr14 (not in MANE)`);
    // the reader names the extra columns (gene symbol shown as the label)
    say(ig.length > 0 && ig.every(h => h.label === 'IGHM' && h.attrs && h.attrs.transcriptType), `${k}: column names read (label IGHM)`);
  }
  console.log(JSON.stringify({ counts }));
  process.exit(fails ? 1 : 0);
})().catch(e => { console.log('FAIL ' + (e && e.message || e)); process.exit(1); });
