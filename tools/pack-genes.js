// Packs data/mane.grch38.tsv into the compact form bundled in LGV.html (read by unpackGenes in
// src/genes.js): text fields in one block, numbers as variable-length integers in separate streams
// (positions as differences, exon and intron lengths apart), the whole compressed with deflate.
// The same transcripts, field for field, as the table (checked by tools/test-core.js).
'use strict';
const zlib = require('zlib');
function packGenes(tsv) {
  const rows = tsv.split('\n').filter(l => l && l[0] !== '#').map(l => l.split('\t'));
  const ver = (tsv.match(/^#(?:CGV|LGV) MANE \S+ (v[\d.]+)/m) || [, ''])[1];
  const streams = { start: [], len: [], cdsS: [], cdsL: [], nex: [], exl: [], intr: [], enst: [], nm: [], hgnc: [] };
  const vi = (a, n) => { if (!Number.isSafeInteger(n) || n < 0) throw new Error('pack: bad number ' + n); while (n > 127) { a.push((n & 127) | 128); n = Math.floor(n / 128); } a.push(n); };
  const chrs = []; let pc = null, p = 0;
  for (const f of rows) {
    if (f.length < 14) throw new Error('pack: short line');
    if (f[1] !== pc) { chrs.push([f[1], 0]); pc = f[1]; p = 0; }
    chrs[chrs.length - 1][1]++;
    const start = +f[3]; vi(streams.start, start - p); p = start;
    vi(streams.len, +f[4]); vi(streams.cdsS, +f[5]); vi(streams.cdsL, +f[6] - +f[5]);
    const es = f[7].split(',').map(Number), el = f[8].split(',').map(Number);
    if (es.length !== el.length || es[0] !== 0) throw new Error('pack: exons');
    vi(streams.nex, es.length); let pos = 0;
    es.forEach((s, i) => { if (i) vi(streams.intr, s - pos); vi(streams.exl, el[i]); pos = s + el[i]; });
    const e = f[9].match(/^ENST(\d{11})\.(\d+)$/), r = f[10].match(/^(N[MR])_(\d+)\.(\d+)$/), h = f[13].match(/^HGNC:(\d+)$/);
    if (!e || !r || (f[13] && !h)) throw new Error('pack: ids ' + f[9] + ' ' + f[10] + ' ' + f[13]);
    vi(streams.enst, +e[1]); vi(streams.enst, +e[2]);
    if (r[2].length !== 6 && r[2].length !== 9) throw new Error('pack: RefSeq digits ' + f[10]);
    vi(streams.nm, +r[2] * 4 + (r[2].length === 9 ? 2 : 0) + (r[1] === 'NR' ? 1 : 0)); vi(streams.nm, +r[3]);
    vi(streams.hgnc, h ? +h[1] : 0);
    if (/[\n\t]/.test(f[0] + f[12])) throw new Error('pack: text');
  }
  const text = [ver, chrs.map(c => c.join('=')).join(','), rows.map(f => f[2]).join(''), rows.map(f => f[11]).join(''),
    ...rows.map(f => f[0]), ...rows.map(f => f[12])].join('\n');
  const parts = [], head = [];
  const tb = Buffer.from(text, 'utf8');
  vi(head, rows.length); vi(head, tb.length);
  for (const k of Object.keys(streams)) vi(head, streams[k].length);
  parts.push(Buffer.from('LGVG1\n'), Buffer.from(head), tb, ...Object.values(streams).map(a => Buffer.from(a)));
  return zlib.deflateRawSync(Buffer.concat(parts), { level: 9, memLevel: 9, windowBits: 15 });
}
module.exports = { packGenes };
