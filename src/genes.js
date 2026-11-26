// CGV: gene annotation — MANE Select and MANE Plus Clinical transcripts on GRCh38, bundled with CGV
// (data/mane.grch38.tsv, built by tools/make-mane.py). Gene search and the gene track work offline:
// no gene name typed by the user ever leaves the computer.
//
// Line format (0-based, half-open): symbol chr strand start length cdsStart cdsEnd exonStarts
// exonLengths ensembl refseq status name hgnc — offsets relative to start, lists comma-separated.
'use strict';
if (typeof contigKey === 'undefined' && typeof require !== 'undefined') Object.assign(globalThis, require('./core.js'));

const GENES = { ready: false, version: '', byChr: new Map(), byName: new Map(), count: 0 };

// The bundled table arrives gzip-compressed and base64-encoded inside a text/plain script element.
async function loadGenes(b64) {
  const bin = atob(b64), u8 = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) u8[i] = bin.charCodeAt(i);
  const ds = new DecompressionStream('gzip'), w = ds.writable.getWriter();
  w.write(u8).catch(() => {}); w.close().catch(() => {});
  parseGenes(await new Response(ds.readable).text());
}

const geneChrKey = name => contigKey(name);
function parseGenes(text) {
  const byChr = new Map(), byName = new Map();
  const add = (k, tx) => { k = k.toUpperCase(); let a = byName.get(k); if (!a) byName.set(k, a = []); a.push(tx); };
  const ints = s => s ? s.split(',').map(Number) : [];
  let n = 0;
  for (const line of text.split('\n')) {
    if (!line || line[0] === '#') { if (line.startsWith('#CGV MANE') || line.startsWith('#LGV MANE')) GENES.version = (line.match(/v[\d.]+/) || [''])[0]; continue; }
    const f = line.split('\t'); if (f.length < 14) continue;
    const start = +f[3], len = +f[4], es = ints(f[7]), el = ints(f[8]);
    if (!Number.isSafeInteger(start) || !(len > 0) || es.length !== el.length || !es.length) continue;
    const exons = new Int32Array(es.length * 2);
    for (let i = 0; i < es.length; i++) { exons[2 * i] = start + es[i]; exons[2 * i + 1] = start + es[i] + el[i]; }
    const tx = { sym: f[0], chr: f[1], strand: f[2] === '-' ? -1 : 1, start, end: start + len, cdsStart: start + +f[5], cdsEnd: start + +f[6],
      exons, enst: f[9], nm: f[10], clinical: f[11] === 'C', name: f[12], hgnc: f[13] };
    const k = geneChrKey(tx.chr);
    let c = byChr.get(k); if (!c) byChr.set(k, c = { list: [], maxLen: 0 });
    c.list.push(tx); if (len > c.maxLen) c.maxLen = len;
    add(tx.sym, tx); add(tx.enst.replace(/\.\d+$/, ''), tx); add(tx.enst, tx); add(tx.nm.replace(/\.\d+$/, ''), tx); add(tx.nm, tx);
    n++;
  }
  for (const c of byChr.values()) c.list.sort((a, b) => a.start - b.start);
  Object.assign(GENES, { ready: true, byChr, byName, count: n });
}

// Transcripts on a contig (any naming: chr7, 7) overlapping [a, b), sorted by start.
function genesIn(contig, a, b) {
  const c = GENES.byChr.get(geneChrKey(contig)); if (!c) return [];
  const L = c.list; let lo = 0, hi = L.length;
  const from = a - c.maxLen;
  while (lo < hi) { const m = (lo + hi) >> 1; if (L[m].start < from) lo = m + 1; else hi = m; }
  const out = [];
  for (let i = lo; i < L.length && L[i].start < b; i++) if (L[i].end > a) out.push(L[i]);
  return out;
}
// Transcripts for a gene symbol, an Ensembl or RefSeq transcript id (with or without version).
// MANE Select first, then MANE Plus Clinical.
function findGene(q) {
  const a = GENES.byName.get(String(q).trim().toUpperCase());
  return a ? a.slice().sort((x, y) => x.clinical - y.clinical) : [];
}
// Symbols starting with a prefix, for suggestions in the locus field.
function geneSuggest(prefix, n) {
  const p = String(prefix).trim().toUpperCase(), out = [];
  if (p.length < 2 || !GENES.ready) return out;
  for (const k of GENES.byName.keys()) if (k.startsWith(p) && !/^(ENST|NM_)/.test(k)) { out.push(GENES.byName.get(k)[0].sym); if (out.length >= n) break; }
  return [...new Set(out)];
}
// Where a position falls in a transcript: { exon, of } or { intron, of } (numbered in the
// transcript's direction), plus whether it is coding; null outside the transcript.
function txPart(tx, pos) {
  if (pos < tx.start || pos >= tx.end) return null;
  const n = tx.exons.length / 2, num = i => tx.strand > 0 ? i + 1 : n - i;
  for (let i = 0; i < n; i++) {
    const s = tx.exons[2 * i], e = tx.exons[2 * i + 1];
    if (pos >= s && pos < e) return { exon: num(i), of: n, coding: pos >= tx.cdsStart && pos < tx.cdsEnd };
    if (i + 1 < n && pos >= e && pos < tx.exons[2 * i + 2]) return { intron: tx.strand > 0 ? i + 1 : n - i - 1, of: n - 1 };
  }
  return null;
}

// ---- Chromosome bands (GRCh38, bundled: data/cytoband.grch38.tsv) for the chromosome map ----
const CYTO = new Map(); // contig key → [{ s, e, name, stain }]
function parseCyto(text) {
  for (const line of text.split('\n')) {
    if (!line || line[0] === '#') continue;
    const f = line.split('\t'); if (f.length < 5) continue;
    const s = +f[1], e = +f[2]; if (!Number.isSafeInteger(s) || !Number.isSafeInteger(e) || e <= s) continue;
    const k = geneChrKey(f[0]); let a = CYTO.get(k); if (!a) CYTO.set(k, a = []);
    a.push({ s, e, name: f[3], stain: f[4] });
  }
}
function cytoIn(contig) { return CYTO.get(geneChrKey(contig)) || null; }
function cytoBandAt(contig, pos) { const a = cytoIn(contig); if (!a) return null; const b = a.find(x => pos >= x.s && pos < x.e); return b ? b.name : null; }

if (typeof module !== 'undefined') module.exports = { CYTO, parseCyto, cytoIn, cytoBandAt, GENES, parseGenes, genesIn, findGene, geneSuggest, txPart };
