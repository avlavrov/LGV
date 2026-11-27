// CGV: gene annotation — MANE Select and MANE Plus Clinical transcripts on GRCh38, bundled with CGV
// (data/mane.grch38.tsv, built by tools/make-mane.py). Gene search and the gene track work offline:
// no gene name typed by the user ever leaves the computer.
//
// Line format (0-based, half-open): symbol chr strand start length cdsStart cdsEnd exonStarts
// exonLengths ensembl refseq status name hgnc — offsets relative to start, lists comma-separated.
'use strict';
if (typeof contigKey === 'undefined' && typeof require !== 'undefined') Object.assign(globalThis, require('./core.js'));

const GENES = { ready: false, version: '', byChr: new Map(), byName: new Map(), count: 0 };

// The bundled table arrives packed (tools/pack-genes.js: deflate, base64) inside a text/plain
// script element.
async function loadGenes(b64) {
  const bin = atob(b64), u8 = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) u8[i] = bin.charCodeAt(i);
  const ds = new DecompressionStream('deflate-raw'), w = ds.writable.getWriter();
  w.write(u8).catch(() => {}); w.close().catch(() => {});
  unpackGenes(new Uint8Array(await new Response(ds.readable).arrayBuffer()));
}

const geneChrKey = name => contigKey(name);
// Collects transcripts (the same objects from the table text or the packed form).
function geneSink() {
  const byChr = new Map(), byName = new Map();
  const add = (k, tx) => { k = k.toUpperCase(); let a = byName.get(k); if (!a) byName.set(k, a = []); a.push(tx); };
  let n = 0;
  return {
    tx(sym, chr, strand, start, len, cds0, cds1, exons, enst, nm, clinical, name, hgnc) {
      const tx = { sym, chr, strand, start, end: start + len, cdsStart: start + cds0, cdsEnd: start + cds1, exons, enst, nm, clinical, name, hgnc };
      const k = geneChrKey(chr);
      let c = byChr.get(k); if (!c) byChr.set(k, c = { list: [], maxLen: 0 });
      c.list.push(tx); if (len > c.maxLen) c.maxLen = len;
      add(sym, tx); add(enst.replace(/\.\d+$/, ''), tx); add(enst, tx); add(nm.replace(/\.\d+$/, ''), tx); add(nm, tx);
      n++;
    },
    done() { for (const c of byChr.values()) c.list.sort((a, b) => a.start - b.start); Object.assign(GENES, { ready: true, byChr, byName, count: n }); }
  };
}
// The table as text (data/mane.grch38.tsv; tests and older builds).
function parseGenes(text) {
  const out = geneSink(), ints = s => s ? s.split(',').map(Number) : [];
  for (const line of text.split('\n')) {
    if (!line || line[0] === '#') { if (line.startsWith('#CGV MANE') || line.startsWith('#LGV MANE')) GENES.version = (line.match(/v[\d.]+/) || [''])[0]; continue; }
    const f = line.split('\t'); if (f.length < 14) continue;
    const start = +f[3], len = +f[4], es = ints(f[7]), el = ints(f[8]);
    if (!Number.isSafeInteger(start) || !(len > 0) || es.length !== el.length || !es.length) continue;
    const exons = new Int32Array(es.length * 2);
    for (let i = 0; i < es.length; i++) { exons[2 * i] = start + es[i]; exons[2 * i + 1] = start + es[i] + el[i]; }
    out.tx(f[0], f[1], f[2] === '-' ? -1 : 1, start, len, +f[5], +f[6], exons, f[9], f[10], f[11] === 'C', f[12], f[13]);
  }
  out.done();
}
// The packed form (tools/pack-genes.js), already inflated: a header, one text block, then
// variable-length integer streams read side by side. Bounds are checked; a damaged table leaves
// the gene track off rather than wrong.
function unpackGenes(u8) {
  const bad = () => { throw new Error('gene table damaged'); };
  if (u8.length < 8 || String.fromCharCode(...u8.subarray(0, 6)) !== 'LGVG1\n') bad();
  let p = 6;
  const vh = () => { let n = 0, m = 1; for (;;) { if (p >= u8.length) bad(); const b = u8[p++]; n += (b & 127) * m; if (b < 128) return n; m *= 128; if (m > 2 ** 49) bad(); } };
  const count = vh(), tlen = vh(), names = ['start', 'len', 'cdsS', 'cdsL', 'nex', 'exl', 'intr', 'enst', 'nm', 'hgnc'], lens = names.map(() => vh());
  if (p + tlen > u8.length) bad();
  const text = new TextDecoder().decode(u8.subarray(p, p + tlen)).split('\n'); p += tlen;
  const S = {}; for (let i = 0; i < names.length; i++) { if (p + lens[i] > u8.length) bad(); S[names[i]] = { a: u8, p, end: p + lens[i] }; p += lens[i]; }
  const v = st => { let n = 0, m = 1; for (;;) { if (st.p >= st.end) bad(); const b = st.a[st.p++]; n += (b & 127) * m; if (b < 128) return n; m *= 128; if (m > 2 ** 49) bad(); } };
  if (text.length !== 4 + 2 * count) bad();
  GENES.version = text[0];
  const chrs = text[1].split(',').map(x => { const i = x.lastIndexOf('='); return [x.slice(0, i), +x.slice(i + 1)]; });
  const strands = text[2], status = text[3];
  if (strands.length !== count || status.length !== count || chrs.reduce((a, c) => a + c[1], 0) !== count) bad();
  const out = geneSink(); let r = 0;
  for (const [chr, n] of chrs) {
    let pos = 0;
    for (let k = 0; k < n; k++, r++) {
      const start = (pos += v(S.start)), len = v(S.len), c0 = v(S.cdsS), c1 = c0 + v(S.cdsL), nx = v(S.nex);
      if (!(len > 0) || !nx || nx > 10000) bad();
      const exons = new Int32Array(nx * 2); let e = start;
      for (let i = 0; i < nx; i++) { if (i) e += v(S.intr); const l = v(S.exl); exons[2 * i] = e; exons[2 * i + 1] = e + l; e += l; }
      const enst = 'ENST' + String(v(S.enst)).padStart(11, '0') + '.' + v(S.enst);
      const nc = v(S.nm), nm = (nc % 2 ? 'NR_' : 'NM_') + String(Math.floor(nc / 4)).padStart(nc & 2 ? 9 : 6, '0') + '.' + v(S.nm);
      const h = v(S.hgnc);
      out.tx(text[4 + r], chr, strands[r] === '-' ? -1 : 1, start, len, c0, c1, exons, enst, nm, status[r] === 'C', text[4 + count + r], h ? 'HGNC:' + h : '');
    }
  }
  out.done();
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

if (typeof module !== 'undefined') module.exports = { CYTO, parseCyto, cytoIn, cytoBandAt, GENES, parseGenes, unpackGenes, genesIn, findGene, geneSuggest, txPart };
