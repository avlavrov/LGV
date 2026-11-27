// Builds LGV's Ensembl gene tracks (bigBed, read by region; tools/bigbed.js) from a GTF:
//   node tools/make-genes.js GENES.gtf.gz OUTDIR [--mane data/mane.grch38.tsv]
// writes OUTDIR/ensembl.all.grch38.bb (every transcript) and OUTDIR/ensembl.canonical.grch38.bb (one
// per gene). One per gene: the transcript tagged Ensembl_canonical; for a GTF without that tag (e.g.
// UCSC's), the gene's MANE Select transcript, else its longest coding part, else its longest
// transcript (the rule is written into the file's description). Contigs: chr1-22, X, Y, M.
// Record: BED 12 (name = transcript id, thick = CDS incl. stop codon) + geneName, geneId,
// transcriptType, tags. Name indexes: transcript id, gene name, gene id (upper case, no version).
'use strict';
const fs = require('fs'), zlib = require('zlib'), path = require('path'), readline = require('readline');
const { writeBigBed } = require('./bigbed.js');
const [IN, OUT] = process.argv.slice(2);
const maneArg = process.argv.indexOf('--mane');
if (!IN || !OUT) { console.error('usage: node tools/make-genes.js GENES.gtf[.gz] OUTDIR [--mane data/mane.grch38.tsv]'); process.exit(2); }
const MANE = maneArg > 0 ? process.argv[maneArg + 1] : path.join(__dirname, '..', 'data', 'mane.grch38.tsv');
const SIZES = { chr1: 248956422, chr2: 242193529, chr3: 198295559, chr4: 190214555, chr5: 181538259, chr6: 170805979, chr7: 159345973, chr8: 145138636, chr9: 138394717, chr10: 133797422, chr11: 135086622, chr12: 133275309, chr13: 114364328, chr14: 107043718, chr15: 101991189, chr16: 90338345, chr17: 83257441, chr18: 80373285, chr19: 58617616, chr20: 64444167, chr21: 46709983, chr22: 50818468, chrX: 156040895, chrY: 57227415, chrM: 16569 };
const chrName = c => { const n = c.startsWith('chr') ? c : 'chr' + c; return n === 'chrMT' ? 'chrM' : n; };
const strip = id => String(id || '').replace(/\.\d+$/, '');

(async () => {
  const maneSelect = new Set(fs.readFileSync(MANE, 'utf8').split('\n').filter(l => l && l[0] !== '#').map(l => l.split('\t')).filter(f => f[11] === 'S').map(f => strip(f[9])));
  const tx = new Map(); let lines = 0, sawCanonTag = false;
  const input = fs.createReadStream(IN).pipe(IN.endsWith('.gz') ? zlib.createGunzip() : new (require('stream').PassThrough)());
  const rl = readline.createInterface({ input, crlfDelay: Infinity });
  for await (const line of rl) {
    if (!line || line[0] === '#') continue;
    lines++;
    const f = line.split('\t'); if (f.length < 9) continue;
    const type = f[2]; if (type !== 'exon' && type !== 'CDS' && type !== 'stop_codon' && type !== 'start_codon') continue;
    const chr = chrName(f[0]); if (!SIZES[chr]) continue;
    const tid = (f[8].match(/transcript_id "([^"]*)"/) || [])[1]; if (!tid) continue;
    let t = tx.get(tid);
    if (!t) {
      // attributes are read once per transcript (the first line of it); tags from every line
      const at = {}; for (const m of f[8].matchAll(/(\w+) "([^"]*)"/g)) if (m[1] !== 'tag' && !(m[1] in at)) at[m[1]] = m[2];
      tx.set(tid, t = { id: tid, chr, strand: f[6], ex: [], c0: Infinity, c1: -Infinity, gene: at.gene_name || at.gene || '', geneId: at.gene_id || '', type: at.transcript_type || at.transcript_biotype || at.gene_type || at.gene_biotype || '', tags: '' });
    }
    for (const m of f[8].matchAll(/tag "([^"]*)"/g)) if (/^(Ensembl_canonical|MANE_Select|MANE_Plus_Clinical|basic|CCDS)$/.test(m[1]) && !t.tags.split(',').includes(m[1])) { t.tags = t.tags ? t.tags + ',' + m[1] : m[1]; if (m[1] === 'Ensembl_canonical') sawCanonTag = true; }
    const s = +f[3] - 1, e = +f[4];
    if (type === 'exon') t.ex.push(s, e);
    else { if (s < t.c0) t.c0 = s; if (e > t.c1) t.c1 = e; }
  }
  const all = [];
  for (const t of tx.values()) {
    if (!t.ex.length) continue;
    const pairs = []; for (let i = 0; i < t.ex.length; i += 2) pairs.push([t.ex[i], t.ex[i + 1]]);
    t.exons = pairs.sort((a, b) => a[0] - b[0]); t.ex = null; t.cds = [t.c0, t.c1];
    t.tags = new Set(t.tags ? t.tags.split(',') : []);
    t.start = t.exons[0][0]; t.end = t.exons[t.exons.length - 1][1];
    if (strip(t.id) && maneSelect.has(strip(t.id))) t.tags.add('MANE_Select');
    t.cdsLen = 0; if (t.cds[1] > t.cds[0]) for (const [a, b] of t.exons) t.cdsLen += Math.max(0, Math.min(b, t.cds[1]) - Math.max(a, t.cds[0]));
    t.len = t.exons.reduce((n, [a, b]) => n + b - a, 0);
    all.push(t);
  }
  // one per gene
  const byGene = new Map();
  for (const t of all) { const k = t.geneId || t.gene || t.id; if (!byGene.has(k)) byGene.set(k, []); byGene.get(k).push(t); }
  const canon = [];
  for (const list of byGene.values()) {
    let pick = sawCanonTag ? list.find(t => t.tags.has('Ensembl_canonical')) : null;
    if (!pick) pick = list.find(t => t.tags.has('MANE_Select'));
    if (!pick) pick = list.slice().sort((a, b) => b.cdsLen - a.cdsLen || b.len - a.len || (a.id < b.id ? -1 : 1))[0];
    canon.push(pick);
  }
  const rec = t => {
    const thick = t.cds[1] > t.cds[0] ? t.cds : [t.start, t.start];
    const tags = [...t.tags].filter(g => /^(Ensembl_canonical|MANE_Select|MANE_Plus_Clinical|basic|CCDS)$/.test(g)).join(',');
    return { chrom: t.chr, start: t.start, end: t.end, rest: [t.id, 0, t.strand === '-' ? '-' : '+', thick[0], thick[1], 0, t.exons.length,
      t.exons.map(([a, b]) => b - a).join(',') + ',', t.exons.map(([a]) => a - t.start).join(',') + ',', t.gene, t.geneId, t.type, tags].join('\t') };
  };
  const rule = sawCanonTag ? 'Ensembl_canonical tag' : 'MANE Select, else longest coding part, else longest transcript';
  const sql = desc => `table lgvGenes\n"${desc}"\n(\nstring chrom; "Chromosome"\nuint chromStart; "Start"\nuint chromEnd; "End"\nstring name; "Transcript"\nuint score; "Unused"\nchar[1] strand; "Strand"\nuint thickStart; "Coding start"\nuint thickEnd; "Coding end"\nuint reserved; "Unused"\nint blockCount; "Exons"\nint[blockCount] blockSizes; "Exon lengths"\nint[blockCount] chromStarts; "Exon starts"\nstring geneName; "Gene"\nstring geneId; "Gene id"\nstring transcriptType; "Transcript type"\nstring tags; "Tags"\n)\n`;
  const opts = (desc) => ({ autoSql: sql(desc), fieldCount: 16, definedFieldCount: 12, itemsPerSlot: 512, blockSize: 256,
    extraIndex: [{ field: 3, key: r => strip(r.rest.split('\t')[0]).toUpperCase() }, { field: 12, key: r => r.rest.split('\t')[9].toUpperCase() }, { field: 13, key: r => strip(r.rest.split('\t')[10]).toUpperCase() }] });
  const src = path.basename(IN);
  fs.mkdirSync(OUT, { recursive: true });
  const chroms = Object.entries(SIZES);
  const A = writeBigBed(all.map(rec), chroms, opts(`LGV genes: all transcripts from ${src}`));
  const C = writeBigBed(canon.map(rec), chroms, opts(`LGV genes: one transcript per gene (${rule}) from ${src}`));
  fs.writeFileSync(path.join(OUT, 'ensembl.all.grch38.bb'), A);
  fs.writeFileSync(path.join(OUT, 'ensembl.canonical.grch38.bb'), C);
  console.log(`${lines} GTF lines; ${all.length} transcripts (${(A.length / 1048576).toFixed(1)} MB), ${canon.length} genes (${(C.length / 1048576).toFixed(1)} MB); one per gene by: ${rule}`);
})().catch(e => { console.error(e); process.exit(1); });
