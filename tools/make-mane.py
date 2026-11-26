#!/usr/bin/env python3
# Usage: python3 tools/make-mane.py MANE.GRCh38.vX.summary.txt.gz MANE.GRCh38.vX.ensembl_genomic.gff.gz data/mane.grch38.tsv
# Builds data/mane.grch38.tsv for LGV from one MANE release (NCBI / EMBL-EBI,
# ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/): the summary gives the transcripts, their RefSeq /
# Ensembl ids, status and gene names; the release's own Ensembl GFF3 gives the exons and UTRs.
# One line per transcript on chr1-22, X, Y (transcripts on patch / alt sequences are left out),
# 0-based half-open. The coding part is the exons minus the UTRs (so it includes the stop codon).
import gzip, re, sys
SUM, GFF, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
ver = re.search(r'v(\d+\.\d+)', SUM).group(1)
op = lambda p: gzip.open(p, 'rt') if p.endswith('.gz') else open(p)
rows = [l.rstrip('\n').split('\t') for l in op(SUM) if not l.startswith('#')]
want = {r[7]: r for r in rows}  # Ensembl_nuc (with version) -> summary row
PRIMARY = {f'chr{i}' for i in range(1, 23)} | {'chrX', 'chrY'}
tx = {}  # transcript id -> {chr, strand, exons[], utrs[]}
for l in op(GFF):
    if l.startswith('#'): continue
    f = l.rstrip('\n').split('\t')
    if len(f) < 9 or f[0] not in PRIMARY: continue
    typ = f[2]
    if typ not in ('exon', 'five_prime_UTR', 'three_prime_UTR'): continue
    attrs = dict(kv.split('=', 1) for kv in f[8].split(';') if '=' in kv)
    t = attrs.get('transcript_id') or attrs.get('Parent')
    if t not in want: continue
    d = tx.setdefault(t, {'chr': f[0], 'strand': f[6], 'exons': [], 'utrs': []})
    (d['exons'] if typ == 'exon' else d['utrs']).append((int(f[3]) - 1, int(f[4])))
out = []; missing = 0
for t, m in want.items():
    d = tx.get(t)
    if not d or not d['exons']: missing += 1; continue
    ex = sorted(d['exons']); s0, e0 = ex[0][0], ex[-1][1]
    pieces = []
    for a, b in ex:
        segs = [(a, b)]
        for ua, ub in d['utrs']:
            n = []
            for x, y in segs:
                if ub <= x or ua >= y: n.append((x, y)); continue
                if x < ua: n.append((x, ua))
                if ub < y: n.append((ub, y))
            segs = n
        pieces += segs
    cs, ce = (min(p[0] for p in pieces), max(p[1] for p in pieces)) if pieces else (e0, e0)
    assert d['strand'] in '+-'
    out.append('\t'.join([m[3], d['chr'], d['strand'], str(s0), str(e0 - s0), str(cs - s0), str(ce - s0),
        ','.join(str(a - s0) for a, b in ex), ','.join(str(b - a) for a, b in ex),
        m[7], m[5], 'C' if m[9] == 'MANE Plus Clinical' else 'S', m[4].replace('\t', ' '), m[2]]))
order = {f'chr{i}': i for i in range(1, 23)}; order['chrX'] = 23; order['chrY'] = 24
out.sort(key=lambda l: (order[l.split('\t')[1]], int(l.split('\t')[3])))
hdr = f'#LGV MANE GRCh38 v{ver} (MANE summary and Ensembl GFF3 of the same release, NCBI/EMBL-EBI). 0-based. symbol chr strand start length cdsStart cdsEnd exonStarts exonLengths ensembl refseq status(S=Select,C=Plus Clinical) name hgnc\n'
open(OUT, 'w').write(hdr + '\n'.join(out) + '\n')
print(len(out), 'written;', missing, 'not on chr1-22/X/Y (patches, alt, chrM) or without exons')
