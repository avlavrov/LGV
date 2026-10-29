#!/usr/bin/env python3
# Usage: python3 tools/make-mane.py <genome-nexus-importer>/data/ data/mane.grch38.tsv
# Builds data/mane.grch38.tsv for CGV from the MANE summary (NCBI/EMBL-EBI) and Ensembl 111 exon
# coordinates (both as exported by genome-nexus-importer). One line per transcript, 0-based half-open.
import gzip, json, sys
R=sys.argv[1]; OUT=sys.argv[2]
mane=[l.rstrip('\n').split('\t') for l in open(R+'common_input/MANE.GRCh38.v1.2.summary.txt') if not l.startswith('#')]
ids={m[7].split('.')[0] for m in mane}
ens={}
with gzip.open(R+'grch38_ensembl111/export/ensembl_biomart_transcripts.json.gz','rt') as f:
    for l in f:
        d=json.loads(l)
        if d['transcript_stable_id'] in ids: ens[d['transcript_stable_id']]=d
chrs={f'NC_0000{i:02d}':f'chr{i}' for i in range(1,23)}; chrs['NC_000023']='chrX'; chrs['NC_000024']='chrY'
out=[]; bad=0
for m in mane:
    acc=m[10].split('.')[0]
    if acc not in chrs: continue
    d=ens[m[7].split('.')[0]]
    ex=sorted((int(e['start'])-1,int(e['end'])) for e in d['exons'])
    utr=[(int(u['start'])-1,int(u['end'])) for u in d.get('utrs') or []]
    s0, e0 = int(m[11])-1, int(m[12])
    if ex[0][0]!=s0 or ex[-1][1]!=e0: bad+=1; continue
    # coding bases = exon bases outside UTRs
    pieces=[]
    for a,b in ex:
        segs=[(a,b)]
        for ua,ub in utr:
            n=[]
            for x,y in segs:
                if ub<=x or ua>=y: n.append((x,y)); continue
                if x<ua: n.append((x,ua))
                if ub<y: n.append((ub,y))
            segs=n
        pieces+=segs
    cs, ce = (min(p[0] for p in pieces), max(p[1] for p in pieces)) if pieces else (e0, e0)
    strand=m[13]
    assert strand in '+-'
    out.append('\t'.join([m[3], chrs[acc], strand, str(s0), str(e0-s0), str(cs-s0), str(ce-s0),
        ','.join(str(a-s0) for a,b in ex), ','.join(str(b-a) for a,b in ex),
        m[7], m[5], 'C' if m[9]=='MANE Plus Clinical' else 'S', m[4].replace('\t',' '), m[2]]))
order={f'chr{i}':i for i in range(1,23)}; order['chrX']=23; order['chrY']=24
out.sort(key=lambda l:(order[l.split('\t')[1]], int(l.split('\t')[3])))
hdr='#CGV MANE GRCh38 v1.2 (MANE summary from NCBI/EMBL-EBI; exons from Ensembl 111). 0-based. symbol chr strand start length cdsStart cdsEnd exonStarts exonLengths ensembl refseq status(S=Select,C=Plus Clinical) name hgnc\n'
open(OUT,'w').write(hdr+'\n'.join(out)+'\n')
print(len(out),'written; bad',bad)
