"""Writes TDF and bigWig test files for LGV: small ones for tools/test-signal.js (test/signal) and,
with "big", large ones for speed checks (test/signal/big, not committed).

TDF follows IGV's TDFWriter layout (version 4, zlib tiles, zoom levels 0-7 of 700 bins per tile,
raw dataset per contig). bigWig follows the UCSC layout (bbiFile.h): header, zoom headers, total
summary, chromosome B+ tree, data sections (bedGraph type) with an R tree, zoom levels with their own
R trees. Values are synthetic.

Run: python3 tools/make-signal-tests.py [big]
"""
import os, sys, struct, zlib, math
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'test', 'signal')
HG38 = {'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555, 'chr5': 181538259,
        'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717, 'chr10': 133797422,
        'chr11': 135086622, 'chr12': 133275309, 'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
        'chr16': 90338345, 'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
        'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895}

# ---------------- TDF ----------------
def cstr(s): return s.encode() + b'\0'

class Tdf:
    def __init__(self, path, track_type, track_line, names, wfs=('mean',), genome='hg38'):
        self.f = open(path, 'wb'); self.pos = 0; self.ds = {}; self.groups = {'/': {}}; self.names = names; self.wfs = wfs
        self.w(b'TDF4' + struct.pack('<i', 4) + struct.pack('<qi', 0, 0))
        h = struct.pack('<i', len(wfs)) + b''.join(cstr(x) for x in wfs) + cstr(track_type)
        tl = track_line.encode()[:1023].ljust(1023, b' ') + b'\0'
        h += tl + struct.pack('<i', len(names)) + b''.join(cstr(n) for n in names) + cstr(genome) + struct.pack('<i', 1)
        self.w(struct.pack('<i', len(h)) + h)
    def w(self, b): self.f.write(b); self.pos += len(b)
    def tile(self, ds, i, data):
        c = zlib.compress(data, 6); self.ds[ds]['tiles'][i] = (self.pos, len(c)); self.w(c)
    def dataset(self, name, tile_width, n):
        self.ds[name] = {'tw': tile_width, 'tiles': [(-1, 0)] * n}
    def close(self):
        di, gi = [], []
        for name, d in self.ds.items():
            p = self.pos
            b = struct.pack('<i', 0) + cstr('FLOAT') + struct.pack('<f', d['tw']) + struct.pack('<i', len(d['tiles']))
            b += b''.join(struct.pack('<qi', a, s) for a, s in d['tiles'])
            self.w(b); di.append((name, p, len(b)))
        for name, g in self.groups.items():
            p = self.pos; b = struct.pack('<i', len(g)) + b''.join(cstr(k) + cstr(v) for k, v in g.items())
            self.w(b); gi.append((name, p, len(b)))
        p = self.pos
        b = struct.pack('<i', len(di)) + b''.join(cstr(n) + struct.pack('<qi', a, s) for n, a, s in di)
        b += struct.pack('<i', len(gi)) + b''.join(cstr(n) + struct.pack('<qi', a, s) for n, a, s in gi)
        self.w(b); self.f.seek(8); self.f.write(struct.pack('<qi', p, len(b))); self.f.close()

def bed_tile(starts, ends, vals):  # vals: (ntracks, n)
    n = len(starts)
    return cstr('bed') + struct.pack('<i', n) + starts.astype('<i4').tobytes() + ends.astype('<i4').tobytes() + struct.pack('<i', vals.shape[0]) + vals.astype('<f4').tobytes()

def fixed_tile(start, span, vals):
    return cstr('fixedStep') + struct.pack('<iif', vals.shape[1], int(start), float(span)) + vals.astype('<f4').tobytes()

def write_tdf(path, chroms, data, track_type='OTHER', track_line='', names=('track',), raw_tile=1000000, wfs=('mean',)):
    """data: {chr: (starts, ends, vals[ntracks, n])}, sorted by start"""
    t = Tdf(path, track_type, track_line, list(names), wfs)
    allv = np.concatenate([v[2].ravel() for v in data.values()]) if data else np.zeros(1)
    allv = allv[~np.isnan(allv)]
    t.groups['/'] = {'genome': 'hg38', 'maxZoom': '7', 'chromosomes': ','.join(data.keys()),
                     'Maximum': '%g' % allv.max(), 'Minimum': '%g' % allv.min(), 'Mean': '%g' % allv.mean()}
    for chr_, (s, e, v) in data.items():
        L = chroms[chr_]
        nt = math.ceil(L / raw_tile); t.dataset(f'/{chr_}/raw', raw_tile, nt)
        ti = s // raw_tile
        for i in np.unique(ti):
            m = ti == i; t.tile(f'/{chr_}/raw', int(i), bed_tile(s[m], e[m], v[:, m]))
        mid = (s + e) // 2
        for z in range(8):
            ntz = 2 ** z; tw = L / ntz; bw = tw / 700
            for wf in wfs:
                t.dataset(f'/{chr_}/z{z}/{wf}', tw, ntz)
            b = np.floor(mid / bw).astype(np.int64)
            for i in np.unique(b // 700):
                m = (b // 700) == i
                out = {wf: np.full((v.shape[0], 700), np.nan, dtype=np.float32) for wf in wfs}
                loc = b[m] - i * 700
                for tr in range(v.shape[0]):
                    vv = v[tr, m]
                    for wf in wfs:
                        o = out[wf][tr]
                        order = np.argsort(loc, kind='stable'); lk = loc[order]; vk = vv[order]
                        u, idx = np.unique(lk, return_index=True)
                        if wf == 'mean': red = np.add.reduceat(vk, idx) / np.diff(np.append(idx, len(vk)))
                        elif wf == 'min': red = np.minimum.reduceat(vk, idx)
                        else: red = np.maximum.reduceat(vk, idx)
                        o[u] = red
                for wf in wfs:
                    t.tile(f'/{chr_}/z{z}/{wf}', int(i), fixed_tile(i * tw, bw, out[wf]))
    t.close()

# ---------------- bigWig ----------------
def write_bigwig(path, chroms, data, items_per_block=1024, zooms=(), le=True):
    """data: {chr: (starts, ends, vals)}; bedGraph sections; R tree with one leaf per block, all in one
    leaf node level under a root (enough for tests: block size large)."""
    E = '<' if le else '>'
    names = list(chroms.keys()); ids = {n: i for i, n in enumerate(names)}
    f = open(path, 'wb')
    hdr_size = 64 + 24 * len(zooms)
    f.write(b'\0' * (hdr_size + 40))
    # chromosome B+ tree: one leaf node
    ct_off = f.tell(); key = max(len(n) for n in names)
    f.write(struct.pack(E + 'IIIIQQ', 0x78CA8C91, len(names), key, 8, len(names), 0))
    f.write(struct.pack(E + 'BBH', 1, 0, len(names)))
    for n in sorted(names): f.write(n.encode().ljust(key, b'\0') + struct.pack(E + 'II', ids[n], chroms[n]))
    # data
    data_off = f.tell(); blocks = []; maxbuf = 0
    nsec = 0; f.write(struct.pack(E + 'Q', 0))
    tot = [0, math.inf, -math.inf, 0.0, 0.0]
    for n in names:
        if n not in data: continue
        s, e, v = data[n]
        for k in range(0, len(s), items_per_block):
            ss, ee, vv = s[k:k + items_per_block], e[k:k + items_per_block], v[k:k + items_per_block]
            sec = struct.pack(E + 'IIIIIBBH', ids[n], int(ss[0]), int(ee[-1]), 0, 0, 1, 0, len(ss))
            rec = np.zeros(len(ss), dtype=[('s', E + 'u4'), ('e', E + 'u4'), ('v', E + 'f4')]); rec['s'] = ss; rec['e'] = ee; rec['v'] = vv
            raw = sec + rec.tobytes(); maxbuf = max(maxbuf, len(raw)); c = zlib.compress(raw, 6)
            blocks.append((ids[n], int(ss[0]), ids[n], int(ee[-1]), f.tell(), len(c))); f.write(c); nsec += 1
            w = (ee - ss).astype(np.float64)
            tot[0] += int(w.sum()); tot[1] = min(tot[1], float(vv.min())); tot[2] = max(tot[2], float(vv.max())); tot[3] += float((vv * w).sum()); tot[4] += float((vv.astype(np.float64) ** 2 * w).sum())
    def rtree(blocks):
        off = f.tell(); n = len(blocks)
        first, last = blocks[0], blocks[-1]
        f.write(struct.pack(E + 'IIQIIIIQII', 0x2468ACE0, max(n, 1), n, first[0], first[1], last[2], last[3], f.tell(), 1, 0))
        # root (inner) -> leaves of up to 256 items
        groups = [blocks[i:i + 256] for i in range(0, n, 256)]
        root = f.tell(); f.write(struct.pack(E + 'BBH', 0, 0, len(groups)))
        child_at = f.tell(); f.write(b'\0' * 24 * len(groups))
        offs = []
        for g in groups:
            offs.append(f.tell()); f.write(struct.pack(E + 'BBH', 1, 0, len(g)))
            for b in g: f.write(struct.pack(E + 'IIIIQQ', *b))
        end = f.tell(); f.seek(child_at)
        for g, o in zip(groups, offs): f.write(struct.pack(E + 'IIIIQ', g[0][0], g[0][1], g[-1][2], g[-1][3], o))
        f.seek(end); return off
    idx_off = rtree(blocks)
    zh = []
    for red in zooms:
        zdata = f.tell(); zb = []
        for n in names:
            if n not in data: continue
            s, e, v = data[n]; b = s // red
            u, idx = np.unique(b, return_index=True); cnt = np.diff(np.append(idx, len(b)))
            w = (e - s).astype(np.float64); vw = v.astype(np.float64)
            recs = np.zeros(len(u), dtype=[('c', E + 'u4'), ('s', E + 'u4'), ('e', E + 'u4'), ('n', E + 'u4'), ('mn', E + 'f4'), ('mx', E + 'f4'), ('sum', E + 'f4'), ('sq', E + 'f4')])
            recs['c'] = ids[n]; recs['s'] = u * red; recs['e'] = np.minimum((u + 1) * red, chroms[n])
            recs['n'] = np.add.reduceat(w, idx); recs['mn'] = np.minimum.reduceat(vw, idx); recs['mx'] = np.maximum.reduceat(vw, idx)
            recs['sum'] = np.add.reduceat(vw * w, idx); recs['sq'] = np.add.reduceat(vw * vw * w, idx)
            for k in range(0, len(recs), items_per_block):
                r = recs[k:k + items_per_block]; raw = r.tobytes(); maxbuf = max(maxbuf, len(raw)); c = zlib.compress(raw, 6)
                zb.append((ids[n], int(r['s'][0]), ids[n], int(r['e'][-1]), f.tell(), len(c))); f.write(c)
        zh.append((red, zdata, rtree(zb)))
    f.seek(0)
    f.write(struct.pack(E + 'IHHQQQHHQQIQ', 0x888FFC26, 4, len(zooms), ct_off, data_off, idx_off, 0, 0, 0, hdr_size, maxbuf, 0))
    for red, d, i in zh: f.write(struct.pack(E + 'IIQQ', red, 0, d, i))
    f.write(struct.pack(E + 'Qdddd', tot[0], tot[1], tot[2], tot[3], tot[4]))
    f.seek(data_off); f.write(struct.pack(E + 'Q', nsec)); f.close()

# ---------------- synthetic data ----------------
rng = np.random.default_rng(7)
def allele_freqs(L, n, cluster=True):
    """variant positions in clusters (like exome targets), allele fractions: het ~0.5, hom ~1"""
    if cluster:
        centers = np.sort(rng.integers(0, L - 5000, max(1, n // 6)))
        pos = np.sort((centers[rng.integers(0, len(centers), n)] + rng.integers(0, 3000, n)).astype(np.int64))
    else:
        pos = np.sort(rng.integers(0, L - 1, n).astype(np.int64))
    pos = np.unique(pos)
    k = rng.random(len(pos))
    v = np.where(k < 0.62, np.clip(rng.normal(0.5, 0.08, len(pos)), 0.02, 0.98), np.where(k < 0.97, np.clip(rng.normal(0.99, 0.01, len(pos)), 0.9, 1), rng.uniform(0.05, 0.25, len(pos))))
    return pos, pos + 1, v.astype(np.float32)

def coverage(L, start, n, step=25, depth=30):
    s = np.arange(start, start + n * step, step, dtype=np.int64); s = s[s < L - step]
    x = np.cumsum(rng.normal(0, 0.6, len(s))); v = np.clip(depth + 6 * np.sin(np.arange(len(s)) / 400) + x - np.convolve(x, np.ones(200) / 200, 'same') + rng.normal(0, 2, len(s)), 0, None)
    return s, s + step, v.astype(np.float32)

if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    small = {'chr7': HG38['chr7'], 'chr17': HG38['chr17']}
    # allele frequencies at variant positions, no track line, type OTHER: the 0-1 rule makes points
    d = {}
    for c in small: s, e, v = allele_freqs(small[c], 3000); d[c] = (s, e, v[None, :])
    write_tdf(os.path.join(OUT, 'af_sample.tdf'), small, d)
    # dense 0-1 values (a GC fraction in 100 bp bins): bars
    d = {}
    for c in small:
        s = np.arange(55_000_000, 55_300_000, 100, dtype=np.int64); v = np.clip(0.45 + 0.1 * np.sin(s / 7000) + rng.normal(0, 0.03, len(s)), 0, 1).astype(np.float32)
        d[c] = (s, s + 100, v[None, :])
    write_tdf(os.path.join(OUT, 'dense01_sample.tdf'), small, d)
    print('small files written')
    if len(sys.argv) > 1 and sys.argv[1] == 'big':
        B = os.path.join(OUT, 'big'); os.makedirs(B, exist_ok=True)
        # whole-genome allele frequencies, ~3.5 million variants (genome data), one file of 3 tracks
        d = {}
        for c, L in HG38.items():
            s, e, v = allele_freqs(L, int(L / 900), cluster=False)
            v2 = np.clip(v + rng.normal(0, 0.03, len(v)), 0, 1).astype(np.float32); v3 = np.clip(1 - v + rng.normal(0, 0.03, len(v)), 0, 1).astype(np.float32)
            d[c] = (s, e, np.vstack([v, v2, v3]))
        write_tdf(os.path.join(B, 'genome_af_3tracks.tdf'), HG38, d, names=('proband', 'mother', 'father'), wfs=('mean', 'min', 'max'))
        print('genome_af_3tracks.tdf')
        # coverage of chr1 at 25 bp windows (~10 million values) as TDF and bigWig
        s, e, v = coverage(HG38['chr1'], 0, HG38['chr1'] // 25)
        write_tdf(os.path.join(B, 'chr1_coverage.tdf'), {'chr1': HG38['chr1']}, {'chr1': (s, e, v[None, :])}, track_type='COVERAGE', raw_tile=500000)
        print('chr1_coverage.tdf')
        write_bigwig(os.path.join(B, 'chr1_coverage.bw'), {'chr1': HG38['chr1']}, {'chr1': (s, e, v)}, zooms=(100, 400, 1600, 6400, 25600, 102400, 409600))
        print('chr1_coverage.bw')
        # exome-like allele frequencies as bigWig, both byte orders
        d = {c: allele_freqs(L, int(L / 15000)) for c, L in HG38.items()}
        write_bigwig(os.path.join(B, 'exome_af.bw'), HG38, d, zooms=(1000, 10000, 100000))
        write_bigwig(os.path.join(B, 'exome_af_be.bw'), HG38, d, zooms=(1000, 10000, 100000), le=False)
        print('exome_af.bw, exome_af_be.bw')
