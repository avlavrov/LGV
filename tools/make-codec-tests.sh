#!/bin/sh
# Makes test/codecs/: htscodecs' own encoders (tests/arith_dynamic, tests/rans4x16pr from
# github.com/samtools/htscodecs, built locally; HTSCODECS_TESTS=<path>) in every mode, and Python's
# bz2 / lzma for bzip2 and xz. The files are committed; this is only needed to remake them.
set -e
O=$(dirname "$0")/../test/codecs; mkdir -p "$O"
samtools view "$(dirname "$0")/../demo/demo.bam" | head -400 | cut -f10,11 | tr '\t' '\n' > "$O/in.txt"
for o in 0 1 4 5 8 9 32 64 65 128 129 192 193; do "$HTSCODECS_TESTS/arith_dynamic" -o $o "$O/in.txt" > "$O/arith.o$o"; done
for o in 0 1 4 5 8 9 64 65 128 129 192 193; do "$HTSCODECS_TESTS/rans4x16pr" -o $o "$O/in.txt" > "$O/r16.o$o"; done
python3 - "$O" <<'PY'
import bz2, lzma, os, random, sys
o = sys.argv[1]; random.seed(1)
data = open(o + '/in.txt', 'rb').read()
rnd = bytes(random.getrandbits(8) for _ in range(60000))
for name, d in (('text', data), ('random', rnd)):
    open(f'{o}/{name}.raw', 'wb').write(d)
    for lvl in (1, 9): open(f'{o}/{name}.bz{lvl}', 'wb').write(bz2.compress(d, lvl))
    for chk, tag in ((lzma.CHECK_CRC32, 'crc32'), (lzma.CHECK_CRC64, 'crc64'), (lzma.CHECK_NONE, 'none')):
        open(f'{o}/{name}.xz-{tag}', 'wb').write(lzma.compress(d, format=lzma.FORMAT_XZ, check=chk))
PY
