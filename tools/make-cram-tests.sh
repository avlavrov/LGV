#!/bin/sh
# Makes test/cram/*.cram from demo/demo.bam. Needs samtools (built from github.com/samtools/samtools
# with htslib); the files are committed, so this is only needed to remake them.
set -e
ST=${SAMTOOLS:-samtools}; D=$(dirname "$0")/../demo; O=$(dirname "$0")/../test/cram
mkdir -p "$O"
$ST view -C -T "$D/demo.fa" -o "$O/v30.cram" --output-fmt-option version=3.0 "$D/demo.bam"
$ST view -C -T "$D/demo.fa" -o "$O/v30_bz_lzma.cram" --output-fmt-option version=3.0 --output-fmt-option use_bzip2=1 --output-fmt-option use_lzma=1 "$D/demo.bam"
$ST view -C -T "$D/demo.fa" -o "$O/v31.cram" --output-fmt-option version=3.1 "$D/demo.bam"
$ST view -C -T "$D/demo.fa" -o "$O/v31_archive.cram" --output-fmt-option version=3.1 --output-fmt-option archive "$D/demo.bam"
# gzip-only CRAM 3.0 variants: reference-based, no reference, embedded reference, small multi-reference slices
for v in "g30:" "g30_noref:--output-fmt-option no_ref=1" "g30_embed:--output-fmt-option embed_ref=1" "g30_small:--output-fmt-option seqs_per_slice=300 --output-fmt-option slices_per_container=4 --output-fmt-option multi_seq_per_slice=1"; do
  n=${v%%:*}; o=${v#*:}
  $ST view -C -T "$D/demo.fa" -o "$O/$n.cram" --output-fmt-option version=3.0 --output-fmt-option use_rans=0 $o "$D/demo.bam"
done
for f in "$O"/*.cram; do $ST index "$f"; done
# htslib's own SAM test cases (needs HTSLIB_TEST=<htslib>/test): sorted BAM → CRAM 3.0 → BAM decoded by htslib
if [ -n "$HTSLIB_TEST" ]; then
  H=$O/htslib; mkdir -p "$H"
  for f in "$HTSLIB_TEST"/auxf#*.sam "$HTSLIB_TEST"/c1#*.sam "$HTSLIB_TEST"/c2#*.sam "$HTSLIB_TEST"/ce#*.sam "$HTSLIB_TEST"/md#*.sam "$HTSLIB_TEST"/xx#*.sam; do
    b=$(basename "$f" .sam); n=$(echo "$b" | tr '#' '_'); p=${b%%#*}
    $ST sort -o "$H/tmp.bam" "$f" && $ST view -C -T "$HTSLIB_TEST/$p.fa" --output-fmt-option version=3.0 --output-fmt-option use_rans=0 -o "$H/$n.cram" "$H/tmp.bam" && $ST index "$H/$n.cram" \
      && $ST view -b -T "$HTSLIB_TEST/$p.fa" -o "$H/$n.bam" "$H/$n.cram" && $ST index "$H/$n.bam" || rm -f "$H/$n".*
    cp "$HTSLIB_TEST/$p.fa" "$HTSLIB_TEST/$p.fa.fai" "$H/"
  done
  rm -f "$H/tmp.bam"
fi
