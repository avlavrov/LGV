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
for f in "$O"/*.cram; do $ST index "$f"; done
