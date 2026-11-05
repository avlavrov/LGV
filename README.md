# LGV — Lavrov Genome Viewer

A fast, light and secure viewer for read alignments (BAM, CRAM) in human genetics. Built with
Claude, an AI model by Anthropic.
It runs in the browser with no server and no third-party dependencies.

**Status:** MVP, version 0.9.0. (Formerly CGV.)

**User guide:** [docs/USER_GUIDE.md](docs/USER_GUIDE.md) (Russian: [docs/USER_GUIDE.ru.md](docs/USER_GUIDE.ru.md)).
Project rules: [docs/PROJECT_RULES.md](docs/PROJECT_RULES.md). Plans: [docs/TODO.md](docs/TODO.md).

## Run

Open `dist/lgv.html` in Chrome, Edge, Firefox or Safari. It starts on the GRCh38 reference
(synthetic demo data when offline). To view your own data, press **Open**:

- **From disk**: choose one or more `.bam` files with their indexes (you can add `.fa` +
  `.fai`), choose a whole folder, or drop files or a folder onto the window. Each `.bam` is
  paired with its index by name: `sample.bam.bai` or `sample.bai`. A web page can read only
  the files it is given, so the index has to be chosen too (selecting the folder is the
  easiest way); a folder with several samples first lists them to pick from. Files are read
  in pieces straight from disk and never uploaded.
- **By link**: one or more HTTPS links to BAM files, one per line. The index is looked for
  next to each BAM (`.bam.bai`, then `.bai`) unless given. The server must allow CORS and
  Range requests. An access token can be given; it is not saved.
- **Several samples** are shown as tracks, one under another, under one ruler and reference,
  all at the same position. Each track has its own coverage, reads and scrolling and is
  closed with its ×. Sorting applies to all of them. Opening your own data replaces the demo.
- **Reference**: the reference sets the genome, as in IGV. Loading one switches the genome
  and closes open samples; the reference itself can be browsed right away. Samples opened
  afterwards are matched to it: contig names are compared without case, `chr` prefix or
  version suffix (and with NCBI header fields), then remaining contigs pair up by a length
  only one sequence on each side has. A sample with no matching contig is refused with
  both lists shown. One `.fa` / `.fasta` file is enough: without a `.fai` LGV indexes it in
  the browser in one pass (plain, bgzip or gzip); with `.fai` (and `.gzi` for `.fa.gz`) it
  opens at once. Public GRCh38, GRCh37 and T2T-CHM13 v2.0 copies (hosted by the IGV team
  on Amazon S3) are offered as presets. Without a reference, the sample's header gives the
  coordinates and reference bases are rebuilt from MD tags where reads cover them.

The interface is in English; Russian is available under **View → Screen → Language**.
More languages are one dictionary in `src/i18n.js`.

## Stage 0 coverage of the MVP spec

| Requirement | State |
| --- | --- |
| F1 BAM by HTTPS link (Range) | done; CRAM in progress: CRAM 3.0 / 3.1 records are read and match BAM read for read; CRAM 3.0 and 3.1 open, including 3.1 archive mode, decoded in Web Workers; each slice's reference MD5 is checked against the loaded reference; every CRAM block codec is read, bzip2 and lzma by LGV's own decoders |
| F2 Local BAM + BAI | done: several samples at once as tracks, index paired by name, folder or drag and drop; in Chrome and Edge choosing just the .bam / .cram is enough — LGV asks once to read its folder (read only, this session) and takes the index from there; with many samples each keeps a minimum height (View → Coverage, 240 px by default) and the whole stack scrolls under the ruler: wheel over coverage or a header, wheel past the end of a sample's reads, the scrollbar on the right, PageUp / PageDown |
| F3 Reference (GRCh38, GRCh37, T2T-CHM13) | done: files or link, plain or bgzip, public presets |
| F4 Access token / presigned URL | done |
| F5 Go to coordinate | done, contig names case-insensitive (chr1 = Chr1 = 1); drag on the ruler to zoom to a region, click it to center; **gene search**: a gene symbol, an ENST or NM id (MANE), with suggestions while typing; several genes or loci separated by spaces open side by side |
| Chromosome map | across the top: the chromosome of each panel with GRCh38 bands and centromere (plain bar for other genomes), a red box where each panel is (numbered like the panel headers); click or drag to move the active panel |
| Strand | the 3' end of every read is pointed, as in IGV, so the strand shows under any coloring |
| Annotation tracks | BED, bedGraph, VCF, GTF, GFF3 (plain, gzip or bgzip) from disk or a link: bands under the gene track; click a feature for its details |
| Signal tracks | TDF and bigWig, read by region (only what a view needs); every track of a TDF as its own row; drawn by the file's UCSC track line (bars, points, windowing function, scale, log, smoothing, colors); a TDF that says nothing is drawn as points when its values are 0–1 at separate positions (allele fractions), else as bars; right click to change |
| Genes | **MANE track** (GRCh38): MANE Select and Plus Clinical transcripts under the reference — coding exons thick, UTRs thin, strand arrows, exon numbers when zoomed in; click a transcript for its ids, span and the exon or intron at that position. Bundled with LGV (MANE v1.2), works offline |
| Start | GRCh38 opens by default (View → "At start open": hg38, hg19, T2T, demo or nothing); without a connection the demo opens instead |
| F6 No hard window limit | windows load on their own up to 64 MB of compressed data, larger on request; downsampling keeps deep regions fast; isolated sequencing errors fade when zoomed out |
| F8 Coverage with allele highlight | done, threshold adjustable; click the coverage for counts at that base as in IGV: A/C/G/T/N with strands, DEL, INS, mean base quality — from every read, before downsampling |
| F9 Smooth read height | done, 1–24 px |
| F10 Coloring + MAPQ as a transparency gradient | right click on reads: strand, pair (orientation and insert size), tag; MAPQ works on top of any of them |
| F11 Group and sort | groups (View): strand, first/second, pair orientation, tag; sort by (right click, at the clicked position): base, MAPQ, strand, insert size |
| F12 Filters | MAPQ, duplicates, secondary, supplementary, QC fail; downsampling |
| Pictures and links | PNG of the view; a link that opens the same view; open LGV from another system with `#ref=…&file=…&locus=…` |
| Navigation | ← → back / forward between places jumped to; − + zoom all panels at once |
| F13 Mate | view as pairs (right click), read details + go to mate; **split screen** as in IGV: right click on a read → "View mate region in split screen" (or the button in the read details) opens a panel centered on the mate, the read and its mate outlined in every panel; several loci typed in the locus field separated by spaces open one panel each; each panel has its own ruler, scroll and zoom, × in its header closes it; at least 6 panels, more as the width allows (60 px each) |
| F14 Read details | IGV-style popup: read, mapping, clipping, mate and tags, then Location and `Base = A @ QV 14` for the clicked base; one read or read + mate side by side; copy; dragged by its top row, then opens where it was left (until the page is closed). Right click on a read: copy its sequence, left or right soft clip, or the pair as FASTA |
| F18 Display profile as JSON | kept in the browser, copied and applied as JSON |
| F23 Screen space for data | one bar of controls, hidden with H |

## How it is built

```
src/core.js    data sources (file, URL), BGZF, BAI, BAM, FASTA (+ .gzi), coverage, layout
src/cram.js    CRAM 3.0 / 3.1: containers, blocks, SAM header, .crai, records
src/compress.js   bzip2 and xz/LZMA2 decoders, written for LGV
src/htscodecs.js  CRAM block codecs adapted from htscodecs: rANS 4x8 / 4x16, arithmetic, FQZcomp, names (BSD; changes in SECURITY.md)
test/cram/     CRAM test files made with samtools (tools/make-cram-tests.sh); htslib/ holds htslib's own test cases
src/app.js     WebGL2 read rendering, Canvas2D ruler and coverage, UI, display profile
src/i18n.js    interface strings: English base, Russian overlay
src/index.html markup and styles
tools/         demo data generator, build, tests
```

BGZF is inflated with the browser's built-in `DecompressionStream('deflate-raw')`; the CRC32 and
length of every block are checked, and every size read from a file is checked before use (see SECURITY.md). There are
no third-party libraries: a smaller attack surface and nothing to keep patched. The local
build (`dist/lgv.html`) carries a strict Content Security Policy: no remote scripts, styles
or fonts, and network access only for data the user points to.

## Development

```
node tools/gen-demo.js                      # demo data in demo/
node tools/test-core.js                     # core checks on the demo data
node tools/test-cram.js                     # CRAM checks against the same data as BAM, with fuzzing
node tools/test-fuzz.js                     # fuzzing of BGZF, BAM, BAI, FASTA/.fai/.gzi, locus field, profile
node tools/test-signal.js                   # TDF and bigWig on files from IGV, igvtools, UCSC tools, bigtools; fuzzing
python3 tools/make-signal-tests.py [big]    # synthetic TDF / bigWig test files (big: large ones for speed)
node tools/audit.js                         # security audit of our own code (see SECURITY.md)
node tools/build.js                         # dist/lgv.html and dist/artifact.html
CGV_L=1000000 node tools/gen-demo.js big    # 1 Mb at 30× for speed tests
```

## Deep data

Targeted panels reach thousands of reads per base. LGV handles them the way IGV does:

- **Downsampling while reading:** at most 100 reads are kept per 50 bp window (both
  adjustable under View → Filters). Coverage still counts every read. Marks above the
  reads show where reads were left out.
- **Streaming reads:** the BAM is read in slabs and reading stops as soon as reads pass the
  window, so the rest of a 16 kb index bin is never fetched or inflated.
- **Block cache:** record start positions in BGZF blocks are remembered, so revisiting or
  panning near a region skips straight to the right place.
- **Sparse coverage and reference:** coverage is kept in 64 kb blocks only where reads are,
  and reference bases are read only under reads (plus the view when it is under 200 kb).
  A 46 Mb window with a few reads uses a few blocks and a few kb of FASTA, not gigabytes.
- **Depth-only coverage for wide views:** above 1 Mb the coverage is counted from CIGAR ranges
  instead of base by base (allele colors cannot be seen at that scale); zooming in below 1 Mb
  loads the window again in detail.
- **Web Workers:** samples are read and inflated in background workers, several in parallel, so
  the page stays responsive while a window loads. Where workers cannot start, reading falls back
  to the page.
- **Load limit by size:** a window loads on its own when the index points to at most 64 MB
  of compressed alignments; larger windows ask first.

## Read layout

- **Stable:** every read keeps its row while the layout lives — panning, reloading the
  window, recoloring, changing read height. Only reads seen for the first time are placed.
  A new layout is made only when something that decides placement changes: contig,
  grouping, sort, filters or downsampling. The same reads and settings always give the
  same picture.
- **Compact:** a new read goes into the lowest row where it fits (2 bp gap), so coverage
  dips fill from the top. Without grouping, the scroll range follows the rows the visible
  part uses, not the deepest spot in the loaded margin.
- **Sorting is a layout, not a row shuffle:** reads covering the sort position take the top
  rows in sort order; all other reads pack around them. The sort stays exact at that
  position and the rest stays compact; panning afterwards keeps it.
- **Pairs:** with *View as pairs* (right click on the reads), both mates share a row, joined by
  a line; where they overlap the overlap is shaded. A pair is placed as one item whose span
  comes from the mate position and its MC tag, so it does not change when the mate loads
  later; a mate not loaded yet is shown by a line to where it lies. Pairs spanning more than
  10 kb (`pairs.maxSpan` in the profile), mates on other contigs, unmapped mates and
  secondary or supplementary alignments are shown as single reads. Downsampling picks reads
  independently, so in deep downsampled regions some pairs show one mate.
- **Soft clips:** shown by default (right click → *Show soft-clipped bases*, or View → Reads):
  clipped bases continue the read in a lighter shade, bases that differ from the reference
  are colored, and the layout leaves room for them. Hidden, a clip is marked by a tick.
- **Downsampling is query-independent:** reads kept in a window never depend on where a
  query began, so reloads do not swap reads.

## Measurements (local disk, headless Chromium without a GPU)

| Data | Window | Reads in window | Load | Build view |
| --- | --- | --- | --- | --- |
| 30× | 2 kb | 1,339 | 26 ms | 7 ms |
| 30× | 1 Mb | 200,024 | 1,034 ms | 801 ms |
| 3000× panel | 32 bp, first visit | 4,527 | 491 ms | 5 ms |
| 3000× panel | 32 bp, nearby | 4,538 | 129 ms | 1 ms |
| 3000× panel | 32 bp, revisit | 4,315 | 64 ms | 1 ms |
| 3000× panel | 20 kb | 784,176 | 2,725 ms | 267 ms |

Over a network, the server's latency adds to the load time.

```
CGV_L=40000 CGV_DEPTH=3000 node tools/gen-demo.js deep   # deep panel-like data for speed tests
```

## License

Not chosen yet.
