# CGV — Convenient, custom, genomic viewer

A fast, light and secure viewer for read alignments (BAM) in human genetics.
It runs in the browser with no server and no third-party dependencies.

**Status:** stage 0, technical prototype.

## Run

Open `dist/cgv.html` in Chrome, Edge, Firefox or Safari. It starts on synthetic demo
data. To view your own data, press **Open**:

- **From disk**: choose a `.bam` and its `.bai` (you can add `.fa` + `.fai`), or drop the
  files onto the window. Files are read in pieces straight from disk and never uploaded.
- **By link**: an HTTPS link to a BAM. The server must allow CORS and Range requests.
  An access token can be given; it is not saved.
- **Reference**: the reference sets the genome, as in IGV. Loading one switches the genome
  and closes open samples; the reference itself can be browsed right away. Samples opened
  afterwards are matched to it: contig names are compared without case, `chr` prefix or
  version suffix (and with NCBI header fields), then remaining contigs pair up by a length
  only one sequence on each side has. A sample with no matching contig is refused with
  both lists shown. One `.fa` / `.fasta` file is enough: without a `.fai` CGV indexes it in
  the browser in one pass (plain, bgzip or gzip); with `.fai` (and `.gzi` for `.fa.gz`) it
  opens at once. Public GRCh38, GRCh37 and T2T-CHM13 v2.0 copies (hosted by the IGV team
  on Amazon S3) are offered as presets. Without a reference, the sample's header gives the
  coordinates and reference bases are rebuilt from MD tags where reads cover them.

The interface is in English; Russian is available under **View → Screen → Language**.
More languages are one dictionary in `src/i18n.js`.

## Stage 0 coverage of the MVP spec

| Requirement | State |
| --- | --- |
| F1 BAM by HTTPS link (Range) | done; CRAM is the next step |
| F2 Local BAM + BAI | done, including drag and drop |
| F3 Reference (GRCh38, GRCh37, T2T-CHM13) | done: files or link, plain or bgzip, public presets |
| F4 Access token / presigned URL | done |
| F5 Go to coordinate | done, contig names case-insensitive (chr1 = Chr1 = 1); drag on the ruler to zoom to a region, click it to center; gene names come with the gene track |
| F6 No hard window limit | windows load on their own up to 64 MB of compressed data, larger on request; downsampling keeps deep regions fast; isolated sequencing errors fade when zoomed out |
| F8 Coverage with allele highlight | done, threshold adjustable |
| F9 Smooth read height | done, 1–24 px |
| F10 Coloring + MAPQ as a transparency gradient | right click on reads: strand, pair (orientation and insert size), tag; MAPQ works on top of any of them |
| F11 Group and sort | groups (View): strand, first/second, pair orientation, tag; sort by (right click, at the clicked position): base, MAPQ, strand, insert size |
| F12 Filters | MAPQ, duplicates, secondary, supplementary, QC fail; downsampling |
| F13 Mate | view as pairs (right click), read details + go to mate; several regions on screen is the next step |
| F14 Read details | IGV-style popup: read, mapping, clipping, mate and tags, then Location and `Base = A @ QV 14` for the clicked base; one read or read + mate side by side; copy |
| F18 Display profile as JSON | kept in the browser, copied and applied as JSON |
| F23 Screen space for data | one bar of controls, hidden with H |

## How it is built

```
src/core.js    data sources (file, URL), BGZF, BAI, BAM, FASTA (+ .gzi), coverage, layout
src/app.js     WebGL2 read rendering, Canvas2D ruler and coverage, UI, display profile
src/i18n.js    interface strings: English base, Russian overlay
src/index.html markup and styles
tools/         demo data generator, build, tests
```

BGZF is inflated with the browser's built-in `DecompressionStream('deflate-raw')`. There are
no third-party libraries: a smaller attack surface and nothing to keep patched. The local
build (`dist/cgv.html`) carries a strict Content Security Policy: no remote scripts, styles
or fonts, and network access only for data the user points to.

## Development

```
node tools/gen-demo.js                      # demo data in demo/
node tools/test-core.js                     # core checks on the demo data
node tools/build.js                         # dist/cgv.html and dist/artifact.html
CGV_L=1000000 node tools/gen-demo.js big    # 1 Mb at 30× for speed tests
```

## Deep data

Targeted panels reach thousands of reads per base. CGV handles them the way IGV does:

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
