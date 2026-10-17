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
| F10 Coloring + MAPQ as a transparency gradient | strand, pair (orientation and insert size), tag; MAPQ works on top of any of them |
| F11 Group and sort | groups: strand, first/second, pair orientation, tag; sort: base, MAPQ, strand, insert size |
| F12 Filters | MAPQ, duplicates, secondary, supplementary, QC fail; downsampling |
| F13 Mate | read details + go to mate; several regions on screen is the next step |
| F14 Read details | all fields and tags, copy |
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
- **Load limit by size:** a window loads on its own when the index points to at most 64 MB
  of compressed alignments; larger windows ask first.

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
