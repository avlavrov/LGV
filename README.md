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
- **Reference**: add a FASTA at any time, separately from the BAM. One `.fa` / `.fasta` file is enough: without a `.fai` CGV indexes it in the browser in one pass (plain, bgzip or gzip). With an index it opens instantly — `.fa` + `.fai`, or
  bgzip-compressed `.fa.gz` + `.fai` + `.gzi`, from disk or by link. Public GRCh38, GRCh37
  and T2T-CHM13 v2.0 copies (hosted by the IGV team on Amazon S3) are offered as presets.
  Without a reference, CGV rebuilds reference bases from MD tags where reads cover them.

The interface is in English; Russian is available under **View → Screen → Language**.
More languages are one dictionary in `src/i18n.js`.

## Stage 0 coverage of the MVP spec

| Requirement | State |
| --- | --- |
| F1 BAM by HTTPS link (Range) | done; CRAM is the next step |
| F2 Local BAM + BAI | done, including drag and drop |
| F3 Reference (GRCh38, GRCh37, T2T-CHM13) | done: files or link, plain or bgzip, public presets |
| F4 Access token / presigned URL | done |
| F5 Go to coordinate | done; gene names come with the gene track |
| F6 No hard window limit | windows up to 3 Mb load at once, larger on request; isolated sequencing errors fade when zoomed out |
| F8 Coverage with allele highlight | done, threshold adjustable |
| F9 Smooth read height | done, 1–24 px |
| F10 Coloring + MAPQ as a transparency gradient | strand, pair (orientation and insert size), tag; MAPQ works on top of any of them |
| F11 Group and sort | groups: strand, first/second, pair orientation, tag; sort: base, MAPQ, strand, insert size |
| F12 Filters | MAPQ, duplicates, secondary, supplementary, QC fail |
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

## Measurements (synthetic 30×, local disk, headless Chromium without a GPU)

| Region | Reads | Read | Build view |
| --- | --- | --- | --- |
| 2 kb | 1,339 | 26 ms | 7 ms |
| 50 kb | 30,063 | 175 ms | 143 ms |
| 1 Mb | 200,024 | 1,034 ms | 801 ms |

Over a network, the server's latency adds to the read time.

## License

Not chosen yet.
