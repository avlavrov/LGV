Write a short article (about 1,200–1,800 words) in English for peerscience.io about LGV, a genome
alignment viewer that runs in a web browser. Style: a concise software / application note. Minimal
discussion and no marketing language; every statement concrete and checkable. Sections: Abstract
(≤150 words), Introduction (3–4 sentences: the need), Implementation, Features, Performance,
Security and privacy, Availability, Limitations and future work, Acknowledgements (AI-assisted
development). Use one table for supported formats and one for performance. Do not compare LGV with
IGV or other viewers by name. Do not invent numbers, users, citations or validation studies: use
only the facts below; mark anything missing as [TODO]. Author and affiliation: [TODO — I will fill
in].

FACTS

What it is
- LGV (Lavrov Genome Viewer): a viewer of read alignments (BAM, CRAM) against a reference genome,
  for medical and human genetics research. Research use only.
- One self-contained HTML file (LGV.html, ~2.6 MB): no installation, no server, no build step for
  the user. Works in Chrome, Edge, Firefox, Safari (Chrome/Edge recommended: folder access to find
  an index next to a file is Chromium-only).
- Files are read in pieces directly from the user's disk or from HTTPS links the user gives
  (CORS + Range requests); nothing is uploaded anywhere. The only connection LGV makes by itself
  is loading the reference chosen to open at start (sequence only).
- Interface in English, Russian as an option. User guide in both languages.
- Built with Claude (Anthropic's AI model) as the coding assistant, in stages, each checked by the
  author on real data (whole genomes ~90 GB, exome and targeted panels).

Implementation
- Own JavaScript code, no runtime dependencies; strict Content Security Policy (no remote scripts,
  styles or fonts). CRAM codecs (rANS 4x8/Nx16, arithmetic, name tokenizer, fqzcomp) adapted from
  htscodecs (BSD 3-Clause); own bzip2 and xz/LZMA2 decoders.
- Reading and decoding run in a pool of Web Workers (one per CPU core, up to 6); reads are drawn
  with WebGL2 instanced rendering; layout and drawing are rebuilt only for rows near the view.
- Coverage is counted from all reads that pass the filters (before downsampling). Above 1 Mb only
  depth is counted (per-bin maxima, 16–64 bp bins); a wide window is split between workers.
  Above 1 Mb depth is first shown at once as an estimate from the BAM index (bytes of reads per
  16 kb from the BAI linear index or CRAI slice sizes, scaled by windows of the same file read in
  detail — the idea of goleft indexcov); up to 10 Mb the exact depth then replaces it.
- Read window (default 200 kb, settable): wider windows show coverage only and offer to load the
  reads with the expected read count and time. Downsampling (default ≤100 reads per 50 bp window)
  is deterministic (read-name hash) and keeps mates together.
- CRAM: every slice is checked against the loaded reference (MD5); BGZF and CRAM blocks are
  CRC-checked; corrupted input stops with a clear message instead of wrong data.

Features
- Reference: GRCh38, GRCh37, T2T-CHM13 presets, or any FASTA (plain, gzip, bgzip; indexed in the
  browser if no .fai). Samples are matched to the reference by contig name (chr prefix and version
  ignored) and length; a sample aligned to another assembly is refused with an explanation.
- Navigation: locus, gene symbol or MANE transcript id; several loci side by side in split panels
  (also "mate in split screen"); chromosome ideogram with GRCh38 bands and their names; contig
  selector; back/forward history.
- Bundled MANE v1.5 (19,373 MANE Select and Plus Clinical transcripts) gene track with exon
  numbers; works offline.
- Reads: mismatches colored, base-quality shading, insertions (click for inserted bases), deletions,
  splice gaps, soft clips shown in the read's color, pointed 3′ ends, pair view, coloring by
  strand / pair orientation and insert size / tag, sorting, grouping, MAPQ transparency with a
  constant outline, filters (MAPQ, duplicates, secondary, supplementary, QC-fail).
- Read card (all fields and tags, mate, base and quality at the clicked position); coverage card
  per position with counts per base and strand and mQV (mean Phred quality of each base).
- Copy read / reverse complement / soft clips / read pair as FASTA.
- Annotation tracks: BED (3–12 columns, narrowPeak/broadPeak), bedGraph, VCF (genotype rows),
  GTF, GFF3 (plain, gzip, bgzip); TDF and bigWig signal tracks read by region, drawn by the UCSC
  track-line rules (graphType, windowingFunction, viewLimits, autoScale, transformFunc, etc.).
  Track-type switches hide all tracks of a kind (BAM off: no reads are read).
- PNG snapshot of the view; a link that reopens the view; URL parameters
  (LGV.html#ref=hg38&file=…&locus=CFTR) to open LGV from a LIMS or report.
- Diagnostic log (file names, loci and timings only; kept in the page).

Performance (local disk, headless Chromium without a GPU; from the project README)
| Data | Window | Reads in window | Load | Build view |
| 30× | 2 kb | 1,339 | 26 ms | 7 ms |
| 30× | 1 Mb | 200,024 | 1,034 ms | 801 ms |
| 3000× panel | 32 bp, first visit | 4,527 | 491 ms | 5 ms |
| 3000× panel | 32 bp, revisit | 4,315 | 64 ms | 1 ms |
| 3000× panel | 20 kb | 784,176 | 2,725 ms | 267 ms |
- Wide views: coverage drawing of a 21 Mb view from per-bin maxima ~3 ms per frame (was 0.6–1.6 s
  base by base); index-based depth estimate appears in <0.3 s for any width.
- [TODO: add timings on the author's machine if wanted.]

Security and privacy
- Every file is treated as untrusted: all counts and offsets are checked; fuzz tests (truncated and
  corrupted BAM, CRAM, BGZF, FASTA, TDF, bigWig, annotation files) must end with a clear message,
  never a hang or wrong data; an audit script checks the build (CSP, no network calls outside the
  user's links, no fixed web addresses except the public reference presets, escaped HTML).
- Nothing is stored except the display profile (in the browser); tokens, links and folder access
  are never stored. No telemetry.

Availability
- Online: https://avlavrov.github.io/LGV/ ; source: https://github.com/avlavrov/LGV (license:
  [TODO]). Single file LGV.html can be saved and used offline (with a local reference).

Limitations and future work
- Browser limits: an index next to a file needs folder permission (Chromium) or choosing both
  files; links need CORS + Range; services without CORS (e.g. UCSC BLAT) cannot be called from the
  page, so BLAT is not included.
- Planned: RNA-seq splice-junction counts and sashimi plots; indexed large annotation files
  (tabix VCF/BED, bigBed); drawings and labels on PNG snapshots; a desktop version (native file
  access) and a company-server version.
