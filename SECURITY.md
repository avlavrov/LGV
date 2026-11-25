# Security in LGV

LGV shows patients' genomic data, so the rule is simple: **data goes only between the source the
user chose and the user's device.** LGV sends nothing anywhere by itself.

## Rules for all code

1. **No third-party packages at run time.** Nothing is installed from npm. Code from others is used
   only as a copy of a fixed, reviewed version inside this repository, with its license, read line by
   line before it is added (see "Code from others" below).
2. **No code from text.** No `eval`, `new Function`, timers with strings, or code loaded at run time.
3. **Text from files is text.** Read names, tags, sample names and file names are shown with
   `textContent` or passed through `esc()`; HTML is built from strings only in a few reviewed places.
4. **Network only to what the user gave.** The only network call is `UrlSource` reading the links the
   user typed, plus the public reference presets: GRCh38 opens at start by default (sequence only,
   from the IGV team's public bucket; View → "At start open" changes it or turns it off). Gene search
   needs no network: the MANE gene table is bundled with LGV, so typed gene names never leave the
   computer. No telemetry, no CDN, no
   fonts from the web. The local build carries a strict Content Security Policy
   (`default-src 'none'`, no `unsafe-eval`, connections only over HTTPS or to this computer).
5. **Storage only for the display profile.** Tokens and links are never saved. Folders the user lets LGV read
   (Chrome, Edge: to find the index next to a chosen file) are read only and kept in memory for the
   session; they are never stored.
6. **Every file is untrusted.** Parsers check every read against the data's bounds, check sizes from
   the file against limits before allocating, verify checksums where the format has them, and turn
   any damage into a clear error instead of a crash or a hang. Heavy reading runs in Web Workers.

## Automatic checks

- `node tools/audit.js` fails on any forbidden pattern above (eval, Function, network calls outside
  `UrlSource`, storage other than the profile, fixed web addresses, HTML built with values not
  passed through `esc()`), checks that the build carries the Content Security Policy, and lists every
  place where HTML is built so a reviewer sees each one.
- `node tools/test-fuzz.js` fuzzes BGZF, the BAM header and records, BAI, FASTA with .fai and .gzi,
  FASTA indexing in the browser, the locus field and a pasted profile (see stage 7 below).
- `node tools/test-cram.js` includes fuzzing: hundreds of randomly damaged CRAM files must each end
  in a clear error, never an unexpected exception or a hang; with checksums switched off the damage
  reaches the record decoder itself.
- Sizes are capped before anything is allocated (`CRAM_LIMITS`): blocks, containers, read length and
  bases per slice, so values that cost no input (a one-symbol Huffman code) cannot inflate memory.
- The audit also checks that the page's scripts declare no top-level name twice.
- CRAM: the CRC32 of every block and container header is always verified (there is no switch to
  turn it off in the page), and each slice's MD5 of its reference stretch is compared with the loaded
  reference, so a different reference cannot silently change the bases shown.
- CRAM and BAM are decoded in Web Workers; the page owns the reference and answers the workers'
  requests for bases, so the workers see only the sample files and the bases they ask for.

- bzip2 and xz/LZMA2 (`src/compress.js`) are LGV's own code: output never grows past the size the
  CRAM block states, distances and table entries are checked, and bzip2 block/stream CRCs and xz
  CRC32/CRC64 are verified. An xz stream written without a check (CRAM's never are: htslib uses
  CRC32) is still covered by the CRAM block's own CRC32.
- `test/codecs`: each codec against its reference encoder (htscodecs' tools in every mode, Python's
  bz2 and lzma), and damaged copies of each.

## Review of all of LGV (stage 7)

Every reader of file data outside CRAM was reviewed line by line and fuzzed the same way as CRAM
(`node tools/test-fuzz.js`, 300 cases per campaign by default, the same seed every run).

- **BGZF** (BAM, bgzip FASTA): every header field is checked (magic, XLEN, BSIZE, ISIZE ≤ 64 KB); the
  CRC32 and length of every inflated block are compared with the block's trailer, as htslib does;
  inflating stops as soon as the output passes the block size, so a tiny "deflate bomb" (51 KB that
  inflates to 50 MB) ends at 64 KB with an error. Before, a damaged block was inflated without any
  check, and damaged data could reach the screen as real reads. Cost: about 20 ms per 1 Mb window.
- **BAM header**: text length, number of contigs, name lengths and contig lengths are checked before
  use (before, a negative length read the header from the wrong place).
- **BAM records**: block_size is capped at 64 MB and must hold the fields it declares; the read name
  must end in NUL; CIGAR operations must be valid; a stored sequence must be as long as the CIGAR
  says (htslib's rule); read length ≤ 8 Mb. Before, a damaged record made LGV buffer the rest of the
  file looking for its end, or loop over a 268-million-base CIGAR operation: the fuzzer ran the old
  code out of memory (checked by switching the record check off).
- **Coverage** visits only the part of an alignment inside the window, so a huge deletion or match
  costs nothing outside it (BAM and CRAM).
- **Tags**: a tag area that ends early or a `B` array longer than the record ends the tag list
  instead of throwing; what was read before is kept.
- **BAI**: every count (contigs, bins, chunks, intervals) is checked against the bytes left before it
  is used (before, a negative count walked backwards through the index).
- **FASTA**: `.fai` lines must be whole numbers with bytes per line = bases + 1 or 2, else an error
  names the line; `.gzi` entries must increase and blocks hold at most 64 KB; an ordinary gzip FASTA
  is unpacked to at most 4 GB.
- **Links**: a server that ignores Range, or sends more than was asked for, cannot make LGV
  download more than it needs: bodies are read in pieces and cut at the requested length; a file
  read whole (an index) is capped at 512 MB; reading past 256 MB of a server without Range support is
  an error that says so.
- **Workers**: a worker that does not answer in 60 s (30 s for opening) is stopped and replaced, the
  user sees an error, and the track opens its file again in the new worker on the next read. The page
  does not retry the same data itself, so a stuck decoder never freezes the view.
- **Locus field**: positions are clamped to the contig, `500-100` is read as `100-500`, and numbers
  too large to be exact are refused.
- **Pasted or stored profile**: besides own keys only, every value must have the default's type, numbers
  are clamped to a range (read height 1–24 px, and so on), choices must be one of the known ones and
  tags two characters; unknown keys are dropped.
- **Audit additions**: file readers may throw only `CgvError`; every reader must have a fuzzing
  campaign in the tests (the audit fails if one is removed).
- Found while fixing: a record without a stored sequence (SEQ `*`, common for secondary alignments)
  was drawn with every base as a mismatch; it now shows no mismatches.

## Data from others

| Data | Source | Version | Terms | Use |
| --- | --- | --- | --- | --- |
| Chromosome bands, GRCh38 (`data/cytoband.grch38.tsv`) | ideogram by Eric Weitz, `data/bands/native/homo-sapiens.json` | 850-band resolution, NCBI | CC0 (public domain) | chromosome map |
| MANE Select and MANE Plus Clinical transcripts, GRCh38 (`data/mane.grch38.tsv`, made by `tools/make-mane.py`) | MANE summary (NCBI / EMBL-EBI) and Ensembl 111 exon coordinates, as exported by genome-nexus-importer (`data/common_input/MANE.GRCh38.v1.2.summary.txt`, `data/grch38_ensembl111/export/ensembl_biomart_transcripts.json.gz`) | MANE v1.2 | public data, no restrictions (NCBI, Ensembl) | gene track and gene search; all 19,268 transcripts on chr1–22, X, Y; every Ensembl transcript version matched the MANE summary |

## Opening from a link

`lgv.html#ref=…&file=…&index=…&locus=…` opens data named in the address. Only https links (and http on
this computer) are read, at most 50 files; access tokens are never taken from or put into a link;
LGV still sends nothing anywhere — the links are only read. The "Link" button builds such a link from
what was opened by link; files from the computer are left out.

## Annotation tracks (BED, bedGraph, VCF, GTF, GFF3)

`src/tracks.js`, own code. A file is read whole, capped at 512 MB unpacked (gzip output stops at the
cap; bgzip blocks are checked like BAM's) and 5 million features. Every coordinate must be a plain
whole number with end ≥ start; BED blocks must be in order and inside the feature; text fields are
cut to a fixed length. A line that does not parse is skipped and counted. Text from the file is
drawn on the canvas or shown through esc(). Fuzzed in tools/test-fuzz.js.

## Signal tracks (TDF, bigWig)

`src/signal.js` (reading) and `src/signalview.js` (drawing), own code; the layouts were taken from
IGV's TDF reader and writer and from UCSC's bbiFile.h, the drawing rules from UCSC's wigTrack.c, no
code copied. The file is read by region, never whole. Every count, offset and size is checked against
the bytes actually read before it is used; a block is capped at 64 MB before and after
decompression; tree walks are capped in depth (24) and nodes (100,000); contig and track counts are
capped. One view reads at most 16 MB of full-resolution data (4 MB over a link) and at most 3,000 TDF
tiles; beyond that the file's summaries are used. Reads run a few at a time; requests the view no
longer needs are dropped. A track line is parsed by a fixed list of known keys, each value checked;
everything else is ignored. Rows are painted into an image buffer, so drawing time does not grow with
the number of values. The drawing choices a user makes are kept in `localStorage` per file name
(`lgv.signal.v1`, at most 100 files). `tools/test-signal.js` checks the readers against files made by
igvtools, UCSC tools and bigtools and fuzzes both.

## Code from others

| Code | From | Version | License | Status |
| --- | --- | --- | --- | --- |
| CRAM block decoders (`src/htscodecs.js`): rANS 4x8, rANS 4x16, arithmetic, FQZcomp, name tokeniser | htscodecs `javascript/` rans.js, rans4x16.js, arith_gen.js, arith_sh.js, byte_model.js, fqzcomp.js, tok3.js, iostream.js (Genome Research Ltd, James Bonfield) | commit 0a487b5 | BSD 3-clause, notice kept in the file | decoding only; read line by line and changed as below |
| Test data only: SAM test cases | htslib `test/` | commit 446b29a | MIT/Expat (test/cram/htslib/LICENSE-htslib) | used by tools/test-cram.js, not shipped |

### Changes made to htscodecs code

- The byte reader checks every read against the end of the data and throws a clear error. The
  original returned `undefined` past the end; the frequency tables are read in loops that stop only
  on a 0 byte, so a cut-off block made the original loop forever (checked: the original
  `rans.js` still runs after 5 s on a 10-byte block; LGV stops at once with an error).
- Frequencies are checked: each at most 4096, their sum at most 4096; codes above the sum and
  contexts without frequencies are errors instead of writing symbol 256 or reading undefined tables.
- A table running past symbol 255, an unknown rANS order and sizes in the block header that differ
  from the CRAM block's are errors.
- rANS 4x16 and arithmetic: run lengths, packed data, nested stripes and output sizes are checked;
  a packed block with more than 16 symbols is an error (the original returned uninitialised memory);
  a range that reaches zero is an error (the renormalising loop would never end); arithmetic models
  of order 1 and FQZcomp's 65,536 quality models are made on first use instead of all up front.
- FQZcomp: the parameter selector count read `gflags.nparam` (always undefined) and a duplicated
  record decoded its flag twice and copied with an undefined length; both now follow htscodecs' C
  decoder. Tables that never end, record lengths outside the block and missing selectors are errors.
- Name tokeniser: references to missing token streams or to names that are not there are errors
  (the original crashed with a TypeError); the name count is checked against the block size.
- Node's `Buffer` replaced by `Uint8Array`; 32-bit values read unsigned.
- Encoders, command-line tools and the `bzip2` package are not taken.

## Findings

| Date | Finding | Fixed |
| --- | --- | --- |
| 2026-10-21 | htscodecs' JavaScript rANS decoder loops forever on a cut-off block (no bounds checks) | LGV's copy checks every read; reported here, upstream report to consider |
| 2026-10-21 | htscodecs' JavaScript decoders on cut-off blocks: rANS 4x16 and arithmetic return a full-size block of wrong data without an error, FQZcomp loops forever, the name tokeniser crashes | LGV's copies stop with a clear error; tested on blocks cut to 12, 40 and 100 bytes |
| 2026-10-21 | A pasted display profile with a `__proto__` key could replace the profile object's prototype (no effect outside it, but not allowed) | only the profile's own keys are taken |
| 2026-10-28 | BAM, BAI, BGZF and FASTA readers trusted sizes and counts from the file: a damaged record could exhaust memory, a damaged BGZF block could show wrong reads, a negative count misread the index | all checked (stage 7); fuzzed in tools/test-fuzz.js |
| 2026-10-28 | Inflating was unbounded: a small damaged or hostile block could expand to gigabytes | output capped at the declared size (BGZF, CRAM gzip, .crai, gzip FASTA) |
| 2026-10-28 | A stuck worker would never answer and the track would wait forever | worker timeout, replaced worker, clear error |
