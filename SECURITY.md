# Security in CGV

CGV shows patients' genomic data, so the rule is simple: **data goes only between the source the
user chose and the user's device.** CGV sends nothing anywhere by itself.

## Rules for all code

1. **No third-party packages at run time.** Nothing is installed from npm. Code from others is used
   only as a copy of a fixed, reviewed version inside this repository, with its license, read line by
   line before it is added (see "Code from others" below).
2. **No code from text.** No `eval`, `new Function`, timers with strings, or code loaded at run time.
3. **Text from files is text.** Read names, tags, sample names and file names are shown with
   `textContent` or passed through `esc()`; HTML is built from strings only in a few reviewed places.
4. **Network only to what the user gave.** The only network call is `UrlSource` reading the links the
   user typed (plus the public reference presets when the user picks one). No telemetry, no CDN, no
   fonts from the web. The local build carries a strict Content Security Policy
   (`default-src 'none'`, no `unsafe-eval`, connections only over HTTPS or to this computer).
5. **Storage only for the display profile.** Tokens and links are never saved.
6. **Every file is untrusted.** Parsers check every read against the data's bounds, check sizes from
   the file against limits before allocating, verify checksums where the format has them, and turn
   any damage into a clear error instead of a crash or a hang. Heavy reading runs in Web Workers.

## Automatic checks

- `node tools/audit.js` fails on any forbidden pattern above (eval, Function, network calls outside
  `UrlSource`, storage other than the profile, fixed web addresses, HTML built with values not
  passed through `esc()`), checks that the build carries the Content Security Policy, and lists every
  place where HTML is built so a reviewer sees each one.
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

- bzip2 and xz/LZMA2 (`src/compress.js`) are CGV's own code: output never grows past the size the
  CRAM block states, distances and table entries are checked, and bzip2 block/stream CRCs and xz
  CRC32/CRC64 are verified. An xz stream written without a check (CRAM's never are: htslib uses
  CRC32) is still covered by the CRAM block's own CRC32.
- `test/codecs`: each codec against its reference encoder (htscodecs' tools in every mode, Python's
  bz2 and lzma), and damaged copies of each.

## Code from others

| Code | From | Version | License | Status |
| --- | --- | --- | --- | --- |
| CRAM block decoders (`src/htscodecs.js`): rANS 4x8, rANS 4x16, arithmetic, FQZcomp, name tokeniser | htscodecs `javascript/` rans.js, rans4x16.js, arith_gen.js, arith_sh.js, byte_model.js, fqzcomp.js, tok3.js, iostream.js (Genome Research Ltd, James Bonfield) | commit 0a487b5 | BSD 3-clause, notice kept in the file | decoding only; read line by line and changed as below |
| Test data only: SAM test cases | htslib `test/` | commit 446b29a | MIT/Expat (test/cram/htslib/LICENSE-htslib) | used by tools/test-cram.js, not shipped |

### Changes made to htscodecs code

- The byte reader checks every read against the end of the data and throws a clear error. The
  original returned `undefined` past the end; the frequency tables are read in loops that stop only
  on a 0 byte, so a cut-off block made the original loop forever (checked: the original
  `rans.js` still runs after 5 s on a 10-byte block; CGV stops at once with an error).
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
| 2026-10-21 | htscodecs' JavaScript rANS decoder loops forever on a cut-off block (no bounds checks) | CGV's copy checks every read; reported here, upstream report to consider |
| 2026-10-21 | htscodecs' JavaScript decoders on cut-off blocks: rANS 4x16 and arithmetic return a full-size block of wrong data without an error, FQZcomp loops forever, the name tokeniser crashes | CGV's copies stop with a clear error; tested on blocks cut to 12, 40 and 100 bytes |
| 2026-10-21 | A pasted display profile with a `__proto__` key could replace the profile object's prototype (no effect outside it, but not allowed) | only the profile's own keys are taken |
