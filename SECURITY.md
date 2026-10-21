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

## Code from others

| Code | From | Version | License | Status |
| --- | --- | --- | --- | --- |
| rANS 4x8 decoder (`src/htscodecs.js`) | htscodecs `javascript/rans.js`, `iostream.js` (Genome Research Ltd, James Bonfield) | commit 0a487b5 | BSD 3-clause, notice kept in the file | decoding only; read line by line and changed as below |
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
- Node's `Buffer` replaced by `Uint8Array`; 32-bit values read unsigned.
- Encoders, command-line tools and the `bzip2` package are not taken.

## Findings

| Date | Finding | Fixed |
| --- | --- | --- |
| 2026-10-21 | htscodecs' JavaScript rANS decoder loops forever on a cut-off block (no bounds checks) | CGV's copy checks every read; reported here, upstream report to consider |
| 2026-10-21 | A pasted display profile with a `__proto__` key could replace the profile object's prototype (no effect outside it, but not allowed) | only the profile's own keys are taken |
