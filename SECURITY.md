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
  in a clear error, never an unexpected exception or a hang.

## Code from others

| Code | From | Version | License | Status |
| --- | --- | --- | --- | --- |
| — | — | — | — | none yet; htscodecs (CRAM codecs, Genome Research Ltd, BSD) is planned, copied at a fixed commit and reviewed |

## Findings

| Date | Finding | Fixed |
| --- | --- | --- |
| 2026-10-21 | A pasted display profile with a `__proto__` key could replace the profile object's prototype (no effect outside it, but not allowed) | only the profile's own keys are taken |
