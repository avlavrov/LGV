# LGV development to-do

Ideas and agreed work for later. Items move to the user guide and README when done.

## Agreed, not started

- **Read card in its own window** (for a second monitor): Chrome and Edge can open an always-on-top
  window with page content (Document Picture-in-Picture); the card would be rendered there and
  updated on every click. About 100 lines; Chrome / Edge only, local `lgv.html` only (claude.ai
  will likely not allow it). A plain pop-up window was rejected: two windows talking to each other
  goes against the security rule "the page does not talk to other windows". (2026-11-01: postponed)
- **IGV read-menu features**: the user marks which to add in the Google Sheet
  "LGV — меню ридов IGV: что добавить".

## Offered, waiting for a decision

- Share link: a link that opens the same loci (and public data links).
- Check how IGV desktop (htsjdk) handles damaged CRAM files (igv.js does not check CRC by default).
- BLAT (removed 2026-11-26: this version sends nothing anywhere). From a web page UCSC's hgBlat
  cannot be read (no CORS); igv.js (IGV-Web) sends BLAT to its own proxy, igv.org/services/
  blatUCSC.php, which forwards to UCSC (companies can point it at their own); JBrowse 2 uses its own
  proxy holding a UCSC API key. Options if it comes back: an own proxy (a company's server or a small
  cloud function), the user's UCSC API key, or the desktop version (no CORS there). Each sends read
  bases outside the computer: needs the user's decision and a Privacy note.
- Wide windows on a slow link: read depth-only windows from the index's own statistics or a
  coverage file (bigWig / TDF next to the BAM) instead of the reads.
- Report to htscodecs upstream: its JavaScript decoders hang or return wrong data on cut-off
  blocks (details in SECURITY.md, Findings).

## Annotation tracks, next steps

- Large indexed files read in pieces: VCF / BED with tabix `.tbi` / `.csi`, bigBed (now: whole file in
  memory, up to 512 MB). TDF and bigWig are read in pieces (done).
- Signal tracks: a TDF whole-genome view (its `/All` datasets); bigWig settings from a track hub or a
  track line given with the link; overlaying several rows in one (UCSC multiWig); a scale shared by
  the rows of a file.
- VCF: more than 8 genotype rows (scrolling), coloring by allele frequency.
- Tracks in the display profile / share link.

## Product stages ahead

- Desktop version (native file access: the index next to a file without asking).
- Company-server version (data stays on the server; viewing through the browser).
