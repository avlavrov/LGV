# CGV user guide

CGV (Convenient, custom, genomic viewer) shows read alignments (BAM, CRAM) against a reference
genome, for medical and human genetics. It runs in a browser from one file, `cgv.html`, with no
installation and no server: files are read in pieces straight from your disk or from links you
give, and nothing is uploaded anywhere.

This guide describes CGV as it is now and is updated with every change. Russian version:
[USER_GUIDE.ru.md](USER_GUIDE.ru.md).

## Contents

1. [Starting](#1-starting)
2. [The screen](#2-the-screen)
3. [Opening data](#3-opening-data)
4. [Moving around](#4-moving-around)
5. [Several regions side by side](#5-several-regions-side-by-side)
6. [Chromosome map](#6-chromosome-map)
7. [Genes (MANE)](#7-genes-mane)
7a. [Your own tracks: BED, bedGraph, VCF, GTF, GFF3](#7a-your-own-tracks-bed-bedgraph-vcf-gtf-gff3)
8. [Coverage](#8-coverage)
9. [Reads: how they are drawn](#9-reads-how-they-are-drawn)
10. [Read details](#10-read-details)
11. [Right-click menu](#11-right-click-menu)
12. [Several samples](#12-several-samples)
13. [View settings](#13-view-settings)
14. [Keys and mouse](#14-keys-and-mouse)
15. [Messages and limits](#15-messages-and-limits)
16. [Privacy and security](#16-privacy-and-security)

## 1. Starting

Open `cgv.html` in Chrome, Edge, Firefox or Safari (Chrome or Edge recommended: some features need
them, see 3.1). At start CGV loads the **GRCh38 / hg38** reference (sequence only, from the IGV
team's public storage) and opens on BRCA1. A short message confirms the reference is loaded.

- What opens at start is set in **View → At start open**: GRCh38, GRCh37, T2T-CHM13, the demo data,
  or nothing.
- Without an internet connection the reference cannot load; CGV then opens the built-in demo data
  (a synthetic contig `demo` with example variants) and says so.
- The interface is in English; Russian is under **View → Language**.

## 2. The screen

From top to bottom:

- **Top bar**: the locus field with **Go**, **← →** (back / forward), **− +** (zoom), the width of the
  view, the loaded reference and samples, the status line (reads shown, load and draw times),
  **Open** and **View**. **H** hides the bar; a button in the corner brings it back.
- **Chromosome map**: the chromosome on screen, with a red box where the view is (section 6).
- **Panel headers**: only with several regions side by side (section 5).
- **Ruler** and **reference sequence**: coordinates are 1-based. The reference shows as
  colored blocks, and as letters when zoomed in far enough.
- **Gene track**: MANE transcripts (section 7).
- **Samples**: each with its coverage and its reads. With several samples, a small label at the top
  of each sample's reads names it (file and sample; it may cover a row or two), with × at the right
  to close it.

The interface is kept compact on purpose: as much of the screen as possible goes to data.
Explanations are kept out of the way, as **tooltips**: rest the pointer on a button, a setting, a
line of a card, the chromosome map, the ruler, the gene track or the coverage, and a short
explanation appears after a moment (the browser's own tooltip — it does not pop up while the
pointer moves, and there are none over the reads).

## 3. Opening data

Press **Open**.

### 3.1 From disk

- **Choose files**: one or more `.bam` / `.cram` files. Each is paired with its index by name:
  `sample.bam.bai` or `sample.bai`, `sample.cram.crai` or `sample.crai` (case ignored). A reference
  (`.fa` / `.fasta` / `.fna`, optionally `.gz`, with `.fai` and `.gzi`) can be chosen at the same
  time.
- **Just the .bam or .cram, without its index** (Chrome, Edge): a web page may read only the files
  you choose, so CGV asks once — **Allow the folder**. The folder picker opens where the file lies;
  confirm it and CGV takes the index from there itself. Later files from the same folder open
  without asking, until the page is closed. The access is read-only and never stored. In Firefox
  and Safari this is not possible: choose the file together with its index, or the whole folder.
- **Choose a folder**: every sample in it is found and paired with its index. A folder with several
  samples first shows a list to pick from (the first four are ticked).
- **Drag and drop** files or a folder onto the screen works the same way.

Opening the same file twice does nothing: it is already open.

### 3.2 By link

Paste one or more HTTPS links to BAM or CRAM files, one per line. The index is looked for next to
each file (`.bam.bai`, then `.bai`; `.cram.crai`, then `.crai`) unless you give it. The server must
allow CORS and Range requests (the usual case for cloud buckets with signed links). An access token
can be given for servers that need one; it is not saved.

### 3.3 Reference

The reference sets the genome. Loading another reference switches the genome and closes
the open samples.

- **Presets**: GRCh38, GRCh37 and T2T-CHM13 v2.0 (public copies hosted by the IGV team).
- **Your own**: a FASTA file from disk or by link. Without a `.fai`, CGV indexes the file in the
  browser in one pass (plain, bgzip or gzip; a progress message shows). With `.fai` (and `.gzi` for
  bgzip) it opens at once.
- Samples are matched to the reference by contig name (case, `chr` prefix and version suffix
  ignored), then by a length only one contig on each side has. A sample made against another build
  (a large contig with the same name and a different length, e.g. chr1 of GRCh37 against GRCh38) is
  refused with an explanation, because its coordinates would be wrong.
- Without a reference, the first sample's header gives the coordinates, and reference bases are
  rebuilt from MD tags where reads cover them.
- **CRAM** stores reads as differences from the reference it was made with. CGV checks each part of
  the file against the loaded reference (MD5) and refuses to show reads built on a different
  reference, rather than show wrong bases.

## 4. Moving around

- **Locus field**: type `chr7:117,559,500-117,560,000`, `chr7:117559590` (a 100 bp window around
  the position), a contig name (`chrX`, `X` and `x` are the same), a **gene symbol** (`CFTR`), or a
  MANE transcript (`NM_000492`, `ENST00000003084`, with or without version). Enter or **Go**.
  Suggestions appear while typing a gene. `500-100` is read as `100-500`; positions beyond the
  contig are moved to its end.
- **Several loci or genes separated by spaces** open side by side (section 5).
- **Drag** the reads left or right to move along the genome.
- **Zoom**: **− +** buttons and keys (all panels at once), Ctrl + wheel or a trackpad pinch, double
  click (zoom in at that point).
- **Drag on the ruler** to zoom to the selected stretch; a click on the ruler centers the view
  there.
- **← →** buttons: back and forward between the places you went to — typed loci and genes,
  examples, "go to mate", "mate in split screen", "show whole gene", clicks on the chromosome map.
  Panning and zooming are not separate steps. Back after several regions side by side brings all of
  them back.

## 5. Several regions side by side

The view splits into panels:

- **Right click a read → View mate region in split screen** (or the button in the read details)
  opens a panel next to the current one, centered on the mate, at the same scale. The read and its
  mate are outlined in every panel, and the new panel scrolls so the mate is in sight.
- **Several loci or genes in the locus field** (separated by spaces) open one panel each; one locus
  turns the view back into one panel.
- Each panel has its own ruler, zoom, loaded window, row layout and scroll. The header shows its
  number and locus; the active panel (the last one clicked) is underlined, and the keys act on it.
- **×** in a panel's header closes it.
- At least 6 panels are always allowed, more as the width allows (about 60 px each).

## 6. Chromosome map

A strip across the top shows the chromosome on screen: for GRCh38 with its bands (dark = G-positive
bands, red = centromere); other genomes get a plain bar. A **red box** marks where the view is; the
active panel's box is filled. With several panels the boxes are numbered like the panel headers,
and panels on different chromosomes get their chromosomes side by side. **Click or drag** on the map
to move the active panel there.

## 7. Genes (MANE)

Under the reference sits the gene track: MANE Select and MANE Plus Clinical transcripts on GRCh38
(19,268 transcripts, MANE v1.2), built into CGV — it works offline, and gene names you type never
leave the computer.

- Thick boxes are coding exons, thin boxes UTRs, the line with arrows is the intron and the strand.
  Plus Clinical transcripts have their own color. Exon numbers (`ex10`) appear inside exons when
  zoomed in, counted in the transcript's direction.
- **Click a transcript**: gene and HGNC id, full name, MANE status, Ensembl and RefSeq ids, gene
  span (+/− strand), number of exons, coding part, and where the clicked position falls: `exon 11 of
  27`, `(UTR)`, or `intron 3 of 26`. **Show whole gene** zooms to it; **Copy** copies the text.
- The track shows only when the coordinates are GRCh38; it can be turned off in **View**.

## 7a. Your own tracks: BED, bedGraph, VCF, GTF, GFF3

Annotation files open the same way as samples — **Open → From disk** (choose, drop, or together with
samples), or **By link** — and are recognized by their name: `.bed` (also `.narrowPeak`,
`.broadPeak`), `.bedgraph` / `.bg`, `.vcf`, `.gtf`, `.gff` / `.gff3`, each also as `.gz` (gzip or
bgzip). No index is needed: the file is read whole (up to 512 MB unpacked).

They appear as bands under the gene track, in the order opened, fixed while the samples scroll.
Each band starts with a thin line holding its file name, and × at the right closes it. Contig
names match with or without `chr`.

- **BED** (3 to 12 columns): features in two rows with their names; blocks are drawn as exons, the
  thick part as coding, a strand gives arrows (or a pointed end), `itemRgb` colors a feature.
- **GTF / GFF3**: exons and CDS are grouped into transcripts (by `transcript_id` or `Parent`) and drawn
  like the gene track; lines without exons (repeats, regions) are drawn as they are.
- **bedGraph**: a bar chart; the scale (largest visible value, `±` when there are negative values) is
  on the name line.
- **VCF**: a mark per variant (red; gray when FILTER is not PASS), and under it a row per sample,
  up to 8 samples, with the sample name: light gray = homozygous reference, blue = heterozygous,
  cyan = homozygous alternative, outline only = no call. With more samples the name line says how
  many there are.
- **Click a feature** for its card: position and strand in the top row with Copy and ×; then the
  track name and the details — for BED name, score, blocks and thick part; for GTF / GFF3 type and
  all attributes; for bedGraph the value; for VCF ID, REF, ALT, QUAL, FILTER, every INFO field and
  each sample's FORMAT fields (up to 30 samples).
- Lines that cannot be read are skipped and counted: a message says how many.
- Tracks stay open when another reference is loaded (they are positions, not reads).

## 8. Coverage

Each sample has a coverage track: gray bars of read depth; the scale is in the corner (`0–31×`).
Where non-reference bases make up at least 20% of reads (**View → Highlight alleles from**), the
bar shows colored parts by base.

- Coverage counts **every** read that passes the filters, also where the reads shown are
  downsampled.
- Above 1 Mb wide, coverage is counted as depth only (no allele colors): it loads several times
  faster, and alleles cannot be seen at that scale anyway.
- **Click the coverage** for the counts at that base. The top row shows the position,
  **Copy** and **×**; then:

  ```
  Total count = 22
  A = 0 (0%, 0+, 0−)
  C (reference) = 12 (55%, 4+, 8−) mQV 37
  G = 10 (45%, 8+, 2−) mQV 34
  T = 0 (0%, 0+, 0−)
  DEL = 0 (0+, 0−)
  INS = 0 (0+, 0−)
  ```

  Each base: number of reads, share of the total, reads on the forward (+) and reverse (−) strand,
  and **mQV**. **Total count** is A + C + G + T + N (deletions
  are counted separately). **DEL**: reads with a deletion at this position. **INS**: reads with an
  insertion right after this base. The counts come from **all reads at that position**, read again
  from the file, so downsampling does not change them; the filters from **View** (duplicates, MAPQ,
  secondary, QC-failed) do apply, as in the view. The counted column is outlined.

  **mQV is not a common standard measure**, so here is exactly what it is: the mean of the
  base qualities (Phred QV, as stored in the file) of that base at this position, over the reads
  counted for it: `mQV = (QV₁ + … + QVₙ) / n`, rounded. Reads without stored qualities are left out
  of the mean (if none has them, mQV is not shown). It helps to tell a real allele, usually read with
  qualities like the reference base, from a cluster of sequencing errors with low qualities. It is
  an arithmetic mean of Phred values, not of error probabilities, so a few very low-quality bases
  pull it down less than they would an error-rate average. Hovering a line of the card shows these
  explanations briefly.

## 9. Reads: how they are drawn

- **Body**: gray by default, or colored as chosen (right click → Color by): strand; pair
  orientation and insert size (red = insert larger than expected, blue = smaller, green = RF,
  teal = FF or RR, orange = mate unmapped, a color per chromosome for mates on another chromosome;
  expected sizes are the 0.5–99.5% range of the loaded pairs); a tag value (for example
  `HP` for haplotypes).
- **Strand**: the 3′ end of each read is pointed, so the strand shows under any coloring.
- **MAPQ**: reads with low mapping quality are more transparent (from MAPQ 0 at 18% to full at
  MAPQ 30 and above; adjustable).
- **Mismatches**: only bases that differ from the reference are colored: A green, C blue, G orange,
  T red; N is not marked. Letters appear inside reads from 7 pixels per base when rows are at least
  9 px tall.
- **Base quality**: below Q5 a mismatch is drawn at 20% strength, rising to full at Q20 (can be
  turned off: View → Fade low-quality bases).
- **Zoomed out**, lone mismatches fade so that sequencing errors do not clutter the picture;
  positions where a non-reference base passes the allele threshold always stay bright.
- **Insertions**: a purple bar, with its length written above from 10 pixels per base.
  **Deletions**: a dark line. **Spliced gaps (N)**: a thin gray line.
- **Soft clips** (shown by default; right click → Show soft-clipped bases): the clipped bases
  continue the read in its own color, slightly thinner so the clip boundary shows; clipped bases
  that differ from the reference are colored like mismatches. With clips hidden, a pink tick marks
  where a read is clipped.
- **Pairs** (right click → View as pairs): both mates of a pair in one row, joined by a thin line;
  where mates overlap the overlap is shaded. A mate outside the loaded window is shown by the line
  to where it lies.
- **Downsampling**: in deep regions at most 100 reads are kept per 50 bp window (adjustable); marks above the reads show where reads were left out. Coverage still counts every read.
- **Rows** never reshuffle on their own: panning, reloading or recoloring keeps every read in its
  row. A new layout is made only when grouping, sorting, filters or downsampling change.

## 10. Read details

**Click a read** for its details: name, sample, read group, length,
flags; mapping (primary, secondary or supplementary, MAPQ), reference span and strand, CIGAR,
clipping; mate (mapped, start, insert size, first or second, pair orientation); all tags (the RG tag
is used for the sample and listed as hidden). At the bottom: the clicked **Location** and **Base**
with its quality, e.g. `Base = C @ QV 37 (soft-clipped)`.

- **Read / Read and mate** at the top switches between one read and both mates side by side (the
  mate is shown when it is loaded in any panel). With both shown, each mate that covers the clicked
  position has its own Location / Base — mates that overlap show only one letter on screen.
- The base the card is about is **outlined** with a thin line in the read (and in the mate).
- **Go to mate** moves the view to the mate; **Mate in split screen** opens it in a new panel;
  **Copy** copies the text.
- **×** at the top right closes the card (also Escape).
- **Drag the card by its top row** to move it; later cards open where you left it (until the page
  is closed), even if part of it is off screen.

## 11. Right-click menu

Right click on the reads:

- **View mate region in split screen** (on a paired read with a mapped mate).
- **Copy read sequence**, **Copy left / right soft-clipped sequence**, **Copy read pair (FASTA)**.
  All are copied as FASTA with the read name: `>name/1`; clips as `>name/1_left_clip` and
  `>name/1_right_clip`; a pair as `>name/1` and `>name/2`. The sequence is as stored in the file, on
  the reference's forward strand. Items that do not apply are grayed out.
- **View as pairs**, **Show soft-clipped bases**.
- **Color by**: none, strand, pair orientation and insert size, tag value (with the tag field).
- **Sort by** at the clicked position: start (unsorted), base, MAPQ, strand, insert size. Reads
  covering that position move to the top in that order. **S** sorts by base at the screen center.

## 12. Several samples

Every sample opened becomes a track under the same ruler, at the same position. Each has its own
coverage, reads and scroll. When they do not fit, each keeps at least its minimum height
(**View → Minimum sample height**, 240 px by default) and the whole stack scrolls under the fixed
ruler and gene track:

- the wheel over a sample's reads scrolls its rows; past their end it scrolls the whole stack;
- the wheel or a vertical drag over coverage scrolls the stack;
- the scrollbar on the right; PageUp / PageDown.

## 13. View settings

**View** opens the settings; they are kept in the browser as a **display profile**, which can be
copied as JSON, pasted from someone else and applied, or reset.

- **Reads**: read height (1–24 px), color, tag for color, fade low-quality bases, show soft clips,
  MAPQ transparency (from / to MQ, opacity at low MQ).
- **Group and sort**: group by strand, first / second of pair, pair orientation or a tag; sort rows.
- **Filters**: minimum MAPQ; hide duplicates, secondary, supplementary, QC-failed reads;
  downsampling (reads per window, window size). Filters and downsampling act while reading, so
  changing them reads the window again.
- **Coverage**: minimum sample height, coverage track height, allele highlight threshold.
- **Screen**: language, hide the top bar, gene track on / off, what opens at start.

## 14. Keys and mouse

| Action | How |
| --- | --- |
| Move along the genome | drag the reads; ← →; Shift + wheel |
| Zoom (all panels) | + −, buttons − +; Ctrl + wheel / pinch; double click (in, at that point) |
| Zoom to a region | drag on the ruler |
| Center on a position | click the ruler |
| Scroll reads | wheel over reads; drag up / down; ↑ ↓ |
| Scroll the stack of samples | wheel over coverage; scrollbar; PageUp / PageDown |
| Back / forward | ← → buttons in the top bar |
| Read details | click a read |
| Counts at a base | click the coverage |
| Transcript details | click the gene track |
| Menu | right click on the reads |
| Sort by base at the center | S |
| Hide / show the top bar | H |
| Close a card, a menu or a message | Escape |

## 15. Messages and limits

- Windows that would read more than 64 MB of compressed data ask first, with the size; **load
  anyway** reads it.
- A damaged file gives a clear message saying what is wrong (for example a checksum mismatch in a
  block) instead of showing wrong data.
- Reading that takes longer than 60 s is stopped with a message; the file may be damaged.
- Messages at the bottom close with × or a click; short ones go by themselves.
- **No index found**: choose the index together with the file, the whole folder, or allow the folder
  when asked (3.1).
- **Another genome build**: the sample was aligned to a different reference than the one loaded
  (3.3).

## 16. Privacy and security

- Data goes only between the files or links you chose and your computer. CGV sends nothing
  anywhere; the only connection it makes by itself is loading the reference chosen to open at
  start (sequence only).
- Gene search and the gene track are built in: typed gene names are not sent anywhere.
- Only the display profile is saved in the browser. Tokens, links and folder access are never
  saved.
- Every file is treated as untrusted: damaged or hostile files end in a clear message, never in
  wrong data, a hang or a crash. Details are in [SECURITY.md](../SECURITY.md).
