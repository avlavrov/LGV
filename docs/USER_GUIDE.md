# LGV user guide

LGV (Lavrov Genome Viewer) shows read alignments (BAM, CRAM) against a reference
genome, for medical and human genetics. It runs in a browser from one file, `lgv.html`, with no
installation and no server: files are read in pieces straight from your disk or from links you
give, and nothing is uploaded anywhere.

This guide describes LGV as it is now and is updated with every change. Russian version:
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
7b. [Signal tracks: TDF and bigWig](#7b-signal-tracks-tdf-and-bigwig)
8. [Coverage](#8-coverage)
9. [Reads: how they are drawn](#9-reads-how-they-are-drawn)
10. [Read details](#10-read-details)
11. [Right-click menu](#11-right-click-menu)
12. [Several samples](#12-several-samples)
13. [View settings](#13-view-settings)
14. [Keys and mouse](#14-keys-and-mouse)
15. [Messages and limits](#15-messages-and-limits)
15a. [Pictures and links](#15a-pictures-and-links)
16. [Privacy and security](#16-privacy-and-security)

## 1. Starting

Open `lgv.html` in Chrome, Edge, Firefox or Safari (Chrome or Edge recommended: some features need
them, see 3.1). At start LGV loads the **GRCh38 / hg38** reference (sequence only, from the IGV
team's public storage) and opens on BRCA1. A short message confirms the reference is loaded.

- What opens at start is set in **View → At start open**: GRCh38, GRCh37, T2T-CHM13, the demo data,
  or nothing.
- Without an internet connection the reference cannot load; LGV then opens the built-in demo data
  (a synthetic contig `demo` with example variants) and says so.
- The interface is in English; Russian is under **View → Language**.

## 2. The screen

From top to bottom:

- **Top bar**: **LGV** (click for About: version, credits, third-party code and data, the
  diagnostic log), the **contig list** (choose a chromosome or contig to see it whole; chromosomes
  first, alternative and unplaced contigs under **Other contigs**), the locus field with **Go**, **← →** (back / forward), **− +** (zoom), **PNG** (save a picture) and **Link**
  (copy a link to this view) — section 15a, the width of the view, the loaded reference and samples,
  the status line (reads shown, load and draw times), **Open** and **View**. **H** hides the bar; a
  button in the corner brings it back.
- **Chromosome map**: the chromosome on screen, with a red box where the view is (section 6).
- **Panel headers**: only with several regions side by side (section 5).
- **Ruler** and **reference sequence**: coordinates are 1-based. The reference shows as
  colored blocks, and as letters when zoomed in far enough.
- **Gene track**: MANE transcripts (section 7).
- **Samples**: each with its coverage and its reads. A small label at the top of each sample's
  reads, also when there is only one, names it (file and sample; it may cover a row or two), with × at the right
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
  you choose, so LGV asks once — **Allow the folder**. The folder picker opens where the file lies;
  confirm it and LGV takes the index from there itself. Later files from the same folder open
  without asking, until the page is closed. The access is read-only and never stored. In Firefox
  and Safari this is not possible: choose the file together with its index, or the whole folder.
- **Choose a folder**: every sample in it is found and paired with its index. A folder with several
  samples first shows a list to pick from (the first four are ticked).
- **Drag and drop** files or a folder onto the screen works the same way.

In Chrome and Edge every file and folder picker opens in the folder used last (until the page is
closed), since data usually lie in the same place.

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
- **Your own**: a FASTA file from disk or by link. Without a `.fai`, LGV indexes the file in the
  browser in one pass (plain, bgzip or gzip; a progress message shows). With `.fai` (and `.gzi` for
  bgzip) it opens at once.
- Samples are matched to the reference by contig name (case, `chr` prefix and version suffix
  ignored), then by a length only one contig on each side has. A sample made against another build
  (a large contig with the same name and a different length, e.g. chr1 of GRCh37 against GRCh38) is
  refused with an explanation, because its coordinates would be wrong.
- Without a reference, the first sample's header gives the coordinates, and reference bases are
  rebuilt from MD tags where reads cover them.
- **CRAM** stores reads as differences from the reference it was made with. LGV checks each part of
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
bands, red = centromere) and the band names (p21.3, q11.21…) inside the bands wide enough for them;
other genomes get a plain bar. A **red box** marks where the view is; the
active panel's box is filled. With several panels the boxes are numbered like the panel headers,
and panels on different chromosomes get their chromosomes side by side. **Click or drag** on the map
to move the active panel there.

## 7. Genes (MANE)

Under the reference sits the gene track: MANE Select and MANE Plus Clinical transcripts on GRCh38
(19,373 transcripts, MANE v1.5), built into LGV — it works offline, and gene names you type never
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
- Lines that cannot be read are skipped and counted: a message says how many. When most lines
  cannot be read, the file is probably not of the format its name says (for example a table of
  `samtools idxstats` saved as `.bed`): the message says so and shows the first such line.
- Tracks stay open when another reference is loaded (they are positions, not reads).

**Track type switches.** When an annotation track is open, the top bar holds a switch per open
kind — `✓ BAM 3`, `✓ VCF`, `✓ BED`, `✓ GFF`, `✓ bedGraph`, `✓ TDF`, `✓ bigWig` (BAM covers CRAM;
the number is how many). A click hides every track of that kind, in all panels, and a second click
brings them back as they were. With BAM off no reads or coverage are read at all, so a whole
chromosome opens at once and the other tracks — a BED of regions with loss of heterozygosity, say —
can be looked at on that scale. Opening a new sample turns BAM back on. The switches last until the
page is closed.

## 7b. Signal tracks: TDF and bigWig

Files of values along the genome — coverage, allele fractions, scores — open like the other tracks
(**Open → From disk** or **By link**) by their name: `.tdf` (made by igvtools) and `.bw` / `.bigwig`.
They are not read whole: only the part of the file a view needs is read, so a file of any size
opens at once. Over a link the server must allow Range requests.

**Several tracks in one file** (a TDF can hold several, for example one per sample) are shown as
separate rows, each with its name at the top left; up to 40 rows, the name line says if there are
more.

**What is read for a view.** Both formats keep the values themselves and precomputed summaries for
coarser scales. Zoomed in, LGV draws the values themselves. Zoomed out it reads the values
themselves as long as the view needs no more than 16 MB of them (4 MB over a link); beyond that it
uses the summaries, and the name line says so: for example `mean per 114 kb window`. A file may have
summaries only for the coarsest scales (igvtools can be told to stop early); between them and the
values themselves the coarsest summaries the file has are used. If the file has
no summary of what is asked for (a TDF made with the mean only, while the windowing function is
maximum), the name line names what is shown instead. When neither fits, the line says `zoom in`.

**How a track is drawn.** The settings of a TDF come from its track line (the `track ...` line of the
wiggle or bedGraph it was made from; a bigWig has none). Every setting of the UCSC wiggle track line
is followed:

| Setting | Meaning | When the file does not say |
| --- | --- | --- |
| `graphType=bar\|points` | bars or points (`line` and `heatmap` are also read) | see below |
| `windowingFunction=maximum\|mean\|minimum\|mean+whiskers` | which value a screen column shows when several fall into it | maximum |
| `autoScale=on\|off`, `viewLimits=lo:hi` | scale from the values in view, or fixed | auto |
| `alwaysZero=on` | the auto scale always includes 0 | off |
| `transformFunc=LOG` | log scale: ln(1 + x), −ln(1 − x) below zero | off |
| `smoothingWindow=2…16` | a moving mean over that many screen columns | off |
| `yLineMark=v`, `yLineOnOff=on` | a dashed line at value v | off |
| `color=r,g,b`, `altColor=r,g,b` | color of values ≥ 0 and of negative values | theme colors |
| `maxHeightPixels=max:default:min` | row height (the default) | 40 px (24 px with more than 3 rows) |
| `visibility=dense` | a heatmap: shade by value | — |

Windowing functions, as UCSC draws them: **maximum** shows the value of the column farthest from zero,
**minimum** the nearest one, **mean** the mean over the bases with data; **mean + whiskers** draws
the mean dark, mean ± 1 SD medium and the extreme value light. A point is one screen pixel; while the points in view are few (fewer than 0.8 × the track's
width in pixels) they are 2 × 2 px, easier to see.

**When a TDF says nothing about how to draw it**: a file of type `ALLELE_FREQUENCY` is drawn as
points; so is a file whose values all lie in 0–1 at separate positions (they cover less than half
of the stretch they span, as variant allele fractions do), on a fixed 0–1 scale. Everything else is
drawn as bars, the UCSC default. A bigWig without settings is drawn by the UCSC defaults (bars,
maximum, auto scale).

**Points from the values themselves show every value**, also when several fall into one screen
column, so a distribution (for example allele fractions around 0.5 and 1) stays visible at any
scale. Points drawn from summaries show one point per window, at its value.

- **Right click** on a signal track: draw as bars, points, line or heatmap; windowing function;
  scale (auto, 0–1, the file's `viewLimits`); always include zero; log scale; row height. The choice
  is remembered for that file in this browser; **Back to the file's settings** forgets it.
- **Click** on a row: a card with the track name and the value at that position (the nearest value
  within a few pixels), and for a screen column with several values or a summary also their mean,
  minimum, maximum and the number of bases with data.
- The name line shows the scale of a single-track file (`0–1`); with several rows each row shows its
  own.
- Contig names match with or without `chr`, as everywhere in LGV.

## 8. Coverage

Each sample has a coverage track: gray bars of read depth; the scale is in the corner (`0–31×`).
Where non-reference bases make up at least 20% of reads (**View → Highlight alleles from**), the
bar shows colored parts by base.

- **Reads are loaded in windows up to 200 kb** (View → Load reads up to, 1–1000 kb). The window
  loaded is a little wider than the screen. In a wider window only the depth is read (quickly, no
  allele colors), the reads area says why, and the message at the bottom offers **Load reads** with
  the number of reads to expect and, once each file has been read before, the expected time (from
  how fast that file was read so far); after that, windows on that contig up to 1.5 times as wide
  load their reads without asking. How many reads come in is predictable: depth × window ÷ read
  length, and with downsampling at most reads per window × window ÷ its width.
- Coverage counts **every** read that passes the filters, also where the reads shown are
  downsampled.
- Above 1 Mb wide, coverage is counted as depth only (no allele colors): it loads several times
  faster, and alleles cannot be seen at that scale anyway.
  Reads are not read at all in such windows (each is below a pixel, and the million reads of an
  exome would take minutes to read, send and lay out); a hint to zoom in stands in their place. The
  same past 300,000 reads in a window. If you zoom in while a wide window is being read, that read is
  stopped and the new view does not wait for it. A wide window is read in parts by several of the
  browser's worker threads at once, and its coverage is kept as the highest depth per short bin
  (16–64 bp), so moving around it draws at once.
- **Above 1 Mb the coverage first appears at once as an estimate from the index** (nothing is read:
  the index tells how many bytes of reads start in each 16 kb, which goes with depth). The scale reads
  `≈ 0–31×`; it is converted into depth by the windows of this file already seen in detail (before
  any, it is in relative units). It is the mean over 16 kb, not the peak. Up to 10 Mb the exact depth
  is then read and replaces it; wider views keep the estimate. Clicking the coverage still gives the
  exact counts at a base.
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
  secondary, QC-failed) do apply, as in the view. The counted column is outlined with a thin gray box.

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
- **MAPQ**: the body of a read with low mapping quality is transparent (50% at MAPQ 0, less and
  less up to solid at MAPQ 40 and above; adjustable up to 100%, an empty read): the read's color
  over the background, never lighter than the background. Every read has the same thin pale gray
  outline at any MAPQ, so a pale or empty read still shows where it lies, even zoomed out.
  Mismatches, insertions and clip marks are not paled: their color depends only on base quality.
- **Mismatches**: only bases that differ from the reference are colored: A green, C blue, G orange,
  T red; N is not marked. Letters appear inside reads from 7 pixels per base when rows are at least
  9 px tall.
- **Base quality**: below Q5 a mismatch is drawn at 20% strength, rising to full at Q20 (can be
  turned off: View → Fade low-quality bases).
- **Zoomed out**, lone mismatches fade so that sequencing errors do not clutter the picture;
  positions where a non-reference base passes the allele threshold always stay bright.
- **Insertions**: a purple bar, with its length written above from 10 pixels per base. **Click the
  bar** (within a few pixels) for the insertion's card: its place between two bases
  (`chr1:21,001^21,002`), length, the inserted bases, their qualities and the read name; the mark is
  outlined. **Copy** copies the inserted bases as FASTA: `>name/2_insertion_chr1:21001^21002`.
  **Deletions**: a dark line. **Spliced gaps (N)**: a thin gray line.
- **Soft clips** (shown by default; right click → Show soft-clipped bases): the clipped bases
  continue the read in its own color, slightly thinner so the clip boundary shows; clipped bases
  that differ from the reference are colored like mismatches. With clips hidden, a pink tick marks
  where a read is clipped.
- **Pairs** (right click → View as pairs): both mates of a pair in one row, joined by a thin gray line (the same gray at any read color and MAPQ);
  where mates overlap the overlap is shaded. A mate outside the loaded window is shown by the line
  to where it lies.
- **Downsampling**: in deep regions at most 100 reads are kept per 50 bp window (adjustable); marks above the reads show where reads were left out. Coverage still counts every read. Which reads are kept is decided by their names, so it is the same on every load, and **pairs stay together**: when a read is kept, its mate in the loaded window is kept too (a window may then hold up to twice the number). A mate opened from a card (**Go to mate**, **Mate in split screen**) is always shown.
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
  **Copy** copies the text. A mate in the loaded window that downsampling left out is named as such
  in the card; both buttons show it.
- **×** at the top right closes the card (also Escape).
- **Drag the card by its top row** to move it; later cards open where you left it (until the page
  is closed), even if part of it is off screen.

## 11. Right-click menu

Right click on the reads:

- **View mate region in split screen** (on a paired read with a mapped mate).
- **Copy read sequence**, **Copy read sequence, reverse complement**, **Copy left / right
  soft-clipped sequence**, **Copy read pair (FASTA)**. All are copied as FASTA with the read name:
  `>name/1`; the reverse complement as `>name/1_revcomp`; clips as `>name/1_left_clip` and
  `>name/1_right_clip`; a pair as `>name/1` and `>name/2`. The sequence is as stored in the file, on
  the reference's forward strand (the reverse complement: the opposite strand, read 5′ to 3′). Items that do not apply are grayed out.
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
- the scrollbar on the right, in its own strip (the panels end before it, nothing lies under it);
  PageUp / PageDown.

## 13. View settings

**View** opens the settings; they are kept in the browser as a **display profile**, which can be
copied as JSON, pasted from someone else and applied, or reset.

- **Reads**: read height (1–24 px), color, tag for color, fade low-quality bases, show soft clips,
  MAPQ transparency (from / to MQ, transparency at low MQ: 0% is a solid read, 100% an empty one
  with only the outline; the same at every scale: a read narrower than 3 pixels keeps its top and
  bottom lines). Every read has a thin pale gray outline; the pointed 3' end is drawn only
  in views up to 5 kb wide, further out reads are rectangles. Mismatches, insertions and clips are
  never paled by MAPQ: their color depends only on the base's own quality.
- **Group and sort**: group by strand, first / second of pair, pair orientation or a tag; sort rows.
- **Filters**: the read window (load reads up to, kb, see section 8); minimum MAPQ; hide duplicates, secondary, supplementary, QC-failed reads;
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
| Picture of the view | PNG button |
| Link to the view | Link button |
| About LGV | click LGV in the top bar |
| Read details | click a read |
| Inserted bases | click an insertion mark |
| Counts at a base | click the coverage |
| Transcript details | click the gene track |
| Menu | right click on the reads |
| Sort by base at the center | S |
| Hide / show the top bar | H |
| Close a card, a menu or a message | Escape |

## 15. Messages and limits

- Windows that would read more than 64 MB of compressed data from links, or 512 MB from files on
  disk (they read much faster), ask first, with the size; **load anyway** reads it. After that, the
  same contig is not asked about again in this panel for windows up to 1.5 times as large. With long
  reads (nanopore, PacBio) every window of a short contig, such as a plasmid, needs about the same
  data, because the reads span it whole: zooming in does not make the size smaller.
- With long reads that overlap, the reads stack thousands of rows deep; only the rows around the
  visible ones are prepared for drawing, the rest as you scroll to them.
- A damaged file gives a clear message saying what is wrong (for example a checksum mismatch in a
  block) instead of showing wrong data.
- Reading that takes longer than 60 s is stopped with a message; the file may be damaged. Each
  sample is read one request at a time, and only the newest view of each panel waits, so fast
  panning and zooming do not pile up reads.
- **About → Copy diagnostic log** copies what took long in this session (sample reads with their
  times, slow frames, the page being busy, background reader restarts): send it with a report of a
  freeze or a slow step. The log stays in the page only; it holds file names and loci, no reads.
- Messages at the bottom close with × or a click; short ones go by themselves.
- **No index found**: choose the index together with the file, the whole folder, or allow the folder
  when asked (3.1).
- **Another genome build**: the sample was aligned to a different reference than the one loaded
  (3.3).

## 15a. Pictures and links

- **PNG** saves the view as it is on screen — all panels, the chromosome map, genes, tracks,
  coverage and reads — at screen resolution, with a line at the bottom naming the loci, the
  samples, the reference and the date. The file is named after the first locus, for example
  `LGV_chr7_117559500-117560000.png`. Cards and menus are not in the picture.
- **Link** copies a link that opens this view: the reference, the samples and tracks that were
  opened by link, and every panel's locus. Files opened from this computer cannot be in a link (a
  message names them), and access tokens are never put in a link.
- **Opening LGV from another system** (a LIMS, a report, a variant table): add parameters after
  `#` (or `?`) to the address of `lgv.html`:

  ```
  lgv.html#ref=hg38&file=https://server/sample.cram&locus=CFTR
  ```

  - `ref` — `hg38`, `hg19`, `chm13`, `demo`, or a link to a FASTA (`.fai` next to it, `.gzi` for
    `.fa.gz`); without it, the reference chosen in View → At start open;
  - `file` — a BAM, CRAM or annotation track link; repeat for several;
  - `index` — the BAM / CRAM index links in the same order, when they are not next to the files;
  - `locus` — loci or genes; several separated by spaces (`%20` or `+`) open side by side.

  Only `https` links are read (and `http` on this computer). Parameters after `#` stay in the
  browser and are not sent to the server that hosts `lgv.html`, so they are the better choice for
  patient data.

## 16. Privacy and security

- Data goes only between the files or links you chose and your computer. LGV sends nothing
  anywhere; the only connection it makes by itself is loading the reference chosen to open at
  start (sequence only).
- Gene search and the gene track are built in: typed gene names are not sent anywhere.
- Only the display profile is saved in the browser. Tokens, links and folder access are never
  saved.
- Every file is treated as untrusted: damaged or hostile files end in a clear message, never in
  wrong data, a hang or a crash. Details are in [SECURITY.md](../SECURITY.md).
