# Reading the DAP-seq MultiQC report

This guide explains three sections of `multiqc_report.html`: General Statistics, Motif presence and MEME-ChIP motifs. It applies to any report from this pipeline.

The examples here use made-up names and numbers: protein `TF1` tested on DNA library `Lib` in replicates `rep1` and `rep2`, with control `Ctrl_Lib_rep1`.

Two terms are used throughout:
- **Library:** one sequenced sample (a pair of FASTQ files if paired-end).
- **Peak set:** the peaks from one run of the peak caller, GEM.

## General Statistics

The table covers three kinds of output: FASTQ files (the sequenced reads), BAM files (those reads aligned to the genome) and peak sets (the binding sites called from the aligned reads). Each kind fills its own columns, so a blank cell means that column doesn't apply to that row, not that a step failed.

MultiQC merges everything that shares a name into one row, so one row can hold more than one kind of output. That's why `TF1_Lib_rep1` has both peak and alignment columns: the replicate's peak set and its unfiltered BAM have the same name. Its peaks were called from the filtered reads, but their numbers are on this row, not on `TF1_Lib_rep1.filtered`.

### What each row is

| Row | What it is | Columns filled |
|---|---|---|
| `TF1_Lib` | **Pooled peak set:** GEM run on all replicates together, with all of their controls. It isn't a library, so it has no trimming or alignment numbers. | Peaks, FRiP, GEM reads, GEM dups |
| `TF1_Lib_rep1` | **One replicate:** the peak set called from this replicate alone, plus the unfiltered BAM with all of the library's alignments. | The columns above, % Aligned, and samtools on all alignments |
| `TF1_Lib_rep1.filtered` | **Filtered BAM:** the same alignments after the filters below. Peak calling, FRiP and coverage tracks use these reads, but the peak numbers are on the `TF1_Lib_rep1` row. | samtools only |
| `TF1_Lib_rep1_1` and `_2` | **Read 1 and read 2 FASTQ files**, after trimming. | Trimmed bases, Dups, GC, Seqs |
| `Ctrl_Lib_rep1` and its variants | **Control library:** the same rows as a replicate, but no peak columns, because controls aren't peak-called. | samtools, % Aligned, FASTQ columns |

The filters are mapping quality of at least 30 (the default), the chromosomes chosen for the run, and, for paired-end libraries, proper pairs in which both mates pass, of which only read 1 is kept.

Some reports look slightly different:
- **Single-end libraries** have no `_1` and `_2` rows. Their trimming and FastQC numbers sit in the `TF1_Lib_rep1` row.
- **The pooled row** appears only when replicates were marked as a group. Controls never have one.

### Columns in one row can count different reads

In a replicate's row, the samtools columns count both mates with no filtering, while FRiP and GEM reads use the filtered reads. For a paired-end library, the counts connect like this:

| Stage | Where to find it | Example |
|---|---|---|
| Read 1 and read 2 after trimming | `_1` and `_2` rows, **Seqs** | 2.0M each |
| All alignment records, both mates | `TF1_Lib_rep1`, **Total seqs** | 4.0M |
| After the filters, read 1 of each passing pair | `TF1_Lib_rep1.filtered`, **Total seqs** | 1.4M |
| After GEM's duplicate filter | `TF1_Lib_rep1`, **GEM reads** | 1.2M |

The filtered count is below half of Total seqs because it keeps one read per DNA fragment and also drops fragments that fail the other filters, including pairs where one mate didn't align or the mates didn't align as a proper pair. GEM has no paired-end mode and would otherwise count each fragment twice.

### What the columns mean

**Peak-set columns**

| Column | Meaning |
|---|---|
| Peaks | Number of peaks (binding events) GEM called. |
| FRiP | Fraction of reads in peaks: % of the filtered reads that overlap a peak. GEM's peaks are narrow (201 bp), so this runs lower than with broad-peak callers. |
| GEM reads | Reads GEM called peaks from: the filtered reads after GEM's own duplicate filter. |
| GEM dups | % of the filtered reads that filter removed. It caps each position at the read count expected from the surrounding ~20 bp, so it removes most duplicates in background and few in strong peaks. |

**Alignment columns**

| Column | Meaning |
|---|---|
| Total seqs | Records in the alignment file, in millions. On a paired-end replicate's row this counts both mates, including unaligned ones. |
| Reads mapped | Reads that aligned, in millions. It appears twice because two samtools commands report it; the values match. |
| % Mapped | Reads mapped as a % of Total seqs. Always 100% on `.filtered` rows, because unaligned reads were removed. |
| % Aligned | Bowtie2's overall alignment rate, before any filter. |
| % Proper pairs | % of reads aligned with their mate at the expected orientation and distance (paired-end only). |
| Error rate | Mismatches per aligned base. |
| Non-primary | Extra alignments reported for a read, in millions. |

**FASTQ columns**

| Column | Meaning |
|---|---|
| Seqs | Reads in this FASTQ file after trimming, in millions. |
| Trimmed bases | % of bases removed by adapter and quality trimming. |
| Dups | % of reads whose sequence is duplicated (FastQC). This is separate from GEM dups, which works on aligned positions. |
| GC | Average GC content, in %. |

The **Configure columns** button above the table reveals more: GEM ctrl reads and GEM ctrl dups (the control's numbers for this peak set, on the experiment's row), % MapQ 0 reads, Mean insert, Reads, % Reads mapped, Avg len, Median len and Failed.

## Motif presence

This table asks how many peaks contain a given motif, across all peaks in the set. Each row is one peak set and one motif, such as `TF1_Lib_rep1 CACGTRG`. Rows are sorted alphabetically, not by significance.

Each peak is scanned in a 201 bp window: its summit ± 100 bp. Peaks too close to a chromosome end for a full window are left out.

Each window is then scanned a second time after shuffling, to show how often the motif turns up by chance. Shuffling scrambles the order of the window's bases, which breaks up real binding sites, so matches found in shuffled windows are mostly chance. The motif itself isn't changed. The shuffle keeps the count of every pair of neighbouring bases (the "dinucleotide content")

### The three kinds of motif

| Name looks like | Where it comes from | How it's matched |
|---|---|---|
| A consensus, such as `CACGTRG` | Given in the config as a motif we expect. It can use IUPAC codes (`R` = A or G, `Y` = C or T, `N` = any). | Exact match on either strand, with no mismatches. |
| A motif name, such as `TF1_matrix` | A position weight matrix provided via config. The name comes from its motif file. | FIMO at p < 1e-4, on either strand. |
| Letters ending in a tool and rank, such as `NCACGTGN-MEME-1` | The top-ranked motif from this peak set's own MEME-ChIP run. | FIMO at p < 1e-4, on either strand. |

The letters in the third kind are only a label; FIMO scans the full matrix. That row is whichever motif MEME-ChIP ranked first, so it could be a known motif from the database or a de novo motif.

### What the columns mean

"Strongest" means ranked by q-value, then by signal. Every percentage except 2+ copies is the % of those peaks with at least one match.

| Column | Meaning |
|---|---|
| Peaks | Peaks scanned. |
| Top 100 | % of the 100 strongest peaks with at least one match. |
| Top 1,000 | % of the 1,000 strongest peaks with at least one match. |
| All peaks |  % of all peaks with at least one match. |
| Weakest 10% | % of the weakest tenth of the peaks (just above the caller's cutoff) with at least one match. |
| Shuffled |  The same as All peaks, but scanning the shuffled windows instead of the real ones. This is the chance level. |
| 2+ copies | % of all peaks with two or more matches. |

Every peak's matches, with copies, offset from the summit and strand, are in `motif_qc/<peak set>.motif_peaks.tsv` in the results folder.

## MEME-ChIP motifs

Each row summarises one full MEME-ChIP run on one peak set: its best de novo motif and its best known motif. Every motif from the run is in `meme_chip/<peak set>/meme-chip.html` in the results folder.

### How the table is generated

1. **Pick the peaks.** Peaks are ranked by q-value, then signal, and the top ones are kept (600 by default).
2. **Cut the sequences.** Each peak's summit ± 100 bp (201 bp) is taken from the reference genome.
3. **Find new motifs.** MEME and STREME search the central 100 bp of each sequence for motifs 6 to 15 bp wide, without using any known motif.
4. **Compare with known motifs.** Tomtom compares each new motif with a motif database chosen for the run. Separately, CentriMo scans the full 201 bp with each database motif and tests whether it concentrates at the peak centres.
5. **Rank everything.** MEME-ChIP lists all these motifs together, sorted by E-value.
6. **Fill the table.** The report takes the first de novo motif and the first known motif from that list.

### What the columns mean

| Column | Meaning |
|---|---|
| Sequences | Sequences given to MEME-ChIP (the top peaks). |
| De novo motif | Logo of the most significant de novo motif. Hover over it to see its consensus, which is also a hidden Consensus column. |
| Found by | The tool that found it: MEME or STREME. |
| Sites | Sequences containing the motif, as counted by that tool. |
| E-value | The de novo motif's significance. Smaller is more significant. |
| Closest known motif | The database motif most similar to the de novo motif (Tomtom). `none` means no match at E <= 1, or no database was used. |
| Top known motif | Logo of the database motif most enriched at the peak centres (CentriMo). |
| Known motif ID | That motif's ID and name in the database. |
| Known E-value | CentriMo's E-value for that motif. |

