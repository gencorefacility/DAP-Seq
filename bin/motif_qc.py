#!/usr/bin/env python3
"""Motif QC for one sample: MultiQC custom-content JSON and two TSV tables.

Inputs are summit-centred peak windows (all equal width, FASTA header = peak rank), the same windows
dinucleotide-shuffled (header = rank + "_shuf") and the ranked peaks they were cut from (line number =
rank). Motifs come from two places, in this order:
  --consensus  IUPAC strings (e.g. CAGGTAG): exact matches on either strand, scanned here
  --fimo P     motif matrices scanned with FIMO --text: P.motifs.tsv (id, alt name), P.real.tsv and
               P.shuffled.tsv (repeatable; one per motif file, in report order)

Writes, per sample (MultiQC merges files that share an id into one table):
  <sample>.motif_summary_mqc.json    per motif: % of peaks with >= 1 match (top 100, top 1,000, all,
                                     weakest 10%), shuffled, 2+ copies
  <sample>.memechip_mqc.json         MEME-ChIP's top de novo and top known motif (from its summary.tsv),
                                     with sequence logos drawn from combined.meme
  <sample>.motif_summary.tsv         the summary table, plus peaks with 0 / 1 / 2 / 3+ copies
  <sample>.motif_peaks.tsv           one row per peak: coordinates, and per motif its copies and their
                                     offsets from the summit (with strand)

Standard library only.
"""

import argparse
import html
import json
import math
import os
import re

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "Y": "[CT]", "S": "[CG]", "W": "[AT]",
         "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]"}
COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")

# Sequence logos: letter outlines in a 100 x 100 box (y down), stretched to each letter's height
LOGO_GLYPHS = {
    "A": "M0 100L38 0H62L100 100H77L68 74H32L23 100ZM39 55H61L50 24Z",   # even-odd: the hole
    "C": "M88.8 17.9A48 50 0 1 0 88.8 82.1L73.4 69.3A28 30 0 1 1 73.4 30.7Z",
    "G": "M88.8 17.9A48 50 0 1 0 99.8 52H56V64H76.8A28 30 0 1 1 73.4 30.7Z",
    "T": "M0 0H100V20H61V100H39V20H0Z",
}
LOGO_COLOURS = {"A": "#109648", "C": "#255c99", "G": "#f7b32b", "T": "#d62839"}


def peak_rank(name):
    return int(name[:-len("_shuf")] if name.endswith("_shuf") else name)


def read_fasta(path):
    """[(rank, sequence)]; the header is the peak rank (shuffled: rank + "_shuf")."""
    seqs = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                seqs.append([peak_rank(line[1:].split()[0]), []])
            elif line:
                seqs[-1][1].append(line.upper())
    return [(r, "".join(s)) for r, s in seqs]


def scan_consensus(seqs, motif):
    """Exact matches on either strand, as [(rank, match centre (1-based), strand)] like read_fimo.
    A match on both strands at the same place (palindromes) counts once, strand '.'."""
    fwd, rev = (re.compile("(?=" + "".join(IUPAC[c] for c in p) + ")")
                for p in (motif, motif.translate(COMPLEMENT)[::-1]))
    hits = []
    for rank, seq in seqs:
        f = {m.start() for m in fwd.finditer(seq)}
        r = {m.start() for m in rev.finditer(seq)}
        strand = {s: "." if s in r else "+" for s in f}
        strand.update((s, "-") for s in r - f)
        hits += [(rank, s + (len(motif) + 1) / 2, strand[s]) for s in sorted(strand)]
    return hits


def read_fimo(path):
    """{motif_id: [(rank, match centre (1-based), strand), ...]}"""
    hits = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6 or f[0].startswith(("motif_id", "#")):
                continue
            hits.setdefault(f[0], []).append((peak_rank(f[2]), (int(f[3]) + int(f[4])) / 2, f[5]))
    return hits


def read_ranked(path):
    """{rank: narrowPeak fields}; rank = line number, as in the windows' headers"""
    with open(path) as fh:
        return {i: line.rstrip("\n").split("\t") for i, line in enumerate(fh, 1)}


def read_meme_matrices(path):
    """{motif id: [[pA, pC, pG, pT] per position]} from a MEME-format motif file"""
    matrices, motif, rows = {}, None, 0
    if not os.path.exists(path):
        return matrices
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if f[:1] == ["MOTIF"]:
                motif, rows = f[1], 0
            elif line.startswith("letter-probability") and motif:
                rows = int(re.search(r"w=\s*(\d+)", line)[1])
                matrices[motif] = []
            elif rows and len(f) == 4:
                matrices[motif].append([float(x) for x in f])
                rows -= 1
    return matrices


def logo(matrix, label, height=30):
    """A sequence logo as inline SVG, so the report needs no image files. One column per position, as
    tall as its information content (at most 2 bits); each base's share of the column is its
    probability, tallest on top. No commas (MultiQC turns commas in table cells into spacer <span>s)
    and no element ids (many logos share one page); the label is the hover text."""
    col, full = 10, 24   # SVG units: column width, height of 2 bits
    letters = []
    for i, probs in enumerate(matrix):
        bits = 2 + sum(p * math.log2(p) for p in probs if p > 0)
        y = full
        for p, base in sorted(zip(probs, "ACGT")):
            h = p * bits * full / 2
            if h < 0.3:   # under half a pixel at this size
                continue
            y -= h
            letters.append(f'<path d="{LOGO_GLYPHS[base]}" fill="{LOGO_COLOURS[base]}" fill-rule="evenodd" '
                           f'transform="translate({i * col + 0.5:g} {y:.2f}) scale({(col - 1) / 100:g} {h / 100:.4f})"/>')
    width = len(matrix) * col
    return (f'<svg xmlns="http://www.w3.org/2000/svg" role="img" width="{round(height * width / full)}" '
            f'height="{height}" viewBox="0 0 {width} {full}"><title>{html.escape(label.replace(",", " "))}</title>'
            f'{"".join(letters)}</svg>')


def memechip_row(directory):
    """MEME-ChIP's top de novo motif (MEME/STREME/DREME) and top known motif (CentriMo, from -db) as a
    table row. summary.tsv lists every motif, most significant first; E-values stay strings
    (1.4e-1742 is below the smallest float). Logos come from combined.meme, whose MOTIF ids are
    summary.tsv's MOTIF_INDEX; a motif without a matrix there shows its consensus instead."""
    with open(f"{directory}/summary.tsv") as fh:
        lines = [line.rstrip("\n").split("\t") for line in fh if line.strip() and not line.startswith("#")]
    motifs = [dict(zip(lines[0], f)) for f in lines[1:]] if lines else []
    with open(f"{directory}/seqs-centered") as fh:   # MEME-ChIP's input sequences, centred
        row = {"sequences": sum(line.startswith(">") for line in fh)}
    matrices = read_meme_matrices(f"{directory}/combined.meme")

    def motif_logo(m):
        matrix = matrices.get(m["MOTIF_INDEX"]) or matrices.get(m["MOTIF_ID"])
        return logo(matrix, m["CONSENSUS"]) if matrix else m["CONSENSUS"]

    denovo = [m for m in motifs if m["MOTIF_SOURCE"] in ("MEME", "STREME", "DREME")]
    known = [m for m in motifs if m["MOTIF_SOURCE"] not in ("MEME", "STREME", "DREME")]
    if denovo:
        m = denovo[0]
        row.update(denovo=motif_logo(m), denovo_consensus=m["CONSENSUS"], denovo_tool=m["MOTIF_SOURCE"],
                   denovo_sites=int(m["SITES"]),
                   denovo_evalue=m["E-VALUE"], denovo_match=m.get("MOST_SIMILAR_MOTIF", "").strip() or "none")
    if known:
        m = known[0]
        alt = m.get("ALT_ID", "").strip()
        row.update(known_logo=motif_logo(m),
                   known=f"{m['MOTIF_ID']} ({alt})" if alt and alt != m["MOTIF_ID"] else m["MOTIF_ID"],
                   known_evalue=m["E-VALUE"])
    return row


def custom_content(cc_id, section, description, plot_type, data, pconfig=None, headers=None):
    d = {"id": cc_id, "section_name": section, "description": description,
         "plot_type": plot_type, "data": data}
    if pconfig:
        d["pconfig"] = pconfig
    if headers:
        d["headers"] = headers
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--consensus", action="append", default=[], help="IUPAC consensus (repeatable)")
    ap.add_argument("--fimo", action="append", default=[], help="prefix of a FIMO scan set (repeatable)")
    ap.add_argument("--windows", required=True, help="peak-window FASTA, header = rank")
    ap.add_argument("--shuffled-windows", required=True, help="the same windows, shuffled")
    ap.add_argument("--peaks", required=True, help="the ranked narrowPeak the windows were cut from")
    ap.add_argument("--halfwidth", type=int, required=True)
    ap.add_argument("--top-n", type=int, default=0, help="peaks sent to MEME-ChIP (noted in the table)")
    ap.add_argument("--memechip", required=True, help="MEME-ChIP output dir (summary.tsv, seqs-centered)")
    args = ap.parse_args()

    # [(label, real matches, shuffled matches)]
    windows, shuffled = read_fasta(args.windows), read_fasta(args.shuffled_windows)
    motifs = [(c, scan_consensus(windows, c), scan_consensus(shuffled, c)) for c in args.consensus]
    for prefix in args.fimo:
        real, shuf = read_fimo(f"{prefix}.real.tsv"), read_fimo(f"{prefix}.shuffled.tsv")
        with open(f"{prefix}.motifs.tsv") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[0]:
                    alt = f[1] if len(f) > 1 and f[1] else f[0]
                    motifs.append((alt, real.get(f[0], []), shuf.get(f[0], [])))

    # Percentages are over the windows actually scanned (peaks too close to a contig end for a
    # full window are left out), strongest first
    ranks = sorted(r for r, _ in windows)
    n = len(ranks)
    tenth = -(-n // 10)   # n / 10, rounded up
    weakest = ranks[-tenth:] if n else []

    def pct(with_match, subset):
        return round(100 * sum(1 for r in subset if r in with_match) / len(subset), 2) if subset else 0

    summary, offsets = {}, []
    centre = args.halfwidth + 1   # summit, 1-based within the window
    with open(f"{args.sample}.motif_summary.tsv", "w") as fh:
        fh.write("peak_set\tmotif\tpeaks\tpct_top_100\tpct_top_1000\tpct_all\tpct_weakest_10pct\tpct_shuffled\t"
                 "pct_2plus_copies\tpeaks_0_copies\tpeaks_1_copy\tpeaks_2_copies\tpeaks_3plus_copies\n")
        for alt, real_hits, shuf_hits in motifs:
            per_peak = {}
            for r, c, strand in real_hits:
                per_peak.setdefault(r, []).append((c - centre, strand))
            offsets.append(per_peak)
            copies = [len(per_peak.get(r, [])) for r in ranks]
            with_shuf = {r for r, *_ in shuf_hits}
            row = {"peaks": n, "top_100": pct(per_peak, ranks[:100]), "top_1000": pct(per_peak, ranks[:1000]),
                   "all": pct(per_peak, ranks), "weakest_10": pct(per_peak, weakest),
                   "shuffled": pct(with_shuf, ranks),
                   "copies_2plus": round(100 * sum(1 for k in copies if k >= 2) / n, 2) if n else 0}
            summary[f"{args.sample} {alt}"] = row
            fh.write("\t".join(map(str, [args.sample, alt] + list(row.values()) + [
                sum(1 for k in copies if k == 0), sum(1 for k in copies if k == 1),
                sum(1 for k in copies if k == 2), sum(1 for k in copies if k >= 3)])) + "\n")

    # One row per ranked peak; motif columns are empty for a peak that had no full window
    scanned = set(ranks)
    with open(f"{args.sample}.motif_peaks.tsv", "w") as fh:
        fh.write("\t".join(["rank", "chrom", "start", "end", "summit", "signal", "qvalue"]
                           + [f"{alt}_{col}" for alt, *_ in motifs for col in ("copies", "offsets")]) + "\n")
        for r, f in sorted(read_ranked(args.peaks).items()):
            row = [r, f[0], f[1], f[2], int(f[1]) + int(f[9]) if int(f[9]) >= 0 else "", f[6], f[8]]
            for per_peak in offsets:
                hits = sorted(per_peak.get(r, []))
                row += [len(hits), ",".join(f"{o:g}({s})" for o, s in hits)] if r in scanned else ["", ""]
            fh.write("\t".join(map(str, row)) + "\n")

    top = (f" MEME-ChIP's motifs were learned from the top {args.top_n:,} peaks." if args.top_n else "")
    top_peaks = f"top {args.top_n:,} peaks" if args.top_n else "all peaks"
    pct_col = {"suffix": "%", "min": 0, "max": 100, "format": "{:,.1f}", "scale": "Purples"}
    out = {
        "motif_summary": custom_content(
            "motif_summary", "Motif presence",
            f"One row per peak set and motif: the % of peaks with at least one match within summit "
            f"+/- {args.halfwidth} bp. Consensus strings: exact match on either strand; motif matrices "
            "(--qc_motifs files, MEME-ChIP): FIMO match. Every peak's matches (copies, offset from the "
            "summit, strand): motif_qc/&lt;peak set&gt;.motif_peaks.tsv." + top,
            "table", summary,
            pconfig={"id": "motif_summary_table", "title": "Motif presence", "col1_header": "Peak set, motif"},
            headers={
                "peaks": {"title": "Peaks", "format": "{:,.0f}", "scale": False,
                          "description": "Peaks scanned (peaks too close to a contig end for a full "
                                         "window are left out)"},
                "top_100": dict(pct_col, title="Top 100", description="% of the 100 strongest peaks "
                                                                      "with at least one match"),
                "top_1000": dict(pct_col, title="Top 1,000", description="% of the 1,000 strongest peaks "
                                                                         "with at least one match"),
                "all": dict(pct_col, title="All peaks", description="% of all peaks with at least one match"),
                "weakest_10": dict(pct_col, title="Weakest 10%", description="% of the weakest tenth of "
                                   "the peaks (just above the caller's cutoff) with at least one match"),
                "shuffled": dict(pct_col, title="Shuffled", scale="Greys",
                                 description="% of the same sequences, dinucleotide-shuffled, with at "
                                             "least one match: the chance level"),
                "copies_2plus": dict(pct_col, title="2+ copies", description="% of peaks with two or "
                                                                             "more matches"),
            }),
        "memechip": custom_content(
            "memechip_motifs", "MEME-ChIP motifs",
            f"Motif discovery (MEME-ChIP) on the {top_peaks} (summit +/- {args.halfwidth} bp). De novo (MEME, "
            "STREME; central 100 bp of each sequence by default): the most significant motif, how many of the "
            "sequences contain it, and the known motif it most resembles (Tomtom against --memechip_db). "
            "Known (CentriMo; full sequences): the motif from --memechip_db most enriched at the peak centres. "
            "Hover a logo for its IUPAC consensus. All of MEME-ChIP's motifs: "
            "meme_chip/&lt;peak set&gt;/meme-chip.html.",
            "table", {args.sample: memechip_row(args.memechip)},
            pconfig={"id": "memechip_motifs_table", "title": "MEME-ChIP motifs", "col1_header": "Peak set"},
            headers={
                "sequences": {"title": "Sequences", "format": "{:,.0f}", "scale": False,
                              "description": "Sequences given to MEME-ChIP (the top peaks)"},
                "denovo": {"title": "De novo motif", "scale": False,
                           "description": "Logo of MEME-ChIP's most significant de novo motif"},
                "denovo_consensus": {"title": "Consensus", "scale": False, "hidden": True,
                                     "description": "The de novo motif's IUPAC consensus "
                                                    "(W = A/T, S = C/G, D = not C, N = any, ...)"},
                "denovo_tool": {"title": "Found by", "scale": False,
                                "description": "Which MEME-ChIP step found it (MEME or STREME)"},
                "denovo_sites": {"title": "Sites", "format": "{:,.0f}", "scale": "Purples",
                                 "description": "Sequences with a site, as reported by MEME-ChIP (MEME "
                                                "counts only the sequences it searched: -meme-searchsize, "
                                                "100,000 bp by default)"},
                "denovo_evalue": {"title": "E-value", "scale": False,
                                  "description": "The de novo motif's E-value"},
                "denovo_match": {"title": "Closest known motif", "scale": False,
                                 "description": "Most similar motif in --memechip_db (Tomtom); "
                                                "'none' = no significant match or no database"},
                "known_logo": {"title": "Top known motif", "scale": False,
                               "description": "Logo of the known motif (--memechip_db) most enriched at "
                                              "peak centres (CentriMo); blank without a database"},
                "known": {"title": "Known motif ID", "scale": False,
                          "description": "ID (and name) of that known motif in --memechip_db"},
                "known_evalue": {"title": "Known E-value", "scale": False,
                                 "description": "CentriMo E-value of the top known motif"},
            }),
    }
    for name, d in out.items():
        with open(f"{args.sample}.{name}_mqc.json", "w") as fh:
            json.dump(d, fh, indent=1)
    print(f"{args.sample}: {n} windows, motifs: {', '.join(m[0] for m in motifs) or 'none'}")


if __name__ == "__main__":
    main()
