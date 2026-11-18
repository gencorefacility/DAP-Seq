#!/usr/bin/env python3
"""GEM QC for one peak set: MultiQC custom-content JSON.

Reads GEM's output dir (<sample>/):
  <sample>.GEM_log.txt   GEM's stdout (the GEM process keeps it): read counts before and after
                         GEM's duplicate filter, as "<cond>_IP  Bases: n  HitCounts: n" lines

Writes:
  <sample>.gem_reads_mqc.json   General Statistics: reads GEM used, % its filter removed

Standard library only.
"""

import argparse
import json
import os
import re

# A count line, e.g. "sample1_IP  <tab>Bases: 900<tab>HitCounts: 1000.0"
#   group 1 = GEM condition (sample1), 2 = IP or CTRL, 3 = HitCounts (reads; Bases = distinct
#   read-start positions, not used)
COUNT = re.compile(r"(\S+)_(IP|CTRL)\s+Bases:\s*\d+\s+HitCounts:\s*([\d.]+)")


def read_counts(path, condition):
    """{"IP"/"CTRL": (before, after)} for one GEM condition; lines of any other condition are
    ignored. "after" = the lines following "Filter duplicate reads."; without that line (--nrf)
    nothing was removed."""
    # before / after GEM's duplicate filter: {"IP": reads, "CTRL": reads}. filtered = whether the
    # "Filter duplicate reads." line, which GEM prints between the two sets of counts, was seen yet
    before, after, filtered = {}, {}, False
    with open(path) as fh:
        for line in fh:
            if line.startswith("Filter duplicate reads"):
                filtered = True
            m = COUNT.match(line)
            # count lines go into before until that line, into after from then on
            if m and m[1] == condition:
                (after if filtered else before)[m[2]] = float(m[3])
    # no count lines: stop, so a missing count is never reported as 0
    if not before:
        raise SystemExit(f"no read counts for condition {condition} "
                         f"('{condition}_IP  Bases: ...  HitCounts: ...') in {path}")
    # no filter line: GEM ran with --nrf, so nothing was removed (0%)
    if not filtered:
        after = before
    # e.g. {"CTRL": (2000.0, 1500.0), "IP": (1000.0, 900.0)}; no CTRL without a control
    return {kind: (before[kind], after[kind]) for kind in sorted(before)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--dir", required=True, help="GEM's output dir (<sample>/)")
    args = ap.parse_args()

    # The GEM process names the condition after the peak set (--expt<id>), so it is --sample
    log = os.path.join(args.dir, f"{args.sample}.GEM_log.txt")
    counts = read_counts(log, args.sample)
    # MultiQC columns, named <key>_reads / <key>_removed (they must match "headers" below):
    # IP -> gem_reads, gem_removed; CTRL -> gem_ctrl_reads, gem_ctrl_removed
    stats = {}
    for kind, key in (("IP", "gem"), ("CTRL", "gem_ctrl")):
        if kind in counts:
            before, after = counts[kind]
            stats[f"{key}_reads"] = int(after)
            stats[f"{key}_removed"] = round(100 * (before - after) / before, 2) if before else 0

    filt = ("GEM's duplicate filter caps how many reads may start at the same position, based on "
            "how many reads are at nearby positions.")
    gem_reads = {
        "id": "gem_reads", "section_name": "GEM reads", "plot_type": "generalstats",
        "data": {args.sample: stats},
        "headers": {
            "gem_reads": {"title": "GEM reads", "format": "{:,.0f}", "scale": "Oranges",
                          "description": "Experiment reads GEM called peaks from: the filtered "
                                         "BAM(s) after GEM's duplicate filter. " + filt},
            "gem_removed": {"title": "GEM dups", "suffix": "%", "min": 0, "max": 100, "scale": "OrRd",
                            "description": "% of the experiment's filtered reads that GEM's "
                                           "duplicate filter removed"},
            "gem_ctrl_reads": {"title": "GEM ctrl reads", "format": "{:,.0f}", "scale": "Oranges",
                               "hidden": True,
                               "description": "Control reads GEM used, after its duplicate filter"},
            "gem_ctrl_removed": {"title": "GEM ctrl dups", "suffix": "%", "min": 0, "max": 100,
                                 "scale": "OrRd", "hidden": True,
                                 "description": "% of the control's filtered reads that GEM's "
                                                "duplicate filter removed"},
        },
    }
    with open(f"{args.sample}.gem_reads_mqc.json", "w") as fh:
        json.dump(gem_reads, fh, indent=1)
    # log line, e.g. "sample1: CTRL 2,000 -> 1,500, IP 1,000 -> 900"
    # (k = IP/CTRL, b = before, a = after)
    print(f"{args.sample}: " + ", ".join(f"{k} {int(b):,} -> {int(a):,}" for k, (b, a) in counts.items()))


if __name__ == "__main__":
    main()
