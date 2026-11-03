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

COUNT = re.compile(r"(\S+)_(IP|CTRL)\s+Bases:\s*\d+\s+HitCounts:\s*([\d.]+)")


def read_counts(path):
    """{"IP"/"CTRL": (before, after)}, summed over conditions. "after" = the lines following
    "Filter duplicate reads."; without that line (--nrf) nothing was removed."""
    before, after, filtered = {}, {}, False
    with open(path) as fh:
        for line in fh:
            if line.startswith("Filter duplicate reads"):
                filtered = True
            m = COUNT.match(line)
            if m:
                (after if filtered else before)[(m[1], m[2])] = float(m[3])
    if not before:
        raise SystemExit(f"no read counts ('<cond>_IP  Bases: ...  HitCounts: ...') in {path}")
    if not filtered:
        after = before
    return {kind: (sum(v for (_c, k), v in before.items() if k == kind),
                   sum(v for (_c, k), v in after.items() if k == kind))
            for kind in sorted({k for _c, k in before})}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--dir", required=True, help="GEM's output dir (<sample>/)")
    args = ap.parse_args()

    counts = read_counts(os.path.join(args.dir, f"{args.sample}.GEM_log.txt"))
    stats = {}
    for kind, key in (("IP", "gem"), ("CTRL", "gem_ctrl")):
        if kind in counts:
            before, after = counts[kind]
            stats[f"{key}_reads"] = int(after)
            stats[f"{key}_removed"] = round(100 * (before - after) / before, 2) if before else 0

    filt = ("GEM's duplicate filter caps each position's read count at the Poisson p = 0.001 level "
            "expected from the surrounding ~20 bp: in background it removes nearly all duplicates, "
            "in strong peaks few.")
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
    print(f"{args.sample}: " + ", ".join(f"{k} {int(b):,} -> {int(a):,}" for k, (b, a) in counts.items()))


if __name__ == "__main__":
    main()
