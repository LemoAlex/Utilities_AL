#!/usr/bin/env python3
"""
merge_blast_hsps.py

Merges HSPs from a BLAST tabular output (outfmt 7) file by subject sequence,
combining overlapping or nearby regions (within a user-defined range).
Coordinates are normalized from smaller to larger.

Usage:
    python merge_blast_hsps.py blast_output.txt --range 100 > merged_output.tsv
"""

import sys
import argparse
from collections import defaultdict

def parse_blast_outfmt7(file_path):
    """Parse BLAST outfmt 7 file, ignoring comment lines."""
    records = []
    with open(file_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 12:
                continue
            qid, sid = parts[0], parts[1]
            sstart, send = int(parts[8]), int(parts[9])
            start, end = sorted((sstart, send))
            records.append((qid, sid, start, end))
    return records

def merge_intervals(intervals, merge_range):
    """Merge intervals that overlap or are within merge_range bp."""
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda x: x[0])
    merged = [intervals[0]]
    for current in intervals[1:]:
        prev = merged[-1]
        # Check if current start is close enough to previous end
        if current[0] <= prev[1] + merge_range:
            merged[-1] = (prev[0], max(prev[1], current[1]))
        else:
            merged.append(current)
    return merged

def main(blast_file, merge_range):
    # Group by query-subject pair
    hits = defaultdict(list)
    for qid, sid, start, end in parse_blast_outfmt7(blast_file):
        hits[(qid, sid)].append((start, end))

    # Merge and print
    print("query\tsubject\tmerged_start\tmerged_end")
    for (qid, sid), intervals in hits.items():
        merged = merge_intervals(intervals, merge_range)
        for start, end in merged:
            print(f"{qid}\t{sid}\t{start}\t{end}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge BLAST HSPs within a given range.")
    parser.add_argument("blast_output", help="BLAST outfmt 7 file")
    parser.add_argument("--range", "-r", type=int, default=0,
                        help="Maximum distance (bp) between HSPs to merge (default: 0 = only overlapping/adjacent)")
    args = parser.parse_args()

    main(args.blast_output, args.range)

