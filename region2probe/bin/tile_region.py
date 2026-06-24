#!/usr/bin/env python3
"""
Created on Mon Oct 20 13:24:40 2025

@author: zshong
"""

import sys
import argparse

def select_probes_by_gap(input_file, output_file, target_gap):
    current_probe = None
    candidate = None
    best_error = float("inf")

    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            line = line.strip()
            if not line:
                continue

            # Fast string parsing: "GENE::chr:start-end"
            try:
                parts = line.split(":")
                chrom = parts[2]
                start_str, end_str = parts[3].split("-")
                start, end = int(start_str), int(end_str)
            except (IndexError, ValueError):
                continue

            # 1. Anchor: Always pick the very first valid probe in the file
            if current_probe is None:
                current_probe = (chrom, start, end, line)
                outfile.write(f"{chrom}\t{start}\t{end}\t{line}\n")
                continue

            # 2. Strict non-overlapping rule check
            # (The gap between them cannot be negative)
            actual_gap = start - current_probe[2]
            if actual_gap <= 0:
                continue

            # How close this actual gap is to the user's target gap
            current_error = abs(actual_gap - target_gap)

            # 3. If starting a brand new candidate evaluation for this window
            if candidate is None:
                candidate = (chrom, start, end, line)
                best_error = current_error
                continue

            # 4. Evaluate if this probe creates a gap closer to the target
            if current_error <= best_error:
                candidate = (chrom, start, end, line)
                best_error = current_error
            else:
                # If the gap error starts increasing, we've passed the sweet spot.
                # Commit the best candidate found for this window.
                outfile.write(f"{candidate[0]}\t{candidate[1]}\t{candidate[2]}\t{candidate[3]}\n")
                current_probe = candidate
                
                # Check if the current line can be reused for the next window
                next_actual_gap = start - current_probe[2]
                if next_actual_gap >= 0:
                    candidate = (chrom, start, end, line)
                    best_error = abs(next_actual_gap - target_gap)
                else:
                    candidate = None

        # Post-loop: Commit the final remaining candidate at End-of-File
        if candidate:
            outfile.write(f"{candidate[0]}\t{candidate[1]}\t{candidate[2]}\t{candidate[3]}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="select probes with given bp gap")
    parser.add_argument("--output", required=True, help="abs path to output bed file")
    parser.add_argument("--input", required=True, help="abs path to intput probe ids txt file")
    parser.add_argument("--bp", required=False, type=int, default=5, help="ideal bp gap between probes")
    args = parser.parse_args()
    
    if "--help" in sys.argv:
        parser.print_help()
        sys.exit()

    select_probes_by_gap(args.input, args.output, args.bp)


