#!/usr/bin/env python3

"""
Output invalid probes from BLAT PSL output according to these rules:

For each probe sequence:
output probes that violate the following criterias:
    - ≤3 contiguous perfect matches (≥20 bp each)
    - no block >50 bp with >90% identity
* takes into account original intended sequence for each probe
    
Usage:
.../filter_psl.py input.psl output.psl
"""
    
import os
import sys
import argparse

# === FUNCTIONS ===
def parse_psl_line(line):
    cols = line.strip().split('\t')
    if len(cols) < 21:
        return None
    return {
        'match': int(cols[0]),
        'misMatch': int(cols[1]),
        'repMatch': int(cols[2]),
        'qNumInsert': int(cols[4]),
        'qBaseInsert': int(cols[5]),
        'tNumInsert': int(cols[6]),
        'tBaseInsert': int(cols[7]),
        'strand': cols[8],
        'qName': cols[9],
        'qSize': int(cols[10]),
        'qStart': int(cols[11]),
        'qEnd': int(cols[12]),
        'tName': cols[13],
        'tSize': int(cols[14]),
        'tStart': int(cols[15]),
        'tEnd': int(cols[16]),
        'blockCount': int(cols[17]),
        'blockSizes': [int(x) for x in cols[18].strip(',').split(',') if x]
    }

def is_offtarget(hit):
    if hit['match'] == hit['qSize'] and hit['misMatch'] == 0:
        return False
    align_len = hit['match'] + hit['misMatch'] + hit['repMatch'] \
                + hit['qBaseInsert'] + hit['tBaseInsert']
    identity = hit['match'] / align_len if align_len else 0
    long_blocks = sum(1 for b in hit['blockSizes'] if b >= 20)
    if long_blocks > 3:
        return True
    if align_len >= 50 and identity >= 0.90:
        return True
    return False

def filter_psl(psl_file, outfile):
    with open(psl_file) as f:
        lines = [l for l in f if l.strip() and not l.startswith('psLayout') and not l.startswith('match')]
    flagged = []
    for line in lines:
        hit = parse_psl_line(line)
        if not hit:
            continue
        if is_offtarget(hit):
            flagged.append(hit)
    if outfile:
        with open(outfile, 'w') as out:
            out.write("qName\ttName\ttStart\ttEnd\tblockSizes\tlongBlocks\tmatch\tmisMatch\n")
            for h in flagged:
                long_blocks = sum(1 for b in h['blockSizes'] if b >= 20)
                out.write(f"{h['qName']}\t{h['tName']}\t{h['tStart']}\t{h['tEnd']}\t"
                          f"{','.join(map(str,h['blockSizes']))}\t{long_blocks}\t"
                          f"{h['match']}\t{h['misMatch']}\n")
        
def main():
    parser = argparse.ArgumentParser(description="filter blat psl output file for ideal probe seqs")
    parser.add_argument("--path_output", required=True, help="abs path to output file")
    parser.add_argument("--path_input", required=True, help="abs path to intput file")
    args = parser.parse_args()
    
    if "--help" in sys.argv:
        parser.print_help()
        sys.exit()

    filter_psl(args.path_input, args.path_output)
    
if __name__ == "__main__":
    main()