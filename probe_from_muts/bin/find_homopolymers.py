#!/usr/bin/env python

from Bio import SeqIO
import re
import sys
import argparse


def find_homopolymers(input, output, bp):
    homopolymer_pattern = re.compile(
        rf'(A{{{bp+1},}}|T{{{bp+1},}}|G{{{bp+1},}}|C{{{bp+1},}}|N{{{bp+1},}})',
        re.IGNORECASE)

    with open(output, "w") as out_handle:
        for record in SeqIO.parse(input, "fasta"):
            seq_str = str(record.seq)
            if homopolymer_pattern.search(seq_str):
                out_handle.write(record.id + "\n")

    print(f"homopolymers > {bp} bp written to:", output)


def main():
    parser = argparse.ArgumentParser(description="filter blat psl output file for ideal probe seqs")
    parser.add_argument("--path_output", required=True, help="abs path to output txt file")
    parser.add_argument("--path_input", required=True, help="abs path to intput fa file")
    parser.add_argument("--bp", required=False, type=int, default=5, help="threshold for bp in homopolymer")
    args = parser.parse_args()
    
    if "--help" in sys.argv:
        parser.print_help()
        sys.exit()

    find_homopolymers(args.path_input, args.path_output, args.bp)
    
if __name__ == "__main__":
    main()