# -*- coding: utf-8 -*-
"""
Scatter plots of read depth vs VAF with BG error filtering

Example:
python plot_bg_error.py \
  -i noncoding_lowcomp_removed.tsv \
  -o noncoding_bg50 \
  --vaf-pctl 95 \
  --depth-pctl 50 \
  --bg-cutoff 50 \
  --save-fig /groups/wyattgrp/users/zshong/projects/tmp_outputs
"""

import argparse
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(
        description="Plot read depth vs VAF with BG error filtering"
    )

    p.add_argument(
        "-i", "--input",
        required=True,
        help="Input TSV file"
    )

    p.add_argument(
        "-o", "--tsv_out",
        required=False,
        default=None,
        help="Output prefix for filtered TSV"
    )

    p.add_argument(
        "--vaf-pctl",
        type=float,
        default=95,
        help="VAF percentile cutoff (default: 95)"
    )

    p.add_argument(
        "--depth-pctl",
        type=float,
        default=50,
        help="Read depth percentile cutoff (default: 50)"
    )

    p.add_argument(
    "--bg-cutoff",
    type=float,
    default=None,  # None means “don’t filter by BG error”
    help="Optional BG error cutoff (VAF / BG allele%%). If not set, BG error filter is skipped."
    )

    p.add_argument(
        "--max-vaf",
        type=float,
        default=100,
        help="Maximum VAF to keep (default: 100)"
    )

    p.add_argument(
    "--save-fig",
    type=str,
    default=None,
    help="Optional file path to save the figure (e.g., output.png or output.pdf). If not provided, figure will just be shown interactively."
    )

    return p.parse_args()


def main():
    args = parse_args()

    # Load data
    df = pd.read_csv(args.input, sep="\t")

    # Coerce numeric columns
    df["read_depth"] = pd.to_numeric(df["read_depth"], errors="coerce")
    df["VAF"] = pd.to_numeric(df["VAF"], errors="coerce")

    # Clean BG allele%
    df["BG allele%"] = (
        df["BG allele%"]
        .astype(str)
        .str.replace("%", "", regex=False)
        .replace("", np.nan)
        .astype(float)
    )

    # Basic filters
    df = df[df["VAF"] <= args.max_vaf]
    df = df[df["Chrom"] != "chrM"]
    df = df[df["Chrom"].str.contains("chr", na=False)]

    # BG error calculation
    df["bg error"] = df["VAF"] / df["BG allele%"]
    df.loc[df["BG allele%"] == 0, "bg error"] = np.nan

    if args.bg_cutoff is not None:
        df = df[
            (df["bg error"].isna()) |
            (df["bg error"] >= args.bg_cutoff)
        ]

    samples = sorted(df["Sample"].dropna().unique())

    # Plotting
    n_cols = 5
    n_rows = math.ceil(len(samples) / n_cols)

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(4 * n_cols, 3 * n_rows),
        sharex=True,
        sharey=True
    )
    axes = axes.flatten()
    red_points = []

    for ax, sample in zip(axes, samples):
        sub = df[df["Sample"] == sample]

        if sub.empty:
            ax.axis("off")
            continue

        depth_cut = sub["read_depth"].quantile(args.depth_pctl / 100)
        vaf_cut = sub["VAF"].quantile(args.vaf_pctl / 100)

        selected = (
            (sub["read_depth"] >= depth_cut) &
            (sub["VAF"] >= vaf_cut)
        )

        red_points.append(sub.loc[selected])

        ax.scatter(
            sub.loc[~selected, "read_depth"],
            sub.loc[~selected, "VAF"],
            s=3, alpha=0.3, color="blue"
        )

        ax.scatter(
            sub.loc[selected, "read_depth"],
            sub.loc[selected, "VAF"],
            s=3, alpha=0.7, color="red"
        )

        ax.axvline(depth_cut, color="blue", lw=1)
        ax.axhline(vaf_cut, color="blue", lw=1)

        ax.text(
            0.98, 0.95,
            f"depth: {depth_cut:.1f}\nVAF: {vaf_cut:.1f}",
            transform=ax.transAxes,
            ha="right", va="top",
            fontsize=8,
            color="blue",
            bbox=dict(facecolor="white", alpha=0.7, edgecolor="none")
        )

        ax.text(
            0.98, 0.25,
            f"red: {selected.sum()}\nblue: {(~selected).sum()}",
            transform=ax.transAxes,
            ha="right", va="bottom",
            fontsize=9
        )

        ax.set_title(sample, fontsize=10)
        ax.set_xlim(-200, 2000)

    # Remove unused panels
    for ax in axes[len(samples):]:
        ax.axis("off")

    fig.suptitle(
        f"VAF ≥ {args.vaf_pctl}%-tile | "
        f"Depth ≥ {args.depth_pctl}%-tile | "
        f"BG error ≥ {args.bg_cutoff}\n{args.input} | "
        f"\n {args.input}",
        fontsize=12
    )

    fig.text(0.5, 0.04, "Read depth", ha="center")
    fig.text(0.04, 0.5, "VAF%", va="center", rotation="vertical")

    plt.tight_layout(rect=[0.05, 0.05, 1, 0.93])

    if args.save_fig:
        plt.savefig(args.save_fig, dpi=300, bbox_inches='tight')
        print(f"Saved figure to: {args.save_fig}")
    else:
        plt.show()


    # Save filtered red points if tsv_out specified
    if args.tsv_out:
        red_df = pd.concat(red_points, ignore_index=True)
        out_tsv = f"{args.tsv_out}_selected.tsv"
        red_df.to_csv(out_tsv, sep="\t", index=False)
    

    print(f"Saved selected variants to: {out_tsv}")


if __name__ == "__main__":
    main()
