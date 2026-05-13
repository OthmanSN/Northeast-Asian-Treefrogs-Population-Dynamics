#!/usr/bin/env python3
"""Plot a Dsuite f-branch matrix as a heatmap.

Example:
    python scripts/03_plot_fbranch_heatmap.py \
        --input results/dsuite/fbranch_out.txt \
        --output results/dsuite/fbranch_heatmap.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot Dsuite f-branch heatmap.")
    parser.add_argument("--input", required=True, help="Path to Dsuite fbranch output file.")
    parser.add_argument("--output", default="fbranch_heatmap.png", help="Output image path.")
    parser.add_argument("--dpi", type=int, default=300, help="Output resolution.")
    return parser.parse_args()


def load_fbranch_matrix(path: Path) -> pd.DataFrame:
    fb = pd.read_csv(path, sep="\t")

    if "branch" not in fb.columns:
        raise ValueError("Expected a 'branch' column in the f-branch output.")

    matrix = (
        fb.drop(columns=["branch_descendants"], errors="ignore")
        .set_index("branch")
        .apply(pd.to_numeric, errors="coerce")
    )
    return matrix


def plot_heatmap(matrix: pd.DataFrame, output: Path, dpi: int) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(8, 6))
    ax = sns.heatmap(
        matrix,
        annot=True,
        fmt=".3f",
        cmap="Reds",
        linewidths=0.7,
        linecolor="grey",
        cbar_kws={"label": "f-branch value"},
        mask=matrix.isna(),
        annot_kws={"size": 9},
    )
    ax.set_title("f-branch introgression matrix", fontsize=14)
    ax.set_xlabel("Recipient taxon")
    ax.set_ylabel("Donor branch")
    plt.tight_layout()
    plt.savefig(output, dpi=dpi)
    plt.close()


def main() -> None:
    args = parse_args()
    matrix = load_fbranch_matrix(Path(args.input))
    plot_heatmap(matrix, Path(args.output), args.dpi)
    print(f"Saved heatmap to {args.output}")


if __name__ == "__main__":
    main()
