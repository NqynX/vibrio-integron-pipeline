#!/usr/bin/env python3

"""
Plot integron cassette architecture from reconstructed cassette table.
"""


import argparse
import pandas as pd
import matplotlib.pyplot as plt


def parse_arguments():

    parser = argparse.ArgumentParser(
        description="Plot integron cassette array"
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Cassette array final TSV"
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output figure"
    )

    return parser.parse_args()



def plot_integron(df, output):

    fig, ax = plt.subplots(
        figsize=(14, 3)
    )

    # Integron boundaries
    integron_start = df["integron_start"].iloc[0]
    integron_end = df["integron_end"].iloc[0]

    integrase_start = df["integrase_start"].iloc[0]
    integrase_end = df["integrase_end"].iloc[0]


    # Draw integron backbone

    ax.plot(
        [integron_start, integron_end],
        [0, 0],
        linewidth=2
    )


    # Draw cassette blocks

    for _, row in df.iterrows():

        start = row["start"]
        length = row["cassette_length"]

        if row["protein_count"] > 0:
            height = 0.5
        else:
            height = 0.25


        ax.add_patch(
            plt.Rectangle(
                (start, -height/2),
                length,
                height
            )
        )


    # Draw integrase

    ax.add_patch(
        plt.Rectangle(
            (
                integrase_start,
                0.3
            ),
            integrase_end - integrase_start,
            0.3
        )
    )


    ax.text(
        integrase_start,
        0.8,
        "intI",
        fontsize=10
    )


    ax.set_xlim(
        integron_start - 5000,
        integron_end + 5000
    )


    ax.set_ylim(
        -1,
        1
    )


    ax.set_xlabel(
        "Chromosomal position (bp)"
    )

    ax.set_yticks([])

    ax.set_title(
        "Vibrio cholerae N16961 chromosomal integron cassette array"
    )


    plt.tight_layout()

    plt.savefig(
        output,
        dpi=300,
        bbox_inches="tight"
    )


def main():

    args = parse_arguments()

    df = pd.read_csv(
        args.input,
        sep="\t"
    )

    plot_integron(
        df,
        args.output
    )


if __name__ == "__main__":
    main()
