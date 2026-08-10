#!/usr/bin/env python3

"""
Plot wrapped integron cassette architecture.

Creates a readable schematic for large cassette arrays.
"""


import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle



def parse_arguments():

    parser = argparse.ArgumentParser(
        description="Plot integron cassette array"
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True
    )

    return parser.parse_args()



def plot_integron(df, output):

    fig, ax = plt.subplots(
        figsize=(14, 6)
    )


    integron_start = df["integron_start"].iloc[0]
    integron_end = df["integron_end"].iloc[0]


    int_start = df["integrase_start"].iloc[0]


    # layout parameters
    per_row = 50
    box_spacing = 1.2
    row_spacing = 1.5


    max_rows = (
        len(df) + per_row - 1
    ) // per_row


    for i, (_, row) in enumerate(df.iterrows()):

        row_number = i // per_row
        column_number = i % per_row


        x = column_number * box_spacing
        y = -row_number * row_spacing


        # scale height only
        if row["protein_count"] > 0:
            height = 0.5
        else:
            height = 0.25


        ax.add_patch(
            Rectangle(
                (x, y),
                1,
                height
            )
        )


        # label every 10 cassettes

        if row["cassette_number"] % 10 == 0:

            ax.text(
                x,
                y + 0.65,
                f"C{row['cassette_number']}",
                fontsize=7,
                ha="center"
            )


    # integrase annotation

    ax.text(
        0,
        1,
        "intI",
        fontsize=12,
        ha="left"
    )


    ax.text(
        0,
        0.75,
        f"Integron: {integron_start:,}-{integron_end:,} bp",
        fontsize=9
    )


    ax.set_title(
        "Vibrio cholerae N16961 chromosomal integron cassette array",
        fontsize=14
    )


    ax.set_xlim(
        -2,
        per_row * box_spacing
    )


    ax.set_ylim(
        -max_rows * row_spacing,
        2
    )


    ax.axis("off")


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
