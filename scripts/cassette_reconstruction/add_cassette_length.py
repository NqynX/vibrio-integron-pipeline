#!/usr/bin/env python3

import argparse
import pandas as pd


def main():

    parser = argparse.ArgumentParser(
        description="Add cassette length column"
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

    args = parser.parse_args()


    df = pd.read_csv(
        args.input,
        sep="\t"
    )


    df["cassette_length"] = (
        df["end"] - df["start"]
    )


    df.to_csv(
        args.output,
        sep="\t",
        index=False
    )


    print(
        f"Added cassette lengths for {len(df)} cassettes"
    )


if __name__ == "__main__":
    main()
