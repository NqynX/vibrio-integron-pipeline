#!/usr/bin/env python3

import argparse
import pandas as pd


def parse_arguments():

    parser = argparse.ArgumentParser(
        description="Reconstruct cassette arrays from IntegronFinder output"
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Path to IntegronFinder .integrons file"
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output TSV file"
    )

    return parser.parse_args()


def extract_cassettes(df):

    cassette_records = []

    # Keep only biological elements
    df = df[df["type_elt"].isin(["attC", "protein"])]

    # Process each detected integron separately
    for integron_id, group in df.groupby("ID_integron"):

        # Ensure genomic order
        group = group.sort_values(
            by="pos_beg"
        ).reset_index(drop=True)

        cassette_number = 1

        # Scan through all features
        for i in range(len(group)):

            current = group.iloc[i]

            # Cassette starts after an attC site
            if current["type_elt"] != "attC":
                continue


            # Find the next attC site
            next_attc_index = None

            for j in range(i + 1, len(group)):

                if group.iloc[j]["type_elt"] == "attC":
                    next_attc_index = j
                    break


            # No downstream attC
            if next_attc_index is None:
                continue


            next_attc = group.iloc[next_attc_index]


            # Features between the two attC sites
            cassette_features = group.iloc[
                i + 1:next_attc_index
            ]


            # Keep only predicted proteins
            proteins = cassette_features[
                cassette_features["type_elt"] == "protein"
            ]


            # Ignore empty intervals
            if proteins.empty:
                continue


            cassette_records.append(
                {
                    "integron_id": integron_id,
                    "replicon": current["ID_replicon"],
                    "cassette_number": cassette_number,
                    "start": int(proteins["pos_beg"].min()),
                    "end": int(proteins["pos_end"].max()),
                    "left_attC": current["element"],
                    "right_attC": next_attc["element"],
                    "protein_count": len(proteins)
                }
            )

            cassette_number += 1


    return pd.DataFrame(cassette_records)


def main():

    args = parse_arguments()


    # Read IntegronFinder table
    df = pd.read_csv(
        args.input,
        sep="\t",
        comment="#"
    )


    cassette_df = extract_cassettes(df)


    cassette_df.to_csv(
        args.output,
        sep="\t",
        index=False
    )


    print(
        f"Extracted {len(cassette_df)} predicted cassettes"
    )


if __name__ == "__main__":
    main()
