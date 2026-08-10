#!/usr/bin/env python3

import argparse
import pandas as pd


def parse_arguments():

    parser = argparse.ArgumentParser(
        description="Extract cassette arrays from IntegronFinder output"
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input IntegronFinder .integrons file"
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output cassette array TSV"
    )

    return parser.parse_args()



def extract_cassettes(df):

    cassette_records = []

    # Keep only relevant features
    df = df[df["type_elt"].isin(["attC", "protein"])]

    # Analyse each integron independently
    for integron_id, group in df.groupby("ID_integron"):

        group = group.sort_values(
            by="pos_beg"
        ).reset_index(drop=True)


        cassette_number = 1


        for i in range(len(group)):

            current = group.iloc[i]


            # Cassette boundaries start with attC
            if current["type_elt"] != "attC":
                continue


            # Find next attC site
            next_attc_index = None

            for j in range(i + 1, len(group)):

                if group.iloc[j]["type_elt"] == "attC":
                    next_attc_index = j
                    break


            if next_attc_index is None:
                continue


            next_attc = group.iloc[next_attc_index]


            # Features between attC boundaries
            internal_features = group.iloc[
                i + 1:next_attc_index
            ]


            proteins = internal_features[
                internal_features["type_elt"] == "protein"
            ]


            cassette_records.append(
                {
                    "integron_id": integron_id,
                    "replicon": current["ID_replicon"],
                    "cassette_number": cassette_number,
                    "start": int(current["pos_end"]),
                    "end": int(next_attc["pos_beg"]),
                    "left_attC": current["element"],
                    "right_attC": next_attc["element"],
                    "protein_count": len(proteins)
                }
            )


            cassette_number += 1


    return pd.DataFrame(cassette_records)



def main():

    args = parse_arguments()


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
