#!/usr/bin/env python3

import argparse
import pandas as pd



def parse_arguments():

    parser = argparse.ArgumentParser(
        description="Add integron metadata to cassette arrays"
    )

    parser.add_argument(
        "-i",
        "--integrons",
        required=True,
        help="IntegronFinder .integrons file"
    )

    parser.add_argument(
        "-c",
        "--cassette",
        required=True,
        help="Cassette array TSV"
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output TSV"
    )

    return parser.parse_args()



def extract_integron_metadata(df):

    metadata = []

    # Only keep actual integrase predictions
    integrases = df[
        (df["type_elt"] == "protein") &
        (df["annotation"] == "intI")
    ]


    for integron_id, group in df.groupby("ID_integron"):

        integron_features = group[
            group["type_elt"].isin(
                ["attC", "protein"]
            )
        ]


        record = {
            "integron_id": integron_id,
            "replicon": group["ID_replicon"].iloc[0],
            "integron_start": int(
                integron_features["pos_beg"].min()
            ),
            "integron_end": int(
                integron_features["pos_end"].max()
            )
        }


        # Add integrase information if available
        intI = integrases[
            integrases["ID_integron"] == integron_id
        ]


        if not intI.empty:

            record["integrase_start"] = int(
                intI["pos_beg"].iloc[0]
            )

            record["integrase_end"] = int(
                intI["pos_end"].iloc[0]
            )

            record["integrase_strand"] = int(
                intI["strand"].iloc[0]
            )

        else:

            record["integrase_start"] = None
            record["integrase_end"] = None
            record["integrase_strand"] = None


        metadata.append(record)


    return pd.DataFrame(metadata)



def main():

    args = parse_arguments()


    integron_df = pd.read_csv(
        args.integrons,
        sep="\t",
        comment="#"
    )


    cassette_df = pd.read_csv(
        args.cassette,
        sep="\t"
    )


    metadata_df = extract_integron_metadata(
        integron_df
    )


    merged = cassette_df.merge(
        metadata_df,
        on=[
            "integron_id",
            "replicon"
        ],
        how="left"
    )


    merged.to_csv(
        args.output,
        sep="\t",
        index=False
    )


    print(
        f"Added metadata to {len(merged)} cassette entries"
    )



if __name__ == "__main__":
    main()
