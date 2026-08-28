#!/usr/bin/env python3
"""
Extract cassette nucleotide sequences from a genome FASTA.

Input:  cassette TSV with columns replicon, cassette_number, start, end
        plus the genome FASTA those coordinates refer to.
Output: one combined multi-FASTA, optionally one FASTA per cassette.

Coordinates are treated as 1-based and inclusive, so extracted length is
(end - start + 1). The upstream cassette_length column used (end - start),
so extracted lengths may be 1 bp longer. Mismatches are reported, not hidden.
"""

import argparse
import os
import sys

import pandas as pd
from Bio import SeqIO


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Extract cassette nucleotide sequences from a genome FASTA."
    )
    parser.add_argument("--cassettes", required=True, help="Cassette table (TSV).")
    parser.add_argument("--genome", required=True, help="Genome FASTA file.")
    parser.add_argument("--combined", required=True, help="Combined multi-FASTA output.")
    parser.add_argument("--outdir", default=None, help="Directory for per-cassette FASTA files.")
    parser.add_argument("--sample-name", default="sample", help="Label used in FASTA headers.")
    parser.add_argument("--flank", type=int, default=0, help="Extra bases each side.")
    parser.add_argument("--split", action="store_true", help="Write one FASTA per cassette.")
    parser.add_argument("--force", action="store_true", help="Allow overwriting output.")
    return parser.parse_args()


def load_genome_sequences(genome_path):
    """Return a dictionary mapping sequence ID to sequence."""
    sequences = {}
    for record in SeqIO.parse(genome_path, "fasta"):
        sequences[record.id] = record.seq
    return sequences


def extract_one_cassette(sequence, row, flank):
    """Return the nucleotide sequence for a single cassette row."""
    start_index = int(row["start"]) - 1 - flank
    end_index = int(row["end"]) + flank

    start_index = max(0, start_index)
    end_index = min(len(sequence), end_index)

    return sequence[start_index:end_index]


def build_fasta_header(sample_name, row, extracted_length, flank):
    """Build a descriptive FASTA header line for one cassette."""
    cassette_label = "cassette_{:03d}".format(int(row["cassette_number"]))
    sequence_name = "{}|{}|{}".format(sample_name, row["replicon"], cassette_label)

    fields = [
        "start={}".format(int(row["start"])),
        "end={}".format(int(row["end"])),
        "extracted_length={}".format(extracted_length),
        "flank={}".format(flank),
    ]
    if "protein_count" in row:
        fields.append("proteins={}".format(int(row["protein_count"])))
    if "integron_id" in row:
        fields.append("integron={}".format(row["integron_id"]))

    return ">{} {}".format(sequence_name, " ".join(fields))


def main():
    args = parse_arguments()

    if os.path.exists(args.combined) and not args.force:
        print("ERROR: {} exists. Use --force to overwrite.".format(args.combined),
              file=sys.stderr)
        sys.exit(1)

    genome_sequences = load_genome_sequences(args.genome)
    print("Replicons in genome: {}".format(list(genome_sequences.keys())))

    cassette_table = pd.read_csv(args.cassettes, sep="\t")
    print("Cassette rows read: {}".format(len(cassette_table)))

    if args.split:
        if args.outdir is None:
            print("ERROR: --split requires --outdir.", file=sys.stderr)
            sys.exit(1)
        os.makedirs(args.outdir, exist_ok=True)

    written_count = 0
    missing_replicons = []
    length_mismatches = 0

    with open(args.combined, "w") as combined_handle:
        for _, row in cassette_table.iterrows():
            replicon_id = str(row["replicon"])

            if replicon_id not in genome_sequences:
                missing_replicons.append(replicon_id)
                continue

            cassette_sequence = extract_one_cassette(
                genome_sequences[replicon_id], row, args.flank
            )
            extracted_length = len(cassette_sequence)

            if "cassette_length" in row and args.flank == 0:
                if int(row["cassette_length"]) != extracted_length:
                    length_mismatches += 1

            header_line = build_fasta_header(
                args.sample_name, row, extracted_length, args.flank
            )

            combined_handle.write(header_line + "\n")
            combined_handle.write(str(cassette_sequence) + "\n")

            if args.split:
                filename = "cassette_{:03d}.fasta".format(int(row["cassette_number"]))
                path = os.path.join(args.outdir, filename)
                with open(path, "w") as single_handle:
                    single_handle.write(header_line + "\n")
                    single_handle.write(str(cassette_sequence) + "\n")

            written_count += 1

    print("Sequences written: {}".format(written_count))
    print("Combined FASTA: {}".format(args.combined))

    if missing_replicons:
        print("WARNING: replicons not in genome: {}".format(sorted(set(missing_replicons))))

    if length_mismatches:
        print("NOTE: {} rows where extracted length != cassette_length column "
              "(expected, end-start vs end-start+1).".format(length_mismatches))


if __name__ == "__main__":
    main()

