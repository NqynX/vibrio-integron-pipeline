#!/usr/bin/env python3
"""Build a cassette table TSV from IntegronFinder .integrons output.

Cassette boundaries = regions between consecutive attC sites
within each integron/CALIN cluster.
"""

import csv
import sys

def main():
    if len(sys.argv) < 3:
        print("Usage: python build_cassette_table.py <integrons_file> <output_tsv>")
        sys.exit(1)

    infile = sys.argv[1]
    outfile = sys.argv[2]

    # Parse: group attC sites by integron ID
    integron_attcs = {}
    with open(infile) as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if not row[0].startswith('integron'):
                continue
            if len(row) < 8:
                continue
            if not row[2].startswith('attc_'):
                continue
            integron_id = row[0]
            replicon = row[1]
            pos_beg = int(row[3])
            pos_end = int(row[4])
            if integron_id not in integron_attcs:
                integron_attcs[integron_id] = {'replicon': replicon, 'sites': []}
            integron_attcs[integron_id]['sites'].append((pos_beg, pos_end, row[2]))

    # Sort attC sites by position within each integron, build cassettes
    cassettes = []
    cassette_num = 1

    for int_id in sorted(integron_attcs.keys()):
        info = integron_attcs[int_id]
        sites = sorted(info['sites'], key=lambda x: x[0])
        replicon = info['replicon']

        print(f"{int_id}: {len(sites)} attC sites on {replicon}")

        # Cassette = region from end of attC_n to start of attC_{n+1}
        # This captures the gene between two attC sites
        for i in range(len(sites) - 1):
            start = sites[i][1] + 1    # 1 bp after end of current attC
            end = sites[i + 1][0] - 1  # 1 bp before start of next attC
            length = end - start + 1
            cassettes.append({
                'replicon': replicon,
                'cassette_number': cassette_num,
                'start': start,
                'end': end,
                'length': length,
                'integron_id': int_id,
                'upstream_attc': sites[i][2],
                'downstream_attc': sites[i + 1][2],
            })
            cassette_num += 1

    # Write TSV
    cols = ['replicon', 'cassette_number', 'start', 'end', 'length', 'integron_id', 'upstream_attc', 'downstream_attc']
    with open(outfile, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for c in cassettes:
            f.write('\t'.join(str(c[col]) for col in cols) + '\n')

    print(f"\nTotal cassettes: {len(cassettes)}")
    print(f"Output: {outfile}")

if __name__ == '__main__':
    main()
