#!/usr/bin/env python3
"""Plot a zoomed-in linear map of an integron cassette array.

Shows: integrase (intI), attC sites (triangles), gene cassettes (arrows),
with genomic coordinates and scale bar.
"""

import argparse
import csv
import sys

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch


def parse_integrons(filepath):
    """Parse IntegronFinder .integrons file into structured data."""
    integrases = []
    attc_sites = []
    proteins = []
    with open(filepath) as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if not row[0].startswith('integron'):
                continue
            if len(row) < 8:
                continue
            int_id = row[0]
            pos_beg = int(row[3])
            pos_end = int(row[4])
            strand = int(row[5])
            evalue = row[6]
            elt_type = row[7]
            annotation = row[8]
            
            if 'intI' in annotation:
                integrases.append({'id': int_id, 'start': pos_beg, 'end': pos_end, 'strand': strand})
            elif row[2].startswith('attc_'):
                attc_sites.append({'id': int_id, 'name': row[2], 'start': pos_beg, 'end': pos_end,
                                   'evalue': float(evalue) if evalue != 'NA' else None})
            elif elt_type == 'protein':
                proteins.append({'id': int_id, 'start': pos_beg, 'end': pos_end, 'strand': strand})

    return integrases, attc_sites, proteins


def plot_zoomed(integrons_file, region_start, region_end, output, title=None):
    """Draw a zoomed linear map of one integron region."""
    integrases, attc_sites, proteins = parse_integrons(integrons_file)

    # Filter to one integron (pick the one overlapping the region)
    int_ids = set()
    for a in attc_sites:
        if region_start <= a['start'] <= region_end:
            int_ids.add(a['id'])
    if not int_ids:
        print(f"No attC sites found in region {region_start}-{region_end}")
        sys.exit(1)
    
    target_id = sorted(int_ids)[0]
    attcs = sorted([a for a in attc_sites if a['id'] == target_id], key=lambda x: x['start'])
    intgs = [i for i in integrases if i['id'] == target_id]
    prots = [p for p in proteins if p['id'] == target_id]
    
    # Determine plot bounds
    all_starts = [a['start'] for a in attcs] + [p['start'] for p in prots]
    all_ends = [a['end'] for a in attcs] + [p['end'] for p in prots]
    if intgs:
        all_starts.append(intgs[0]['start'])
        all_ends.append(intgs[0]['end'])
    
    x_min = min(all_starts) - 500
    x_max = max(all_ends) + 500
    region_len = x_max - x_min

    # Setup figure
    fig, ax = plt.subplots(figsize=(16, 4))
    fig.subplots_adjust(left=0.05, right=0.95, top=0.85, bottom=0.25)

    # Scale function: genomic coord -> plot x
    margin = 0.05
    def gx(pos):
        return margin + (pos - x_min) / region_len * (1 - 2 * margin)

    # Draw backbone line
    ax.plot([0, 1], [0.5, 0.5], color='#333333', linewidth=2, zorder=1)

    # Draw integrase
    if intgs:
        intl = intgs[0]
        x1, x2 = gx(intl['start']), gx(intl['end'])
        w = x2 - x1
        color = '#E74C3C'
        rect = mpatches.FancyBboxPatch((x1, 0.55), w, 0.2, boxstyle="round,pad=0.01",
                                        facecolor=color, edgecolor='#333', linewidth=1, zorder=3)
        ax.add_patch(rect)
        ax.text((x1 + x2) / 2, 0.65, 'intI', ha='center', va='center',
                fontsize=10, fontweight='bold', color='white', zorder=4)

    # Draw protein genes (between attC sites)
    for p in prots:
        x1, x2 = gx(p['start']), gx(p['end'])
        w = x2 - x1
        if w < 0.001:
            continue
        color = '#3498DB'
        y_base = 0.55
        # Arrow direction based on strand
        if p['strand'] == 1:
            arrow = mpatches.FancyArrow(x1, y_base + 0.1, w, 0, width=0.15,
                                         head_width=0.2, head_length=min(w * 0.15, 0.015),
                                         facecolor=color, edgecolor='#333', linewidth=0.5, zorder=3)
        else:
            arrow = mpatches.FancyArrow(x2, y_base + 0.1, -w, 0, width=0.15,
                                         head_width=0.2, head_length=min(w * 0.15, 0.015),
                                         facecolor=color, edgecolor='#333', linewidth=0.5, zorder=3)
        ax.add_patch(arrow)

    # Draw attC sites
    for i, a in enumerate(attcs):
        xc = gx((a['start'] + a['end']) / 2)
        # Triangle pointing down
        tri_x = [xc - 0.006, xc + 0.006, xc]
        tri_y = [0.5, 0.5, 0.35]
        ax.fill(tri_x, tri_y, color='#2ECC71', edgecolor='#333', linewidth=0.5, zorder=3)
        # Label
        ax.text(xc, 0.28, f"attC\n{i+1}", ha='center', va='top', fontsize=6, color='#2ECC71')

    # Scale bar
    scale_bp = round(region_len / 5 / 100) * 100  # ~5 divisions, round to 100bp
    if scale_bp < 100:
        scale_bp = 100
    scale_x1 = 0.7
    scale_x2 = scale_x1 + (scale_bp / region_len) * (1 - 2 * margin)
    ax.plot([scale_x1, scale_x2], [0.12, 0.12], color='black', linewidth=2)
    ax.plot([scale_x1, scale_x1], [0.10, 0.14], color='black', linewidth=1.5)
    ax.plot([scale_x2, scale_x2], [0.10, 0.14], color='black', linewidth=1.5)
    ax.text((scale_x1 + scale_x2) / 2, 0.07, f'{scale_bp:,} bp', ha='center', va='top', fontsize=9)

    # Coordinate labels
    ax.text(0.0, 0.12, f'{x_min:,}', ha='left', va='top', fontsize=7, color='#666')
    ax.text(1.0, 0.12, f'{x_max:,}', ha='right', va='top', fontsize=7, color='#666')

    # Title
    int_type = "complete integron" if intgs else "CALIN"
    if title is None:
        title = f"V. cholerae N16961 chrII — {int_type} ({target_id})"
    ax.set_title(title, fontsize=13, fontweight='bold', pad=10)

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor='#E74C3C', edgecolor='#333', label='Integrase (intI)'),
        mpatches.Patch(facecolor='#3498DB', edgecolor='#333', label='Gene cassette ORF'),
        mpatches.Patch(facecolor='#2ECC71', edgecolor='#333', label='attC site'),
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=8, framealpha=0.9)

    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(0.0, 0.95)
    ax.axis('off')

    plt.savefig(output, dpi=300, bbox_inches='tight')
    print(f"Saved: {output}")


def main():
    parser = argparse.ArgumentParser(description="Plot zoomed integron map")
    parser.add_argument("-i", "--integrons", required=True, help="IntegronFinder .integrons file")
    parser.add_argument("-s", "--start", type=int, required=True, help="Region start (bp)")
    parser.add_argument("-e", "--end", type=int, required=True, help="Region end (bp)")
    parser.add_argument("-o", "--output", required=True, help="Output PNG path")
    parser.add_argument("-t", "--title", default=None, help="Custom title")
    args = parser.parse_args()
    
    plot_zoomed(args.integrons, args.start, args.end, args.output, args.title)


if __name__ == '__main__':
    main()
