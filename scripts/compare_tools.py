#!/usr/bin/env python3
"""Compare attC detections from IntegronFinder and HattCI."""

import sys
import csv

def parse_integronfinder(filepath):
    """Parse IntegronFinder .integrons file, return list of (seqid, start, end, integron_id)."""
    sites = []
    with open(filepath) as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if row[0].startswith('#') or not row[0].startswith('integron'):
                continue
            if len(row) < 8:
                continue
            seqid = row[1]
            element = row[2]
            pos_beg = int(row[3])
            pos_end = int(row[4])
            if element.startswith('attc_'):
                sites.append((seqid, pos_beg, pos_end, row[0]))
    return sites

def parse_hattci(filepath, min_vscore=0):
    """Parse HattCI .out file, return list of (seqid, start, end, vscore)."""
    sites = []
    with open(filepath) as f:
        for i, line in enumerate(f):
            if i < 2:  # skip blank line + header
                continue
            if line.startswith('-'):
                break  # separator line
            parts = line.split()
            if len(parts) < 5:
                continue
            try:
                vscore = float(parts[4])
                start = int(parts[2])
                end = int(parts[3])
                seqid = parts[1]
                if vscore > min_vscore:
                    sites.append((seqid, start, end, vscore))
            except (ValueError, IndexError):
                continue
    return sites

def write_bed(sites, filepath, name_prefix):
    """Write sites to BED format (0-based, half-open)."""
    with open(filepath, 'w') as f:
        for site in sites:
            seqid, start, end = site[0], site[1], site[2]
            # Convert from 1-based closed to 0-based half-open for BED
            bed_start = start - 1
            label = f"{name_prefix}_{start}-{end}"
            f.write(f"{seqid}\t{bed_start}\t{end}\t{label}\n")

def compute_overlap(sites_a, sites_b, max_dist=50):
    """Find matching sites between two lists using proximity matching.
    
    A match is counted if start positions are within max_dist bp of each other.
    Each site can only match once (greedy, sorted by distance).
    Returns (matched_a, matched_b, matches_list).
    """
    matches = []
    used_a = set()
    used_b = set()
    
    # Collect all candidate pairs
    candidates = []
    for i, sa in enumerate(sites_a):
        for j, sb in enumerate(sites_b):
            dist = abs(sa[1] - sb[1])  # compare start positions
            if dist <= max_dist:
                candidates.append((dist, i, j))
    
    # Sort by distance, greedily assign matches
    candidates.sort()
    for dist, i, j in candidates:
        if i not in used_a and j not in used_b:
            matches.append((sites_a[i], sites_b[j], dist))
            used_a.add(i)
            used_b.add(j)
    
    matched_a = len(used_a)
    matched_b = len(used_b)
    return matched_a, matched_b, matches

def main():
    if len(sys.argv) < 3:
        print("Usage: python compare_tools.py <integronfinder.integrons> <hattci.out> [min_vscore]")
        sys.exit(1)
    
    if_file = sys.argv[1]
    hattci_file = sys.argv[2]
    min_vscore = float(sys.argv[3]) if len(sys.argv) > 3 else 0
    
    # Parse both tools
    if_sites = parse_integronfinder(if_file)
    hattci_sites = parse_hattci(hattci_file, min_vscore)
    
    print(f"IntegronFinder attC sites: {len(if_sites)}")
    print(f"HattCI hits (Vscore > {min_vscore}): {len(hattci_sites)}")
    print()
    
    # Write BED files
    write_bed(if_sites, "integronfinder_attc.bed", "IF")
    write_bed(hattci_sites, "hattci_attc.bed", "HC")
    print("Wrote: integronfinder_attc.bed")
    print("Wrote: hattci_attc.bed")
    print()
    
    # Compute overlap at different distance thresholds
    for max_dist in [10, 25, 50, 100]:
        ma, mb, matches = compute_overlap(if_sites, hattci_sites, max_dist)
        tp = ma  # true positives = IF sites that matched a HattCI hit
        
        precision = tp / len(hattci_sites) if len(hattci_sites) > 0 else 0
        recall = tp / len(if_sites) if len(if_sites) > 0 else 0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
        
        print(f"--- Proximity threshold: {max_dist} bp ---")
        print(f"  IF sites matched: {ma}/{len(if_sites)} (recall={recall:.3f})")
        print(f"  HattCI hits matched: {mb}/{len(hattci_sites)} (precision={precision:.3f})")
        print(f"  F1 score: {f1:.3f}")
        print(f"  HattCI-only hits: {len(hattci_sites) - mb}")
        print(f"  IF-only sites: {len(if_sites) - ma}")
        print()
    
    # Show top matches with distances
    _, _, matches = compute_overlap(if_sites, hattci_sites, 50)
    matches.sort(key=lambda x: x[2])
    print("=== Closest matches (top 10) ===")
    for sa, sb, dist in matches[:10]:
        print(f"  IF {sa[1]}-{sa[2]}  vs  HattCI {sb[1]}-{sb[2]}  Vscore={sb[3]:.2f}  dist={dist}bp")

if __name__ == '__main__':
    main()
