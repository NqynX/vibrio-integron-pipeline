#!/usr/bin/env python3
"""Permutation distance vs genome-wide ANI: is order conserved beyond clonality?"""
import csv, os, random, statistics
from collections import defaultdict

PROJ = "/scratch/user/uqcngu19/vibrio-integron-pipeline"
A = f"{PROJ}/analysis"
random.seed(0)

def gene_index(pid):
    t = pid.split("|")[-1].split("_")[-1]
    return int(t) if t.isdigit() else -1

fam = {}
for line in open(f"{A}/clust_cluster.tsv"):
    rep, mem = line.rstrip("\n").split("\t")[:2]
    fam[mem] = rep
sp = {r["accession"]: r["species"] for r in
      csv.DictReader(open(f"{PROJ}/results/per_genome_with_species.csv"))}

arrays, ipos = defaultdict(list), {}
for r in csv.DictReader(open(f"{A}/protein_meta.csv")):
    g = gene_index(r["protein_id"])
    if g < 0: continue
    if r["is_integrase"].strip().lower() == "true":
        ipos[r["accession"]] = g; continue
    if r["protein_id"] in fam:
        arrays[r["accession"]].append((g, fam[r["protein_id"]]))

ordered, counts = {}, {}
for acc, it in arrays.items():
    it.sort()
    if acc in ipos and it and ipos[acc] > it[-1][0]: it = it[::-1]
    v = [f for _, f in it]
    if len(v) < 20: continue
    ordered[acc] = v
    c = defaultdict(int)
    for f in v: c[f] += 1
    counts[acc] = c
print(f"{len(ordered)} genomes with >=20 ordered cassettes")

def mtf(perm):
    pos = {v: i for i, v in enumerate(perm)}
    m = len(perm)
    while m > 1 and pos[m-2] < pos[m-1]: m -= 1
    return m - 1

_rc = {}
def rand_exp(n):
    if n not in _rc:
        _rc[n] = statistics.mean(mtf(random.sample(range(n), n)) for _ in range(100))
    return _rc[n]

# mash distances
md = {}
for line in open(f"{A}/mash/dists.tsv"):
    p = line.split("\t")
    if len(p) < 3: continue
    x = os.path.basename(p[0]).replace(".fna","")
    y = os.path.basename(p[1]).replace(".fna","")
    if x == y or x not in ordered or y not in ordered: continue
    md[(x, y) if x < y else (y, x)] = float(p[2])
print(f"{len(md)} genome pairs with mash distance\n")

BINS = [(0.000,0.0005,"ANI >99.95 (near-clonal)"),
        (0.0005,0.002,"99.8-99.95"),
        (0.002,0.005,"99.5-99.8"),
        (0.005,0.01, "99.0-99.5"),
        (0.01, 0.02, "98-99"),
        (0.02, 0.05, "95-98"),
        (0.05, 1.0,  "<95 (different species)")]
res = {b[2]: dict(n=0, shared=[], obs=[], exp=[], zero=0) for b in BINS}

for (x, y), d in md.items():
    cx, cy = counts[x], counts[y]
    s = {f for f in cx if cx[f] == 1 and cy.get(f) == 1}
    if len(s) < 10: continue
    oA = [f for f in ordered[x] if f in s]
    oB = [f for f in ordered[y] if f in s]
    idx = {f: i for i, f in enumerate(oB)}
    perm = [idx[f] for f in oA]
    n = len(perm)
    dist = mtf(perm)
    for lo, hi, lab in BINS:
        if lo <= d < hi:
            r = res[lab]
            r["n"] += 1; r["shared"].append(n)
            r["obs"].append(dist/(n-1)); r["exp"].append(rand_exp(n)/(n-1))
            if dist == 0: r["zero"] += 1
            break

print(f"{'ANI band':<26}{'pairs':>7}{'med core':>10}{'obs':>8}{'random':>8}"
      f"{'ratio':>8}{'identical':>11}")
print("-"*80)
for _, _, lab in BINS:
    r = res[lab]
    if not r["n"]: continue
    o, e = statistics.mean(r["obs"]), statistics.mean(r["exp"])
    print(f"{lab:<26}{r['n']:>7}{statistics.median(r['shared']):>10.0f}"
          f"{o:>8.3f}{e:>8.3f}{o/e:>8.3f}{r['zero']/r['n']*100:>10.1f}%")
print("-"*80)

with open(f"{A}/cayley_vs_ani.csv","w",newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["genome_a","genome_b","species_a","species_b","mash_dist",
                "n_shared","mtf_distance","random_mean"])
    for (x, y), d in md.items():
        cx, cy = counts[x], counts[y]
        s = {f for f in cx if cx[f] == 1 and cy.get(f) == 1}
        if len(s) < 10: continue
        oA = [f for f in ordered[x] if f in s]
        oB = [f for f in ordered[y] if f in s]
        idx = {f: i for i, f in enumerate(oB)}
        perm = [idx[f] for f in oA]
        w.writerow([x, y, sp.get(x,"?"), sp.get(y,"?"), f"{d:.5f}",
                    len(perm), mtf(perm), f"{rand_exp(len(perm)):.2f}"])
print("wrote analysis/cayley_vs_ani.csv")
