#!/usr/bin/env python3
"""Move-to-front permutation distances between Vibrio cassette arrays."""
import csv, random, statistics, sys
from collections import defaultdict

PROJ = "/scratch/user/uqcngu19/vibrio-integron-pipeline"
A = f"{PROJ}/analysis"
SPECIES = sys.argv[1] if len(sys.argv) > 1 else "Vibrio cholerae"
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
    acc = r["accession"]
    if sp.get(acc) != SPECIES: continue
    g = gene_index(r["protein_id"])
    if g < 0: continue
    if r["is_integrase"].strip().lower() == "true":
        ipos[acc] = g; continue
    if r["protein_id"] in fam:
        arrays[acc].append((g, fam[r["protein_id"]]))

ordered = {}
for acc, it in arrays.items():
    it.sort()
    if acc in ipos and it and ipos[acc] > it[-1][0]: it = it[::-1]
    ordered[acc] = [f for _, f in it]

accs = [a for a, v in ordered.items() if len(v) >= 20]
print(f"{SPECIES}: {len(accs)} genomes with >=20 ordered cassettes")
print(f"  oriented by integrase: {sum(1 for a in accs if a in ipos)}\n")
if len(accs) < 2: sys.exit("not enough genomes")

def mtf(perm):
    pos = {v: i for i, v in enumerate(perm)}
    m = len(perm)
    while m > 1 and pos[m-2] < pos[m-1]: m -= 1
    return m - 1

_rc = {}
def rand_exp(n):
    if n not in _rc:
        _rc[n] = statistics.mean(mtf(random.sample(range(n), n)) for _ in range(200))
    return _rc[n]

def compare(a, b):
    ca, cb = defaultdict(int), defaultdict(int)
    for f in a: ca[f] += 1
    for f in b: cb[f] += 1
    s = {f for f in ca if ca[f] == 1 and cb.get(f) == 1}
    if len(s) < 10: return None
    oA = [f for f in a if f in s]
    oB = [f for f in b if f in s]
    idx = {f: i for i, f in enumerate(oB)}
    perm = [idx[f] for f in oA]
    return len(perm), mtf(perm), rand_exp(len(perm))

pairs = [(accs[i], accs[j]) for i in range(len(accs)) for j in range(i+1, len(accs))]
if len(pairs) > 5000: pairs = random.sample(pairs, 5000)
rows = [r for r in (compare(ordered[x], ordered[y]) for x, y in pairs) if r]
if not rows: sys.exit("no comparable pairs")

ns  = [r[0] for r in rows]
obs = [r[1]/(r[0]-1) for r in rows]
exp = [r[2]/(r[0]-1) for r in rows]
print(f"comparable pairs        : {len(rows)} of {len(pairs)}")
print(f"shared single-copy core : median {statistics.median(ns):.0f} "
      f"(range {min(ns)}-{max(ns)})\n")
print(f"OBSERVED normalised distance : mean {statistics.mean(obs):.3f}  "
      f"median {statistics.median(obs):.3f}")
print(f"RANDOM expectation           : mean {statistics.mean(exp):.3f}")
print(f"ratio observed/random        : {statistics.mean(obs)/statistics.mean(exp):.3f}")
print(f"identical order (distance 0) : {sum(1 for r in rows if r[1]==0)}/{len(rows)}")

with open(f"{A}/cayley_{SPECIES.replace(' ','_')}.csv","w",newline="") as fh:
    w = csv.writer(fh); w.writerow(["n_shared","mtf_distance","random_mean"])
    w.writerows(rows)
print(f"\nwrote analysis/cayley_{SPECIES.replace(' ','_')}.csv")
