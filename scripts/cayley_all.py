#!/usr/bin/env python3
"""All-vs-all: how comparable are arrays within vs between species?"""
import csv, random, statistics
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

ordered = {}
for acc, it in arrays.items():
    it.sort()
    if acc in ipos and it and ipos[acc] > it[-1][0]: it = it[::-1]
    ordered[acc] = [f for _, f in it]

accs = [a for a, v in ordered.items() if len(v) >= 20]
print(f"{len(accs)} genomes with >=20 ordered cassettes\n")

def mtf(perm):
    pos = {v: i for i, v in enumerate(perm)}
    m = len(perm)
    while m > 1 and pos[m-2] < pos[m-1]: m -= 1
    return m - 1

_rc = {}
def rand_exp(n):
    if n not in _rc:
        _rc[n] = statistics.mean(
            mtf(random.sample(range(n), n)) for _ in range(100))
    return _rc[n]

sets = {a: defaultdict(int) for a in accs}
for a in accs:
    for f in ordered[a]: sets[a][f] += 1

pairs = [(accs[i], accs[j]) for i in range(len(accs)) for j in range(i+1, len(accs))]
print(f"{len(pairs)} total pairs; sampling 30000")
if len(pairs) > 30000: pairs = random.sample(pairs, 30000)

stats = defaultdict(lambda: dict(n=0, comparable=0, shared=[], obs=[], exp=[]))
for x, y in pairs:
    same = "same species" if sp.get(x) == sp.get(y) else "different species"
    d = stats[same]; d["n"] += 1
    cx, cy = sets[x], sets[y]
    s = {f for f in cx if cx[f] == 1 and cy.get(f) == 1}
    d["shared"].append(len(s))
    if len(s) < 10: continue
    oA = [f for f in ordered[x] if f in s]
    oB = [f for f in ordered[y] if f in s]
    idx = {f: i for i, f in enumerate(oB)}
    perm = [idx[f] for f in oA]
    n = len(perm)
    d["comparable"] += 1
    d["obs"].append(mtf(perm)/(n-1))
    d["exp"].append(rand_exp(n)/(n-1))

print(f"\n{'relationship':<20}{'pairs':>8}{'comparable':>12}{'%':>7}"
      f"{'med shared':>12}{'obs':>8}{'random':>8}{'ratio':>8}")
print("-"*85)
for k, d in stats.items():
    if not d["n"]: continue
    med = statistics.median(d["shared"])
    if d["obs"]:
        o, e = statistics.mean(d["obs"]), statistics.mean(d["exp"])
        print(f"{k:<20}{d['n']:>8}{d['comparable']:>12}"
              f"{d['comparable']/d['n']*100:>6.1f}%{med:>12.0f}{o:>8.3f}{e:>8.3f}{o/e:>8.3f}")
    else:
        print(f"{k:<20}{d['n']:>8}{d['comparable']:>12}"
              f"{d['comparable']/d['n']*100:>6.1f}%{med:>12.0f}       -       -       -")
