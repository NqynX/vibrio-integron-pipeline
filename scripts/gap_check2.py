#!/usr/bin/env python3
import glob, re, statistics
from collections import defaultdict

PROJ = "/scratch/user/uqcngu19/vibrio-integron-pipeline"
GD   = f"{PROJ}/data/complete_vibrionaceae/genomes"

# replicon lengths, for circular-aware gaps
L = {}
for f in glob.glob(f"{GD}/*.fna"):
    name, n = None, 0
    for line in open(f):
        if line.startswith(">"):
            if name: L[name] = n
            name = line[1:].split()[0]; n = 0
        else:
            n += len(line.strip())
    if name: L[name] = n
print(f"{len(L)} replicon lengths loaded\n")

recs = []
for p in glob.glob(f"{PROJ}/results/scale/**/*.integrons", recursive=True):
    acc = re.search(r"results/scale/([^/]+)/", p).group(1)
    hdr, rows = None, []
    for line in open(p):
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if hdr is None: hdr = f; continue
        rows.append(dict(zip(hdr, f)))
    g = defaultdict(list)
    for r in rows:
        if r.get("type_elt","").lower() == "attc":
            g[(r["ID_replicon"], r["ID_integron"])].append(
                (int(float(r["pos_beg"])), int(float(r["pos_end"]))))
    for (rep, iid), v in g.items():
        if len(v) < 3: continue
        v.sort()
        gaps = [v[i+1][0] - v[i][1] for i in range(len(v)-1)]
        gaps = sorted(x for x in gaps if x >= 0)
        if len(gaps) < 2: continue
        wrap = L.get(rep, 0) - v[-1][1] + v[0][0] if rep in L else None
        recs.append(dict(acc=acc, rep=rep, iid=iid, n=len(v),
                         largest=gaps[-1], second=gaps[-2],
                         median=statistics.median(gaps), wrap=wrap))

print("Is the largest gap the origin wrap? (largest ~ replicon_len - array_span)")
wrapped = [r for r in recs if r["wrap"] is not None and r["largest"] > 50000]
print(f"  integrons with a gap > 50 kb: {len(wrapped)}")
for r in sorted(wrapped, key=lambda x: -x["largest"])[:6]:
    print(f"    {r['acc']:<18} n={r['n']:>4} largest={r['largest']:>9} "
          f"wrap_gap={r['wrap']:>9} median={r['median']:>6.0f}")

print("\nSECOND-largest gap = the real internal-gap question")
sec = sorted(r["second"] for r in recs)
n = len(sec)
for q in [50, 75, 90, 95, 99]:
    print(f"  p{q:<4} {sec[int(n*q/100)-1]:>8}")
print(f"  max   {sec[-1]:>8}")
print(f"\nintegrons whose SECOND gap exceeds T (true over-merge exposure):")
for t in [4000, 6000, 8000, 10000]:
    k = sum(1 for r in recs if r["second"] > t)
    print(f"  > {t:>6} : {k:>5} / {len(recs)} ({k/len(recs)*100:5.1f}%)")

print("\nN16961 (literature ground truth: ONE 126 kb array) for comparison:")
for r in recs:
    if r["acc"].startswith("GCF_000006745"):
        print(f"  n={r['n']} largest={r['largest']} second={r['second']} median={r['median']:.0f}")
