#!/usr/bin/env python3
import glob, re, statistics
from collections import defaultdict

PROJ = "/scratch/user/uqcngu19/vibrio-integron-pipeline"
THRESH = [2000, 4000, 6000, 8000, 10000]

allgaps, per_int = [], []
for p in glob.glob(f"{PROJ}/results/scale/**/*.integrons", recursive=True):
    acc = re.search(r"results/scale/([^/]+)/", p).group(1)
    hdr, recs = None, []
    for line in open(p):
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if hdr is None: hdr = f; continue
        recs.append(dict(zip(hdr, f)))
    g = defaultdict(list)
    for r in recs:
        if r.get("type_elt","").lower() == "attc":
            g[(r["ID_replicon"], r["ID_integron"])].append(
                (int(float(r["pos_beg"])), int(float(r["pos_end"]))))
    for k, v in g.items():
        if len(v) < 2: continue
        v.sort()
        gaps = [v[i+1][0] - v[i][1] for i in range(len(v)-1)]
        gaps = [x for x in gaps if x >= 0]
        if not gaps: continue
        allgaps += gaps
        per_int.append((acc, k[1], len(v), max(gaps), statistics.median(gaps)))

allgaps.sort()
n = len(allgaps)
print(f"{n} inter-attC gaps across {len(per_int)} integrons\n")
print("gap percentiles (bp):")
for q in [50, 75, 90, 95, 99, 99.9]:
    print(f"  p{q:<5} {allgaps[int(n*q/100)-1]:>8}")
print(f"  max    {allgaps[-1]:>8}")

print("\nhow many integrons contain a gap larger than T")
print("(these are the ones -dt 10000 merged that a smaller T would split):")
for t in THRESH:
    k = sum(1 for r in per_int if r[3] > t)
    print(f"  gap > {t:>6} : {k:>5} integrons ({k/len(per_int)*100:5.1f}%)")

print("\nintegrons with the largest internal gaps (suspect merges):")
for r in sorted(per_int, key=lambda x: -x[3])[:12]:
    print(f"  {r[0]:<20} {r[1]:<14} attC={r[2]:>4} max_gap={r[3]:>7} median_gap={r[4]:>6.0f}")
