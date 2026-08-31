#!/usr/bin/env python3
import csv, statistics, sys
from collections import defaultdict

PROJ = "/scratch/user/uqcngu19/vibrio-integron-pipeline"
rows = []
for sp in ["Vibrio_cholerae", "Vibrio_parahaemolyticus", "Vibrio_vulnificus"]:
    try:
        for r in csv.DictReader(open(f"{PROJ}/analysis/cayley_{sp}.csv")):
            rows.append((sp, int(r["n_shared"]), int(r["mtf_distance"]),
                         float(r["random_mean"])))
    except FileNotFoundError:
        pass
if not rows:
    sys.exit("run scripts/cayley.py per species first")

print(f"{'shared core':<14}{'pairs':>7}{'obs':>9}{'random':>9}{'ratio':>8}{'dist=0':>9}")
print("-"*58)
for lo, hi, lab in [(10,19,"10-19"),(20,49,"20-49"),(50,99,"50-99"),(100,10**9,"100+")]:
    sub = [r for r in rows if lo <= r[1] <= hi]
    if not sub: continue
    o = statistics.mean(r[2]/(r[1]-1) for r in sub)
    e = statistics.mean(r[3]/(r[1]-1) for r in sub)
    z = sum(1 for r in sub if r[2] == 0)
    print(f"{lab:<14}{len(sub):>7}{o:>9.3f}{e:>9.3f}{o/e:>8.3f}{z:>9}")
print("-"*58)
print(f"{'ALL':<14}{len(rows):>7}"
      f"{statistics.mean(r[2]/(r[1]-1) for r in rows):>9.3f}"
      f"{statistics.mean(r[3]/(r[1]-1) for r in rows):>9.3f}")
