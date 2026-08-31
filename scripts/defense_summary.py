#!/usr/bin/env python3
import csv
from collections import Counter, defaultdict

PROJ = "/scratch/user/uqcngu19/vibrio-integron-pipeline"
A = f"{PROJ}/analysis"

# family -> spread
spread = {}
with open(f"{A}/family_spread.csv") as fh:
    for r in csv.DictReader(fh):
        spread[r["family"]] = (int(r["n_genomes"]), int(r["n_species"]), int(r["n_proteins"]))

# representative -> defence type
dtype = {}
with open(f"{A}/defense/clust_rep_seq_defense_finder_systems.tsv") as fh:
    rd = csv.DictReader(fh, delimiter="\t")
    for r in rd:
        for p in str(r.get("protein_in_syst","")).split(","):
            p = p.strip()
            if p: dtype[p] = r.get("subtype") or r.get("type")

F = len(spread)
tot_prot = sum(v[2] for v in spread.values())
dfam = {f for f in dtype if f in spread}
dprot = sum(spread[f][2] for f in dfam)

print("="*62)
print(f"families                    : {F}")
print(f"defence-system families     : {len(dfam)} ({len(dfam)/F*100:.1f}%)")
print(f"cassette proteins           : {tot_prot}")
print(f"in defence families         : {dprot} ({dprot/tot_prot*100:.1f}%)")
print("="*62)

print("\nDEFENCE SYSTEM TYPES (by number of families)")
for t, c in Counter(dtype[f] for f in dfam).most_common(25):
    print(f"  {c:>4}  {t}")

print("\nIS DEFENCE ENRICHED AMONG WIDESPREAD FAMILIES?")
bins = [(1,1,"1 species (private)"),(2,2,"2"),(3,5,"3-5"),(6,10,"6-10"),(11,10**9,">10 species")]
for lo, hi, lab in bins:
    fam = [f for f,v in spread.items() if lo <= v[1] <= hi]
    d = sum(1 for f in fam if f in dfam)
    print(f"  {lab:<22} n={len(fam):>6}  defence={d:>5} ({d/len(fam)*100 if fam else 0:5.1f}%)")

print("\nTOP 20 MOST WIDESPREAD FAMILIES — defence status")
top = sorted(spread.items(), key=lambda kv: -kv[1][1])[:20]
for f, (g, s, n) in top:
    print(f"  s={s:>3} g={g:>4} n={n:>5}  {dtype.get(f,'-'):<28} {f[:40]}")
