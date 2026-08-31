#!/usr/bin/env python3
import glob, csv, os
from collections import Counter, defaultdict

PROJ = "/scratch/user/uqcngu19/vibrio-integron-pipeline"
sp = {r["accession"]: r["species"] for r in
      csv.DictReader(open(f"{PROJ}/results/per_genome_with_species.csv"))}
cass = {r["accession"]: int(r["cassettes"]) for r in
        csv.DictReader(open(f"{PROJ}/results/per_genome_with_species.csv"))}

types, per_gen, rows = Counter(), {}, []
files = glob.glob(f"{PROJ}/analysis/defense_per_genome/*/*_systems.tsv")
for f in files:
    acc = os.path.basename(os.path.dirname(f))
    n = 0
    with open(f) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            t = r.get("subtype") or r.get("type")
            if t: types[t] += 1; n += 1
    per_gen[acc] = n
    rows.append((acc, sp.get(acc,"?"), cass.get(acc,0), n))

tot = sum(per_gen.values())
withdef = sum(1 for v in per_gen.values() if v > 0)
print("="*60)
print(f"genomes scanned            : {len(per_gen)}")
print(f"genomes with >=1 system    : {withdef} ({withdef/len(per_gen)*100:.1f}%)")
print(f"total defence systems      : {tot}")
print(f"mean systems per genome    : {tot/len(per_gen):.1f}")
print(f"max in one genome          : {max(per_gen.values())}")
print("="*60)
print("\nTOP SYSTEM TYPES")
for t,c in types.most_common(30): print(f"  {c:>5}  {t}")

print("\nBY SPECIES (n>=5 genomes)")
bys = defaultdict(list)
for acc,s,c,n in rows: bys[s].append(n)
for s,v in sorted(bys.items(), key=lambda kv: -len(kv[1])):
    if len(v) < 5: continue
    print(f"  {s[:40]:<42} n={len(v):>4} mean={sum(v)/len(v):>6.1f} max={max(v):>4}")

with open(f"{PROJ}/analysis/defence_per_genome.csv","w",newline="") as fh:
    w = csv.writer(fh); w.writerow(["accession","species","cassettes","defence_systems"])
    w.writerows(sorted(rows, key=lambda r: -r[3]))
print("\nwrote analysis/defence_per_genome.csv")
