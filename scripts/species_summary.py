#!/usr/bin/env python3
import csv, statistics
from collections import defaultdict

PROJ = "/scratch/user/uqcngu19/vibrio-integron-pipeline"

sp = {}
with open(f"{PROJ}/assembly_metadata.tsv") as fh:
    next(fh)
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2:
            sp[p[0].strip()] = p[1].strip()

per = defaultdict(lambda: dict(attC=0, cass=0, complete=0, calin=0))
with open(f"{PROJ}/results/scale_summary.csv") as fh:
    for r in csv.DictReader(fh):
        d = per[r["accession"]]
        d["attC"] += int(r["attC"]); d["cass"] += int(r["cassettes"])
        if r["type"] == "complete": d["complete"] += 1
        elif r["type"] == "CALIN":  d["calin"] += 1

missing = [a for a in per if a not in sp]
print(f"accessions with no species match: {len(missing)}")
if missing[:3]: print("  e.g.", missing[:3])

byspec = defaultdict(list)
for acc, d in per.items():
    byspec[sp.get(acc, "UNKNOWN")].append(d["attC"])

print("\n" + "="*78)
print(f"{'species':<38}{'n':>4}{'%SI':>7}{'mean':>8}{'median':>8}{'max':>7}")
print("="*78)
rank = sorted(byspec.items(), key=lambda kv: -len(kv[1]))
for name, vals in rank:
    if len(vals) < 3: continue
    si = sum(1 for v in vals if v >= 20)
    print(f"{name[:37]:<38}{len(vals):>4}{si/len(vals)*100:>6.0f}%"
          f"{statistics.mean(vals):>8.1f}{statistics.median(vals):>8.0f}{max(vals):>7}")
print("="*78)

print("\ntop 20 arrays with species:")
top = sorted(per.items(), key=lambda kv: -kv[1]["attC"])[:20]
for acc, d in top:
    print(f"  {acc:<20} {sp.get(acc,'?')[:40]:<42} attC={d['attC']:>4} "
          f"complete={d['complete']} CALIN={d['calin']}")

zero = [(acc, sp.get(acc,'?')) for acc, d in per.items() if d["attC"] == 0]
print(f"\n{len(zero)} genomes with NO attC — species breakdown:")
z = defaultdict(int)
for _, s in zero: z[s] += 1
for s, n in sorted(z.items(), key=lambda kv: -kv[1])[:15]:
    print(f"  {n:>3}  {s}")

with open(f"{PROJ}/results/per_genome_with_species.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["accession","species","attC","cassettes","complete","calin"])
    for acc, d in sorted(per.items(), key=lambda kv: -kv[1]["attC"]):
        w.writerow([acc, sp.get(acc,"?"), d["attC"], d["cass"], d["complete"], d["calin"]])
print(f"\nwrote results/per_genome_with_species.csv")
