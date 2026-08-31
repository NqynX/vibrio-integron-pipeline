#!/usr/bin/env python3
import csv
from collections import Counter, defaultdict

PROJ = "/scratch/user/uqcngu19/vibrio-integron-pipeline"
A = f"{PROJ}/analysis"

spread = {}
for r in csv.DictReader(open(f"{A}/family_spread.csv")):
    spread[r["family"]] = (int(r["n_genomes"]), int(r["n_species"]), int(r["n_proteins"]))

# defence
dtype = {}
with open(f"{A}/defense/clust_rep_seq_defense_finder_systems.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        for p in str(r.get("protein_in_syst","")).split(","):
            p = p.strip()
            if p: dtype[p] = r.get("subtype") or r.get("type")

# pfam: best hit per rep
pf, pfdesc = {}, {}
for line in open(f"{A}/pfam_reps.tbl"):
    if line.startswith("#"): continue
    p = line.split(None, 22)
    if len(p) < 14: continue
    rep, dom, ev = p[3], p[0], float(p[12])
    if rep not in pf or ev < pf[rep][1]:
        pf[rep] = (dom, ev)
        pfdesc[rep] = p[22].strip() if len(p) > 22 else ""

F = len(spread)
tot = sum(v[2] for v in spread.values())
has_pf  = {f for f in pf if f in spread}
has_def = {f for f in dtype if f in spread}
annot   = has_pf | has_def
dark    = set(spread) - annot
dark_p  = sum(spread[f][2] for f in dark)

print("="*64)
print(f"families                       : {F}")
print(f"  with a Pfam domain           : {len(has_pf):>6} ({len(has_pf)/F*100:5.1f}%)")
print(f"  with a defence system        : {len(has_def):>6} ({len(has_def)/F*100:5.1f}%)")
print(f"  with either                  : {len(annot):>6} ({len(annot)/F*100:5.1f}%)")
print(f"  DARK (neither)               : {len(dark):>6} ({len(dark)/F*100:5.1f}%)")
print(f"\ncassette proteins              : {tot}")
print(f"  in dark families             : {dark_p:>6} ({dark_p/tot*100:5.1f}%)")
print("="*64)

print("\nCOMMONEST Pfam DOMAINS (by families)")
for d, c in Counter(pf[f][0] for f in has_pf).most_common(25):
    print(f"  {c:>4}  {d}")

print("\nANNOTATION RATE vs FAMILY SPREAD")
for lo, hi, lab in [(1,1,"1 species"),(2,2,"2"),(3,5,"3-5"),(6,10,"6-10"),(11,10**9,">10 species")]:
    fam = [f for f,v in spread.items() if lo <= v[1] <= hi]
    a = sum(1 for f in fam if f in annot)
    print(f"  {lab:<14} n={len(fam):>6}  annotated={a:>5} ({a/len(fam)*100 if fam else 0:5.1f}%)")

print("\nTOP 20 MOST WIDESPREAD FAMILIES — what are they?")
for f,(g,s,n) in sorted(spread.items(), key=lambda kv: -kv[1][1])[:20]:
    lab = pf[f][0] if f in pf else (dtype.get(f) or "UNANNOTATED")
    print(f"  s={s:>3} g={g:>4} n={n:>5}  {lab[:30]:<32} {pfdesc.get(f,'')[:44]}")
