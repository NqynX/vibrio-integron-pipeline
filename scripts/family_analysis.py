#!/usr/bin/env python3
import csv, statistics
from collections import defaultdict, Counter

PROJ = "/scratch/user/uqcngu19/vibrio-integron-pipeline"

# species per accession
sp = {}
with open(f"{PROJ}/results/per_genome_with_species.csv") as fh:
    for r in csv.DictReader(fh):
        sp[r["accession"]] = r["species"]

# integrases to exclude
integrase = set()
with open(f"{PROJ}/analysis/protein_meta.csv") as fh:
    for r in csv.DictReader(fh):
        if r["is_integrase"].strip().lower() == "true":
            integrase.add(r["protein_id"])
print(f"excluding {len(integrase)} integrase proteins\n")

fam_gen = defaultdict(set)      # family -> genomes
fam_spec = defaultdict(set)     # family -> species
fam_count = Counter()           # family -> total proteins
dup = defaultdict(Counter)      # genome -> family -> copies
n = 0
with open(f"{PROJ}/analysis/clust_cluster.tsv") as fh:
    for line in fh:
        rep, mem = line.rstrip("\n").split("\t")[:2]
        if mem in integrase: continue
        acc = mem.split("|")[0]
        fam_gen[rep].add(acc)
        fam_spec[rep].add(sp.get(acc, "?"))
        fam_count[rep] += 1
        dup[acc][rep] += 1
        n += 1

F = len(fam_gen)
print("="*60)
print(f"cassette proteins (non-integrase) : {n}")
print(f"families                          : {F}")
print("="*60)

ng = Counter(len(v) for v in fam_gen.values())
ns = Counter(len(v) for v in fam_spec.values())
print("\nFAMILIES BY NUMBER OF GENOMES")
priv_g = ng[1]
print(f"  in exactly 1 genome : {priv_g:>6} ({priv_g/F*100:5.1f}%)")
for lo, hi, lab in [(2,2,"2"),(3,5,"3-5"),(6,10,"6-10"),(11,50,"11-50"),(51,10**9,">50")]:
    k = sum(v for kk,v in ng.items() if lo <= kk <= hi)
    print(f"  in {lab:<18}: {k:>6} ({k/F*100:5.1f}%)")

print("\nFAMILIES BY NUMBER OF SPECIES")
priv_s = ns[1]
print(f"  in exactly 1 species: {priv_s:>6} ({priv_s/F*100:5.1f}%)")
for lo, hi, lab in [(2,2,"2"),(3,5,"3-5"),(6,10,"6-10"),(11,10**9,">10")]:
    k = sum(v for kk,v in ns.items() if lo <= kk <= hi)
    print(f"  in {lab:<18}: {k:>6} ({k/F*100:5.1f}%)")

print("\nWITHIN-GENOME DUPLICATION")
dupfam = [(g,f,c) for g,d in dup.items() for f,c in d.items() if c > 1]
gwith = len({g for g,_,_ in dupfam})
print(f"  genomes with >=1 duplicated family : {gwith} / {len(dup)} ({gwith/len(dup)*100:.1f}%)")
print(f"  duplicated family instances        : {len(dupfam)}")
if dupfam:
    print(f"  max copies of one family in a genome: {max(c for _,_,c in dupfam)}")
    print(f"  median copies (where >1)           : {statistics.median([c for _,_,c in dupfam]):.0f}")

print("\nMOST WIDESPREAD FAMILIES (genomes / species / total copies)")
for rep in sorted(fam_gen, key=lambda r: -len(fam_gen[r]))[:15]:
    print(f"  {rep[:52]:<54} g={len(fam_gen[rep]):>4} s={len(fam_spec[rep]):>3} n={fam_count[rep]:>5}")

with open(f"{PROJ}/analysis/family_spread.csv","w",newline="") as fh:
    w = csv.writer(fh); w.writerow(["family","n_genomes","n_species","n_proteins"])
    for rep in fam_gen:
        w.writerow([rep, len(fam_gen[rep]), len(fam_spec[rep]), fam_count[rep]])
print(f"\nwrote analysis/family_spread.csv")
