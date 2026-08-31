#!/usr/bin/env python3
import csv, glob, os
from collections import defaultdict

PROJ = "/scratch/user/uqcngu19/vibrio-integron-pipeline"
A = f"{PROJ}/analysis"

# protein -> family
fam = {}
for line in open(f"{A}/clust_cluster.tsv"):
    rep, mem = line.rstrip("\n").split("\t")[:2]
    fam[mem] = rep

integrase = {r["protein_id"] for r in csv.DictReader(open(f"{A}/protein_meta.csv"))
             if r["is_integrase"].strip().lower() == "true"}

# defence-positive PROTEINS from the per-genome run
dprot = set()
for f in glob.glob(f"{A}/defense_per_genome/*/*_systems.tsv"):
    for r in csv.DictReader(open(f), delimiter="\t"):
        for p in str(r.get("protein_in_syst","")).split(","):
            p = p.strip()
            if p: dprot.add(p)
print(f"defence-positive proteins (per-genome run): {len(dprot)}")

# Pfam-positive families
pfam_fam = set()
for line in open(f"{A}/pfam_reps.tbl"):
    if line.startswith("#"): continue
    p = line.split(None, 4)
    if len(p) > 3: pfam_fam.add(p[3])

prots = [p for p in fam if p not in integrase]
dfam = {fam[p] for p in dprot if p in fam}
fams = {fam[p] for p in prots}

ann_f = (pfam_fam | dfam) & fams
ann_p = [p for p in prots if fam[p] in pfam_fam or p in dprot or fam[p] in dfam]

print("="*58)
print(f"cassette proteins (non-integrase) : {len(prots)}")
print(f"families                          : {len(fams)}")
print(f"  Pfam-positive families          : {len(pfam_fam & fams)} ({len(pfam_fam & fams)/len(fams)*100:.1f}%)")
print(f"  defence-positive families       : {len(dfam)} ({len(dfam)/len(fams)*100:.1f}%)")
print(f"  annotated (either)              : {len(ann_f)} ({len(ann_f)/len(fams)*100:.1f}%)")
print(f"  DARK families                   : {len(fams)-len(ann_f)} ({(len(fams)-len(ann_f))/len(fams)*100:.1f}%)")
print(f"\n  annotated proteins              : {len(ann_p)} ({len(ann_p)/len(prots)*100:.1f}%)")
print(f"  DARK proteins                   : {len(prots)-len(ann_p)} ({(len(prots)-len(ann_p))/len(prots)*100:.1f}%)")
print("="*58)
