#!/usr/bin/env python3
import glob, os
import numpy as np
from collections import defaultdict

PROJ = "/scratch/user/uqcngu19/vibrio-integron-pipeline"
REF  = "GCF_000006745.1"

def if_sites(acc):
    out = []
    for p in glob.glob(f"{PROJ}/results/scale/{acc}/**/*.integrons", recursive=True):
        hdr = None
        for line in open(p):
            if line.startswith("#"): continue
            f = line.rstrip("\n").split("\t")
            if hdr is None: hdr = f; continue
            d = dict(zip(hdr, f))
            if d.get("type_elt","").lower() != "attc": continue
            b, e = int(float(d["pos_beg"])), int(float(d["pos_end"]))
            out.append((d["ID_replicon"], min(b,e), max(b,e)))
    return sorted(set(out))

def hat_hits(acc):
    p = glob.glob(f"{PROJ}/results/hattci/{acc}*.hattci.tsv")
    if not p: return []
    out = []
    for line in open(p[0]):
        q = line.split()
        if not q or q[0] == "hit": continue
        try:
            rid = q[1].split("|")[-1] if "|" in q[1] else q[1]
            s, e, v = int(q[2]), int(q[3]), float(q[4])
        except (ValueError, IndexError): continue
        out.append((rid, min(s,e), max(s,e), v))
    return out

def matched(q, truth, slop=30):
    idx = defaultdict(list)
    for r,a,b,*_ in q: idx[r].append((a,b))
    return sum(1 for r,a,b in truth
               if any(y >= a-slop and x <= b+slop for x,y in idx.get(r,[])))

t, h = if_sites(REF), hat_hits(REF)
vs = np.array([x[3] for x in h])
print(f"N16961: {len(t)} IF sites, {len(h)} HattCI hits")
print(f"Vscore range {vs.min():.2f} to {vs.max():.2f}; "
      f"median {np.median(vs):.2f}; {(vs>0).sum()} hits above 0\n")

print(f"{'thr':>7}{'hits':>7}{'match':>7}{'recall':>9}{'prec':>8}{'F1':>8}")
print("-"*46)
grid = list(np.arange(-3, 3.01, 0.25)) + [4,5,6,8,10]
best = (None, -1)
for thr in grid:
    sub = [x for x in h if x[3] >= thr]
    if not sub: continue
    m = matched(sub, t)
    rec, pre = m/len(t), m/len(sub)
    f1 = 2*rec*pre/(rec+pre) if rec+pre else 0
    print(f"{thr:>7.2f}{len(sub):>7}{m:>7}{rec*100:>8.1f}%{pre*100:>7.1f}%{f1:>8.3f}")
    if f1 > best[1]: best = (thr, f1)

for THR, label in [(0.0, "Vscore > 0 (untuned: model beats null)"),
                   (best[0], f"max-F1 on N16961 ({best[0]:.2f})")]:
    accs = sorted({os.path.basename(p).split(".hattci")[0]
                   for p in glob.glob(f"{PROJ}/results/hattci/*.hattci.tsv")})
    TI = TM = TH = 0
    poor = []
    for acc in accs:
        if acc == REF: continue
        tt = if_sites(acc)
        if len(tt) < 20: continue          # skip trivial arrays
        sub = [x for x in hat_hits(acc) if x[3] >= THR]
        m = matched(sub, tt)
        TI += len(tt); TM += m; TH += len(sub)
        if m/len(tt) < 0.5: poor.append((acc, len(tt), m))
    print(f"\n=== {label} ===")
    print(f"  genomes (>=20 attC)  : {sum(1 for a in accs if a!=REF and len(if_sites(a))>=20)}")
    print(f"  IF sites             : {TI}")
    print(f"  recovered by HattCI  : {TM} ({TM/TI*100:.1f}%)")
    print(f"  HattCI hits total    : {TH}  (excess {TH-TM}, {TH/TI:.1f}x)")
    if poor: print(f"  poor agreement (<50%): {poor}")
