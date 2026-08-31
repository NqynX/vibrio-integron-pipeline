#!/usr/bin/env python3
import glob, os, re
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

def recovered(q, truth, slop=30):
    idx = defaultdict(list)
    for r, a, b, *_ in q: idx[r].append((a,b))
    n = 0
    for r, a, b in truth:
        if any(y >= a-slop and x <= b+slop for x,y in idx.get(r, [])): n += 1
    return n

accs = sorted({os.path.basename(p).split(".hattci")[0]
               for p in glob.glob(f"{PROJ}/results/hattci/*.hattci.tsv")})
print(f"{len(accs)} genomes with HattCI output\n")

# sanity: do HattCI hits land where IntegronFinder's array is?
t = if_sites(REF); h = hat_hits(REF)
if t and h:
    rep = t[0][0]
    hs = [x for x in h if x[0] == rep]
    print(f"COORDINATE SANITY CHECK ({REF}, {rep})")
    print(f"  IF array span     : {min(a for _,a,_ in t)} - {max(b for _,_,b in t)}  (n={len(t)})")
    if hs:
        print(f"  HattCI hit span   : {min(a for _,a,_,_ in hs)} - {max(b for _,_,b,_ in hs)}  (n={len(hs)})")
    print()

# calibrate on N16961 only
print("CALIBRATION on N16961 only")
best = 0.0
for thr in np.arange(0, 30, 1.0):
    sub = [x for x in h if x[3] >= thr]
    r = recovered(sub, t)
    rec = r/len(t)*100
    print(f"  Vscore>={thr:>4.0f}  hits={len(sub):>5}  recovered={r:>4}/{len(t)}  ({rec:5.1f}%)")
    if rec >= 90: best = thr
THR = best
print(f"\n>>> frozen threshold: Vscore >= {THR}\n")

print(f"{'genome':<20}{'IF':>6}{'HattCI':>8}{'match':>7}{'%':>7}{'extra':>7}")
print("-"*57)
tot_if = tot_ok = tot_hat = 0
for acc in accs:
    t = if_sites(acc)
    if not t: continue
    sub = [x for x in hat_hits(acc) if x[3] >= THR]
    ok = recovered(sub, t)
    tag = "  <- calibration" if acc == REF else ""
    print(f"{acc:<20}{len(t):>6}{len(sub):>8}{ok:>7}{ok/len(t)*100:>6.1f}%{len(sub)-ok:>7}{tag}")
    if acc != REF:
        tot_if += len(t); tot_ok += ok; tot_hat += len(sub)
print("-"*57)
if tot_if:
    print(f"\nINDEPENDENT SET ({len(accs)-1} genomes, threshold not tuned on them)")
    print(f"  IntegronFinder attC sites : {tot_if}")
    print(f"  recovered by HattCI       : {tot_ok} ({tot_ok/tot_if*100:.1f}%)")
    print(f"  HattCI-only hits          : {tot_hat-tot_ok}")
    print(f"\n  >>> ABSTRACT: HattCI independently recovered {tot_ok/tot_if*100:.0f}% of sites")
