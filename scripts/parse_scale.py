#!/usr/bin/env python3
import glob, os, re, csv, statistics
from collections import defaultdict

PROJ = "/scratch/user/uqcngu19/vibrio-integron-pipeline"
paths = sorted(glob.glob(f"{PROJ}/results/scale/**/*.integrons", recursive=True))
print(f"parsing {len(paths)} .integrons files\n")

rows = []
for p in paths:
    m = re.search(r"results/scale/([^/]+)/", p)
    acc = m.group(1) if m else os.path.basename(p)
    hdr, recs = None, []
    for line in open(p):
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if hdr is None: hdr = f; continue
        recs.append(dict(zip(hdr, f)))
    groups = defaultdict(list)
    for r in recs:
        groups[(r.get("ID_replicon","?"), r.get("ID_integron","?"))].append(r)
    if not groups:
        rows.append(dict(accession=acc, replicon="", integron="", type="none",
                         intI="no", attC=0, cassettes=0, start=0, end=0, span_kb=0))
        continue
    for (rep, iid), mem in groups.items():
        attc = [x for x in mem if x.get("type_elt","").lower()=="attc"]
        inti = [x for x in mem if "inti" in (x.get("annotation","") or "").lower()]
        pos  = [(int(float(x["pos_beg"])), int(float(x["pos_end"])))
                for x in mem if x.get("pos_beg") and x.get("pos_end")]
        lo = min(a for a,_ in pos) if pos else 0
        hi = max(b for _,b in pos) if pos else 0
        rows.append(dict(accession=acc, replicon=rep, integron=iid,
            type=next((x.get("type","") for x in mem if x.get("type")), "?"),
            intI="yes" if inti else "no",
            attC=len(attc), cassettes=max(len(attc)-1,0),
            start=lo, end=hi, span_kb=round((hi-lo)/1000,1)))

out = f"{PROJ}/results/scale_summary.csv"
with open(out,"w",newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)

per = defaultdict(lambda: dict(attC=0, cass=0, integrons=0, complete=0, calin=0, in0=0))
for r in rows:
    d = per[r["accession"]]
    d["attC"] += r["attC"]; d["cass"] += r["cassettes"]
    if r["type"] != "none":
        d["integrons"] += 1
        if r["type"] == "complete": d["complete"] += 1
        elif r["type"] == "CALIN":  d["calin"] += 1
        elif r["type"] == "In0":    d["in0"] += 1

sizes   = [d["attC"] for d in per.values()]
withany = [s for s in sizes if s > 0]
big     = [s for s in sizes if s >= 20]
print("="*58)
print(f"genomes parsed        : {len(per)}")
print(f"with >=1 attC         : {len(withany)} ({len(withany)/len(per)*100:.1f}%)")
print(f"with >=20 attC (SI)   : {len(big)} ({len(big)/len(per)*100:.1f}%)")
print(f"no integron at all    : {len(per)-len(withany)}")
print(f"TOTAL attC            : {sum(sizes)}")
print(f"TOTAL cassettes       : {sum(d['cass'] for d in per.values())}")
print(f"genomes w/ complete   : {sum(1 for d in per.values() if d['complete'])}")
print(f"genomes w/ CALIN only : {sum(1 for d in per.values() if d['calin'] and not d['complete'])}")
if withany:
    print(f"\namong genomes with >=1 attC:")
    print(f"  mean   : {statistics.mean(withany):.1f}")
    print(f"  median : {statistics.median(withany):.0f}")
    print(f"  max    : {max(withany)}")
top = sorted(per.items(), key=lambda kv: -kv[1]["attC"])[:12]
print("\ntop 12 largest arrays:")
for a,d in top:
    print(f"  {a}  attC={d['attC']:4d}  cass={d['cass']:4d}  "
          f"complete={d['complete']} CALIN={d['calin']}")
print("="*58)
print(f"wrote {out}")
