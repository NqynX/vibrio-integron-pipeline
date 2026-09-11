import csv, json, os, sys
from collections import defaultdict
INP=sys.argv[1] if len(sys.argv)>1 else "analysis/cayley_vs_ani.csv"
OUT=sys.argv[2] if len(sys.argv)>2 else "network_data.json"
MIN_SHARED=int(os.environ.get("MIN_SHARED","12"))
MAX_NORM=float(os.environ.get("MAX_NORM","0.55"))
TOP_K=int(os.environ.get("TOP_K","6"))
GENUS={"Vibrio":"V.","Aliivibrio":"A.","Photobacterium":"P.","Grimontia":"G.","Salinivibrio":"S."}
def binom(s):
    t=s.split()
    return f"{GENUS.get(t[0],t[0][:1]+'.')} {t[1]}" if len(t)>=2 else s
def label(s):
    for f,a in GENUS.items(): s=s.replace(f+" ",a+" ")
    return s.split(" str. ")[0].split(" subsp. ")[0][:38]
meta,cand,total={},{},0
for r in csv.DictReader(open(INP,newline="")):
    try:
        nsh=int(r["n_shared"]);mtf=int(r["mtf_distance"]);mash=float(r["mash_dist"]);rnd=float(r["random_mean"])
    except (ValueError,KeyError): continue
    total+=1
    if nsh<MIN_SHARED or rnd<=0: continue
    norm=mtf/rnd
    if norm>MAX_NORM: continue
    a,b=r["genome_a"],r["genome_b"]
    meta.setdefault(a,{"s":binom(r["species_a"]),"l":label(r["species_a"])})
    meta.setdefault(b,{"s":binom(r["species_b"]),"l":label(r["species_b"])})
    k=frozenset((a,b))
    if k in cand and cand[k]["shared"]>=nsh: continue
    cand[k]=dict(source=a,target=b,mtf=mtf,shared=nsh,mash=round(mash,5),norm=round(norm,3))
inc=defaultdict(list)
for e in cand.values(): inc[e["source"]].append(e);inc[e["target"]].append(e)
keepK=set()
for node,es in inc.items():
    for e in sorted(es,key=lambda x:x["mtf"])[:TOP_K]: keepK.add(frozenset((e["source"],e["target"])))
edges=[cand[k] for k in keepK]
used=set()
for e in edges: used.add(e["source"]);used.add(e["target"])
nodes=[dict(id=a,name=meta[a]["l"],acc=a,species=meta[a]["s"]) for a in used]
json.dump(dict(metric="mtf_cayley",
  description=("Vibrio superintegrons connected where cassette ORDER is conserved (move-to-front "
    "distance well below the random-permutation baseline; a path on the Cayley graph of the symmetric "
    f"group). Signal core: n_shared>={MIN_SHARED}, mtf/random<={MAX_NORM}, top-{TOP_K} per genome."),
  filter=dict(min_shared=MIN_SHARED,max_norm=MAX_NORM,top_k=TOP_K),
  nodes=nodes,edges=edges),open(OUT,"w"),separators=(",",":"))
sp=defaultdict(int)
for n in nodes: sp[n["species"]]+=1
print(f"scanned {total} pairs -> kept {len(nodes)} nodes, {len(edges)} edges -> {OUT} ({os.path.getsize(OUT)} bytes)")
print("species:",dict(sorted(sp.items(),key=lambda x:-x[1])[:10]))
