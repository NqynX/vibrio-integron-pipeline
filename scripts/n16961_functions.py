import csv, glob
from collections import defaultdict
PROJ="/scratch/user/uqcngu19/vibrio-integron-pipeline"
ACC="GCF_000006745.1"; RID="NC_002506.1"
pid=lambda el: ACC+"|"+RID+"|"+el
elements=[]; hdr=None
for line in open(PROJ+"/data/example/N16961.integrons"):
    if line.startswith("#"): continue
    f=line.rstrip("\n").split("\t")
    if hdr is None: hdr=f; continue
    d=dict(zip(hdr,f))
    if d["type_elt"].lower()=="protein" and d["annotation"].lower()!="inti":
        elements.append(d["element"])
defence={}
gf=glob.glob(PROJ+"/analysis/defense_per_genome/"+ACC+"/*genes.tsv")
if gf:
    for r in csv.DictReader(open(gf[0]), delimiter="\t"):
        defence[r["hit_id"]]=r.get("subtype") or r.get("type") or "defence"
m2rep={}
for line in open(PROJ+"/analysis/clust_cluster.tsv"):
    rep,mem=line.rstrip("\n").split("\t")[:2]; m2rep[mem]=rep
pfam=defaultdict(list)
for line in open(PROJ+"/analysis/pfam_reps.tbl"):
    if line.startswith("#") or not line.strip(): continue
    p=line.split()
    if len(p)<23: continue
    pfam[p[3]].append((p[0]," ".join(p[22:])))
TA=["pare","pard","phdyefm","phd","yefm","taca","higa","higb","rele","relb",
"vapc","vapb","maze","mazf","pemk","hica","hicb","ccda","ccdb","brnt","brna",
"hipa","hipb","zeta","antitoxin","toxin-antitoxin","abieii"]
AMR=["lactamas","aminoglycoside","chloramphenicol","dihydrofolate","dihydropteroate",
"quinolone","qnr","tetracycline","sulfona","rifampin","macrolide","dfra","qace"]
def cls(dom,desc):
    s=(dom+" "+desc).lower()
    if any(k in s for k in TA): return "ta"
    if any(k in s for k in AMR): return "amr"
    return "other"
rows=[]; counts=defaultdict(int)
for el in elements:
    P=pid(el)
    if P in defence:
        cat,label="defence",defence[P]
    else:
        doms=pfam.get(m2rep.get(P),[])
        if doms:
            cat,label="other",doms[0][0]
            for dom,desc in doms:
                c=cls(dom,desc)
                if c in ("ta","amr"): cat,label=c,dom; break
        else:
            cat,label="hypothetical","no annotation"
    rows.append((el,cat,label)); counts[cat]+=1
out=PROJ+"/data/example/N16961.functions.tsv"
open(out,"w").write("element\tcategory\tlabel\n"+"\n".join("\t".join(r) for r in rows)+"\n")
print("cassette CDS:",len(rows))
for c in ["defence","ta","amr","other","hypothetical"]: print("  %-13s%d"%(c,counts[c]))
print("wrote",out)
