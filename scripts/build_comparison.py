import csv, glob, json, os
from collections import defaultdict, Counter
PROJ="/scratch/user/uqcngu19/vibrio-integron-pipeline"
GENOMES=[("V. cholerae N16961","GCF_000006745.1"),
         ("V. cholerae 2740-80","GCF_045689275.1"),
         ("V. parahaemolyticus","GCF_040254425.1")]
m2rep={}
for line in open(PROJ+"/analysis/clust_cluster.tsv"):
    rep,mem=line.rstrip("\n").split("\t")[:2]; m2rep[mem]=rep
pfam=defaultdict(list)
for line in open(PROJ+"/analysis/pfam_reps.tbl"):
    if line.startswith("#") or not line.strip(): continue
    p=line.split()
    if len(p)>=23: pfam[p[3]].append((p[0]," ".join(p[22:])))
TA=["pare","pard","phdyefm","phd","yefm","taca","higa","higb","rele","relb","vapc",
"vapb","maze","mazf","pemk","hica","hicb","ccda","ccdb","brnt","brna","hipa","hipb",
"zeta","antitoxin","toxin-antitoxin","abieii","mqsa","yafq"]
AMR=["lactamas","aminoglycoside","chloramphenicol","dihydrofolate","dihydropteroate",
"quinolone","qnr","tetracycline","sulfona","rifampin","macrolide","dfra","qace"]
def cls(dom,desc):
    s=(dom+" "+desc).lower()
    if any(k in s for k in TA): return 2
    if any(k in s for k in AMR): return 4
    return 1
famid={}
def fid(rep):
    if not rep: return -1
    if rep not in famid: famid[rep]=len(famid)
    return famid[rep]
def defence_set(acc):
    d={}
    for gf in glob.glob(PROJ+"/analysis/defense_per_genome/"+acc+"/*genes.tsv"):
        for r in csv.DictReader(open(gf),delimiter="\t"):
            d[r["hit_id"]]=r.get("subtype") or r.get("type") or "defence"
    return d
genomes=[]; famfn={}
for name,acc in GENOMES:
    ig=glob.glob(PROJ+"/results/scale/"+acc+"/**/*.integrons",recursive=True)
    if not ig: print("MISSING",acc); continue
    dset=defence_set(acc); hdr=None; rows=[]
    for line in open(ig[0]):
        if line.startswith("#"): continue
        f=line.rstrip("\n").split("\t")
        if hdr is None: hdr=f; continue
        rows.append(dict(zip(hdr,f)))
    cnt=Counter((r["ID_replicon"],r["ID_integron"]) for r in rows if r["type_elt"].lower()=="attc")
    if not cnt: print("no attc",acc); continue
    (repl,intid),_=cnt.most_common(1)[0]
    sel=[r for r in rows if r["ID_replicon"]==repl and r["ID_integron"]==intid]
    pb=lambda r:int(float(r["pos_beg"])); pe=lambda r:int(float(r["pos_end"]))
    a0=min(pb(r) for r in sel); a1=max(pe(r) for r in sel)
    isinti=lambda r:"inti" in (r["annotation"]+r["model"]).lower()
    inti=[r for r in sel if isinti(r)]
    flip=bool(inti) and (sum(pb(r) for r in inti)/len(inti)>(a0+a1)/2)
    out=[]
    for r in sel:
        if r["type_elt"].lower()=="attc": continue
        s,e=pb(r),pe(r); d=-1 if r["strand"]=="-1" else 1
        if flip: s,e=a1-e+a0,a1-s+a0; d=-d
        if isinti(r): out.append([s,e,d,0,-1]); continue
        P=acc+"|"+repl+"|"+r["element"]; fam=fid(m2rep.get(P))
        out.append([s,e,d,1,fam])
        if fam!=-1:
            if P in dset: famfn[fam]=[3,dset[P]]
            elif fam not in famfn:
                doms=pfam.get(m2rep.get(P),[])
                if doms:
                    code=1; lab=doms[0][0]
                    for dom,desc in doms:
                        c=cls(dom,desc)
                        if c in (2,4): code=c; lab=dom; break
                    famfn[fam]=[code,lab]
                else: famfn[fam]=[0,"no annotation"]
    out.sort(key=lambda z:z[0])
    genomes.append({"name":name,"acc":acc,"f":out})
    print(acc, sum(1 for z in out if z[3]==1),"cassettes")
data={"genomes":genomes,"fam":{str(k):v for k,v in famfn.items()}}
json.dump(data,open(PROJ+"/data/example/comparison.json","w"),separators=(",",":"))
for i in range(len(genomes)-1):
    A=set(z[4] for z in genomes[i]["f"] if z[3]==1 and z[4]>=0)
    B=set(z[4] for z in genomes[i+1]["f"] if z[3]==1 and z[4]>=0)
    print("shared families:",genomes[i]["acc"],"vs",genomes[i+1]["acc"],"=",len(A&B))
print("wrote comparison.json",os.path.getsize(PROJ+"/data/example/comparison.json"),"bytes")
