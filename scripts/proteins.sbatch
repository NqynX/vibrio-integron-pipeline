#!/bin/bash
#SBATCH --job-name=cassette-fam
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/fam_%j.out
#SBATCH --error=logs/fam_%j.err
export PATH=/scratch/user/uqcngu19/micromamba_root/envs/integronfinder/bin:$PATH
PROJ=/scratch/user/uqcngu19/vibrio-integron-pipeline
cd $PROJ && mkdir -p analysis

python - <<'PYEOF'
import glob, os, re, csv
from Bio import SeqIO
from Bio.Seq import Seq

PROJ = "/scratch/user/uqcngu19/vibrio-integron-pipeline"
GD   = f"{PROJ}/data/complete_vibrionaceae/genomes"
out  = open(f"{PROJ}/analysis/cassette_proteins.faa", "w")
meta = open(f"{PROJ}/analysis/protein_meta.csv", "w", newline="")
w = csv.writer(meta)
w.writerow(["protein_id","accession","replicon","integron","aa_len","is_integrase"])

n = 0
for p in sorted(glob.glob(f"{PROJ}/results/scale/**/*.integrons", recursive=True)):
    acc = re.search(r"results/scale/([^/]+)/", p).group(1)
    fna = f"{GD}/{acc}.fna"
    if not os.path.exists(fna): continue
    hdr, rows = None, []
    for line in open(p):
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if hdr is None: hdr = f; continue
        d = dict(zip(hdr, f))
        if d.get("type_elt","").lower() == "protein": rows.append(d)
    if not rows: continue
    seqs = {r.id: r.seq for r in SeqIO.parse(fna, "fasta")}
    for d in rows:
        s = seqs.get(d["ID_replicon"])
        if s is None: continue
        b, e = int(float(d["pos_beg"])), int(float(d["pos_end"]))
        sub = s[b-1:e]
        if int(float(d.get("strand", 1) or 1)) == -1: sub = sub.reverse_complement()
        prot = str(Seq(str(sub)).translate(to_stop=True))
        if len(prot) < 30: continue
        pid = f"{acc}|{d['ID_replicon']}|{d['element']}"
        out.write(f">{pid}\n{prot}\n"); n += 1
        w.writerow([pid, acc, d["ID_replicon"], d["ID_integron"], len(prot),
                    "inti" in (d.get("annotation","") or "").lower()])
out.close(); meta.close()
print(f"wrote {n} cassette proteins")
PYEOF

export PATH=/scratch/user/uqcngu19/micromamba_root/envs/vibrio-tools/bin:$PATH
export TMPDIR=/scratch/user/uqcngu19/tmp
mkdir -p $TMPDIR
cd $PROJ/analysis && rm -rf clust tmp
mmseqs easy-cluster cassette_proteins.faa clust $TMPDIR/mmseqs_tmp \
  --min-seq-id 0.5 -c 0.8 --threads ${SLURM_CPUS_PER_TASK} > mmseqs.log 2>&1
echo "clustered:"; wc -l clust_cluster.tsv
echo "families:"; cut -f1 clust_cluster.tsv | sort -u | wc -l
