#!/bin/bash
#SBATCH --job-name=pfam
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=logs/pfam_%j.out
#SBATCH --error=logs/pfam_%j.err
export PATH=/scratch/user/uqcngu19/micromamba_root/envs/vibrio-tools/bin:$PATH
export TMPDIR=/scratch/user/uqcngu19/tmp
PROJ=/scratch/user/uqcngu19/vibrio-integron-pipeline
cd $PROJ/analysis
if [ ! -f Pfam-A.hmm.h3i ]; then
  wget -q https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
  gunzip -f Pfam-A.hmm.gz
  hmmpress Pfam-A.hmm
fi
hmmscan --cut_ga --cpu ${SLURM_CPUS_PER_TASK} \
    --domtblout pfam_reps.tbl Pfam-A.hmm clust_rep_seq.fasta > /dev/null
echo "hits:"; grep -vc '^#' pfam_reps.tbl
echo "reps with a hit:"; grep -v '^#' pfam_reps.tbl | awk '{print $4}' | sort -u | wc -l
