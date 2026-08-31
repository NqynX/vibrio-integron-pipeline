#!/bin/bash
#SBATCH --job-name=defense
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=logs/defense_%j.out
#SBATCH --error=logs/defense_%j.err
export PATH=/scratch/user/uqcngu19/micromamba_root/envs/vibrio-tools/bin:$PATH
export TMPDIR=/scratch/user/uqcngu19/tmp
mkdir -p $TMPDIR
PROJ=/scratch/user/uqcngu19/vibrio-integron-pipeline
cd $PROJ/analysis
rm -rf defense
defense-finder run clust_rep_seq.fasta -o defense --workers ${SLURM_CPUS_PER_TASK}
echo "--- systems found ---"
ls defense/
head -3 defense/*systems.tsv 2>/dev/null
