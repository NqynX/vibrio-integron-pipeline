#!/bin/bash
#SBATCH --job-name=df-array
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --array=1-911%40
#SBATCH --output=logs/df_%A_%a.out
#SBATCH --error=logs/df_%A_%a.err
export PATH=/scratch/user/uqcngu19/micromamba_root/envs/vibrio-tools/bin:$PATH
export TMPDIR=/scratch/user/uqcngu19/tmp
mkdir -p $TMPDIR
PROJ=/scratch/user/uqcngu19/vibrio-integron-pipeline
cd $PROJ
F=$(sed -n "${SLURM_ARRAY_TASK_ID}p" defense_list.txt)
if [ -z "$F" ]; then echo "past end"; exit 0; fi
ACC=$(basename "$F" .faa)
OUT=$PROJ/analysis/defense_per_genome/$ACC
if [ -s "$OUT/${ACC}_defense_finder_systems.tsv" ]; then echo "skip $ACC"; exit 0; fi
mkdir -p "$OUT"
defense-finder run "$F" -o "$OUT" --workers ${SLURM_CPUS_PER_TASK}
echo "done $ACC rc=$?"
