#!/bin/bash
#SBATCH --job-name=if-scale
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --array=1-988%50
#SBATCH --output=logs/if_%A_%a.out
#SBATCH --error=logs/if_%A_%a.err
set -euo pipefail
export PATH=/scratch/user/uqcngu19/micromamba_root/envs/integronfinder/bin:$PATH
PROJ=/scratch/user/uqcngu19/vibrio-integron-pipeline
cd $PROJ
G=$(sed -n "${SLURM_ARRAY_TASK_ID}p" genome_list.txt)
[[ -z "${G:-}" ]] && { echo "past end of list"; exit 0; }
N=$(basename "$G" .fna)
OUT=$PROJ/results/scale/$N
compgen -G "$OUT/Results_Integron_Finder_*/*.integrons" >/dev/null && { echo "skip $N"; exit 0; }
S=$(date +%s)
integron_finder "$G" --outdir "$OUT" --local-max --circ --func-annot \
  --distance-thresh 10000 --keep-palindromes --cpu ${SLURM_CPUS_PER_TASK} --mute
E=$(date +%s)
MB=$(du -sm "$OUT" | cut -f1)
AT=$(awk -F'\t' 'NR>2 && tolower($8)=="attc"' $OUT/Results_Integron_Finder_*/*.integrons 2>/dev/null | wc -l)
echo "TIMING $N seconds=$((E-S)) megabytes=$MB attC=$AT"
