#!/bin/bash
#SBATCH --job-name=if-control
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/control_%j.out
#SBATCH --error=logs/control_%j.err
set -euo pipefail
export PATH=/scratch/user/uqcngu19/micromamba_root/envs/integronfinder/bin:$PATH
PROJ=/scratch/user/uqcngu19/vibrio-integron-pipeline
cd $PROJ
GD=$PROJ/data/complete_vibrionaceae/genomes
mkdir -p results/control
for acc in GCF_000006745.1 GCF_000039765.1 GCF_000196095.1 GCF_000016245.1 GCF_000011805.1; do
  F=$(ls $GD/${acc}*.fna 2>/dev/null | head -1)
  [[ -z "$F" ]] && { echo "MISSING $acc"; continue; }
  O=results/control/$acc
  compgen -G "$O/Results_Integron_Finder_*/*.integrons" >/dev/null && { echo "skip $acc"; continue; }
  echo "=== $acc DEFAULT (no --local-max) ==="
  integron_finder "$F" --outdir "$O" --circ --func-annot \
    --distance-thresh 10000 --keep-palindromes --cpu ${SLURM_CPUS_PER_TASK} --mute
done
