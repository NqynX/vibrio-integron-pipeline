#!/bin/bash
#SBATCH --job-name=hattci
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --array=1-50%25
#SBATCH --output=logs/hat_%A_%a.out
#SBATCH --error=logs/hat_%A_%a.err
PROJ=/scratch/user/uqcngu19/vibrio-integron-pipeline
HAT=$PROJ/tools/HattCI/hattci.out
cd $PROJ
mkdir -p results/hattci
G=$(sed -n "${SLURM_ARRAY_TASK_ID}p" hattci_list.txt)
if [ -z "$G" ]; then echo "past end"; exit 0; fi
N=$(basename "$G" .fna)
O=$PROJ/results/hattci/${N}.hattci.tsv
if [ -s "$O" ]; then echo "skip $N"; exit 0; fi
S=$(date +%s)
$HAT -b -t ${SLURM_CPUS_PER_TASK} "$G" "$O"
E=$(date +%s)
echo "HATTCI $N seconds=$((E-S)) hits=$(($(wc -l < $O)-1)) rc=$?"
