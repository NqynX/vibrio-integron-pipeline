#!/bin/bash
#SBATCH --job-name=vibrio-pipeline
#SBATCH --cpus-per-task=16
export PATH=/scratch/user/uqcngu19/micromamba_root/envs/integronfinder/bin:/scratch/user/uqcngu19/micromamba_root/envs/rgi/bin:$PATH
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/pipeline_%j.out
#SBATCH --error=logs/pipeline_%j.err

INTEGRON_FINDER=/scratch/user/uqcngu19/micromamba_root/envs/integronfinder/bin/integron_finder
RGI=/scratch/user/uqcngu19/micromamba_root/envs/rgi/bin/rgi
BASEDIR=/scratch/user/uqcngu19/vibrio-integron-pipeline
THREADS=16

module load blast+/2.14.1-gompi-2023a

for PAIR in \
  "data/reference_genomes/V_vulnificus_CMCP6/V_vulnificus_CMCP6.fna V_vulnificus_CMCP6" \
  "data/reference_genomes/V_parahaemolyticus_RIMD2210633/V_parahaemolyticus_RIMD2210633.fna V_parahaemolyticus_RIMD2210633" \
  "data/reference_genomes/V_cholerae_O395/V_cholerae_O395.fna V_cholerae_O395" \
  "data/reference_genomes/V_fischeri_ES114/V_fischeri_ES114.fna V_fischeri_ES114"; do

  GENOME=$(realpath $(echo $PAIR | awk '{print $1}'))
  NAME=$(echo $PAIR | awk '{print $2}')
  echo "=== Pipeline for $NAME ==="

  echo "Step 1: IntegronFinder..."
  OUTDIR=$BASEDIR/results/integron_detection/integronfinder/${NAME}
  mkdir -p $OUTDIR
  $INTEGRON_FINDER $GENOME --outdir $OUTDIR 2>&1 | tail -3

  echo "Step 2: Building cassette table..."
  INTEG_FILE=$(ls $OUTDIR/Results_Integron_Finder_*/*.integrons 2>/dev/null | head -1)
  if [ -z "$INTEG_FILE" ]; then echo "No integrons — skipping $NAME"; echo ""; continue; fi
  mkdir -p $BASEDIR/results/cassette_reconstruction
  python3 $BASEDIR/scripts/cassette_reconstruction/build_cassette_table.py $INTEG_FILE $BASEDIR/results/cassette_reconstruction/${NAME}_cassettes.tsv
  CASSETTE_COUNT=$(tail -n +2 $BASEDIR/results/cassette_reconstruction/${NAME}_cassettes.tsv | wc -l)
  echo "Cassettes: $CASSETTE_COUNT"
  if [ "$CASSETTE_COUNT" -eq 0 ]; then echo "No cassettes — skipping $NAME"; echo ""; continue; fi

  echo "Step 3: Extracting cassette sequences..."
  python3 $BASEDIR/scripts/cassette_reconstruction/extract_cassette_seqs.py --cassettes $BASEDIR/results/cassette_reconstruction/${NAME}_cassettes.tsv --genome $GENOME --combined $BASEDIR/results/cassette_reconstruction/${NAME}_cassettes.fasta

  echo "Step 4: BLASTx against Swiss-Prot..."
  mkdir -p $BASEDIR/results/annotation
  blastx -query $BASEDIR/results/cassette_reconstruction/${NAME}_cassettes.fasta -db $BASEDIR/data/databases/uniprot_sprot.fasta -out $BASEDIR/results/annotation/${NAME}_cassettes_blastx_swissprot.tsv -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" -evalue 1e-5 -num_threads $THREADS -max_target_seqs 1 2>&1 | tail -1

  echo "Step 5: RGI..."
  $RGI main -i $BASEDIR/results/cassette_reconstruction/${NAME}_cassettes.fasta -o $BASEDIR/results/annotation/${NAME}_cassettes_rgi -t contig --clean --local -n $THREADS 2>&1 | tail -3

  echo "=== Done: $NAME ==="; echo ""
done

echo "=== ALL PIPELINES COMPLETE ==="
