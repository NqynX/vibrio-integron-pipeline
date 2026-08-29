#!/bin/bash
# Usage: bash scripts/run_pipeline.sh <genome.fna> <sample_name>
GENOME=$1; NAME=$2; BASEDIR=/scratch/user/uqcngu19/vibrio-integron-pipeline; THREADS=8

echo "=== Pipeline for $NAME ==="
echo "Step 1: Running IntegronFinder..."
eval "$(micromamba shell hook --shell bash)"
micromamba activate integronfinder
OUTDIR=$BASEDIR/results/integron_detection/integronfinder/${NAME}
integron_finder $GENOME --outdir $OUTDIR 2>&1 | tail -3

echo "Step 2: Building cassette table..."
INTEG_FILE=$(ls $OUTDIR/Results_Integron_Finder_*/${NAME}*.integrons 2>/dev/null | head -1)
if [ -z "$INTEG_FILE" ]; then echo "ERROR: No .integrons file found"; exit 1; fi
python3 $BASEDIR/scripts/cassette_reconstruction/build_cassette_table.py $INTEG_FILE $BASEDIR/results/cassette_reconstruction/${NAME}_cassettes.tsv
echo "Cassettes found: $(tail -n +2 $BASEDIR/results/cassette_reconstruction/${NAME}_cassettes.tsv | wc -l)"

echo "Step 3: Extracting cassette sequences..."
python3 $BASEDIR/scripts/cassette_reconstruction/extract_cassette_seqs.py $BASEDIR/results/cassette_reconstruction/${NAME}_cassettes.tsv $GENOME $BASEDIR/results/cassette_reconstruction/${NAME}_cassettes.fasta

echo "Step 4: BLASTx against Swiss-Prot..."
module load blast+/2.14.1-gompi-2023a
blastx -query $BASEDIR/results/cassette_reconstruction/${NAME}_cassettes.fasta -db $BASEDIR/data/databases/uniprot_sprot.fasta -out $BASEDIR/results/annotation/${NAME}_cassettes_blastx_swissprot.tsv -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" -evalue 1e-5 -num_threads $THREADS -max_target_seqs 1 2>&1 | tail -1

echo "Step 5: Running RGI..."
micromamba activate rgi
rgi main -i $BASEDIR/results/cassette_reconstruction/${NAME}_cassettes.fasta -o $BASEDIR/results/annotation/${NAME}_cassettes_rgi -t contig --clean --local -n $THREADS 2>&1 | tail -3

echo "=== Pipeline complete for $NAME ==="
echo "Results:"
echo "  IntegronFinder: $OUTDIR/"
echo "  Cassettes:      $BASEDIR/results/cassette_reconstruction/${NAME}_cassettes.tsv"
echo "  Sequences:      $BASEDIR/results/cassette_reconstruction/${NAME}_cassettes.fasta"
echo "  BLASTx:         $BASEDIR/results/annotation/${NAME}_cassettes_blastx_swissprot.tsv"
echo "  RGI:            $BASEDIR/results/annotation/${NAME}_cassettes_rgi.txt"
