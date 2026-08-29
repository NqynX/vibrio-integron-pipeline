#!/bin/bash
# Usage: bash scripts/run_pipeline.sh <genome.fna> <sample_name>
GENOME=$(realpath $1); NAME=$2; BASEDIR=/scratch/user/uqcngu19/vibrio-integron-pipeline; THREADS=8
INTEGRON_FINDER=/scratch/user/uqcngu19/micromamba_root/envs/integronfinder/bin/integron_finder
RGI=/scratch/user/uqcngu19/micromamba_root/envs/rgi/bin/rgi

echo "=== Pipeline for $NAME ==="

echo "Step 1: Running IntegronFinder..."
OUTDIR=$BASEDIR/results/integron_detection/integronfinder/${NAME}
mkdir -p $OUTDIR
 $INTEGRON_FINDER $GENOME --outdir $OUTDIR 2>&1 | tail -3

echo "Step 2: Building cassette table..."
INTEG_FILE=$(ls $OUTDIR/Results_Integron_Finder_*/*.integrons 2>/dev/null | head -1)
if [ -z "$INTEG_FILE" ]; then echo "No integrons found — skipping remaining steps"; exit 0; fi
mkdir -p $BASEDIR/results/cassette_reconstruction
python3 $BASEDIR/scripts/cassette_reconstruction/build_cassette_table.py $INTEG_FILE $BASEDIR/results/cassette_reconstruction/${NAME}_cassettes.tsv
CASSETTE_COUNT=$(tail -n +2 $BASEDIR/results/cassette_reconstruction/${NAME}_cassettes.tsv | wc -l)
echo "Cassettes found: $CASSETTE_COUNT"
if [ "$CASSETTE_COUNT" -eq 0 ]; then echo "No cassettes — skipping remaining steps"; exit 0; fi

echo "Step 3: Extracting cassette sequences..."
python3 $BASEDIR/scripts/cassette_reconstruction/extract_cassette_seqs.py --cassettes $BASEDIR/results/cassette_reconstruction/${NAME}_cassettes.tsv --genome $GENOME --combined $BASEDIR/results/cassette_reconstruction/${NAME}_cassettes.fasta

echo "Step 4: BLASTx against Swiss-Prot..."
module load blast+/2.14.1-gompi-2023a
mkdir -p $BASEDIR/results/annotation
blastx -query $BASEDIR/results/cassette_reconstruction/${NAME}_cassettes.fasta -db $BASEDIR/data/databases/uniprot_sprot.fasta -out $BASEDIR/results/annotation/${NAME}_cassettes_blastx_swissprot.tsv -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" -evalue 1e-5 -num_threads $THREADS -max_target_seqs 1 2>&1 | tail -1

echo "Step 5: Running RGI..."
 $RGI main -i $BASEDIR/results/cassette_reconstruction/${NAME}_cassettes.fasta -o $BASEDIR/results/annotation/${NAME}_cassettes_rgi -t contig --clean --local -n $THREADS 2>&1 | tail -3

echo "=== Pipeline complete for $NAME ==="
echo "  Cassettes: $(tail -n +2 $BASEDIR/results/cassette_reconstruction/${NAME}_cassettes.tsv | wc -l)"
echo "  BLASTx hits: $(wc -l < $BASEDIR/results/annotation/${NAME}_cassettes_blastx_swissprot.tsv)"
echo "  RGI hits: $(grep -c 'CATALOG' $BASEDIR/results/annotation/${NAME}_cassettes_rgi.txt 2>/dev/null || echo 0)"
