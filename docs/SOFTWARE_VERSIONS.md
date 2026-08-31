# Software and database versions
Captured: 2026-08-31T11:27:43Z

## Genome set
Source: NCBI RefSeq, taxon Vibrionaceae, --assembly-level complete --assembly-source RefSeq
Retrieved: 2026-08-31   Genomes: 988
datasets CLI: datasets version: 18.36.0

## Detection
hmmsearch not found in Path please use --hmmsearch option to specify it.
cmsearch not found in Path please use --cmsearch option to specify it.
prodigal not found in Path please use --prodigal option to specify it.
Parameters: --local-max --circ --func-annot --distance-thresh 10000 --keep-palindromes --cpu 4

## HattCI
HattCI 1.0b, built from https://github.com/maribuon/HattCI
commit: 653f632
Parameters: -b (both strands); all reported hits used (no Vscore filter)

## Clustering
MMseqs2 18.8cc5c
Parameters: easy-cluster --min-seq-id 0.5 -c 0.8 --cov-mode 0

## Annotation
DefenseFinder Usage: defense-finder [OPTIONS] COMMAND [ARGS]... Try 'defense-finder -h' for help. 
models: defense-finder-models 3.1.0, CasFinder 3.1.0
HMMER # HMMER 3.4 (Aug 2023); http://hmmer.org/
Pfam-A: HMMER3/f [3.3 | Nov 2019]
AMR: NCBIfam-AMRFinder HMMs bundled with IntegronFinder --func-annot

## Distance
Mash 2.3
Parameters: sketch -s 10000; dist all-vs-all
