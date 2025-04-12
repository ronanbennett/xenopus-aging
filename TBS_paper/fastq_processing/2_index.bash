#!/usr/bin/bash
# run this program with qsub. No additional args needed.

#$ -N ind
#$ -V
#$ -l h_rt=2:00:00,h_data=8G
#$ -pe shared 10

# stop execution if any step fails.
set -e 
set -o pipefail

FASTA_FILE="${SCRATCH}/references/XENTR_10.0_genome.fasta"
OUTDIR="${SCRATCH}/references/"

python3 -m bsbolt Index -G ${FASTA_FILE} -DB ${OUTDIR}
