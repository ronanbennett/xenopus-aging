#!/usr/bin/bash
# run this program with qsub. No additional args needed.

#$ -V
#$ -l h_rt=2:00:00,h_data=2G
#$ -pe shared 10

# stop execution if any step fails.
set -e 
set -o pipefail

FROGNAME="$1"
DATADIR="${SCRATCH}/data"
OUTDIR="${SCRATCH}/outputs"
INDEXDIR="${SCRATCH}/references"
LOGDIR="${SCRATCH}/logs"
THREADS=10

python -m bsbolt CallVariation -OR -I ${OUTDIR}/${FROGNAME}.dedup.sorted.fixmate.bam -DB ${SCRATCH}/references -BQ 10 -O ${OUTDIR}/${FROGNAME}.allgeno
