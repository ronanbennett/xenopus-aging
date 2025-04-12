#!/usr/bin/bash
# run this program with qsub.

#$ -V
#$ -l h_rt=7:00:00,h_data=2G
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

# echo job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo " "

echo 'running trim fastq...'

fastp --in1 ${DATADIR}/${FROGNAME}_R1.fastq.gz \
    --in2 ${DATADIR}/${FROGNAME}_R2.fastq.gz \
    --out1 ${OUTDIR}/${FROGNAME}_R1_trimmed.fastq.gz --out2 ${OUTDIR}/${FROGNAME}_R2_trimmed.fastq.gz \
    --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    --thread 8 --length_required 40 --trim_front1 5 --trim_tail1 5 \
    --trim_front2 5 --trim_tail2 5 --trim_poly_g --qualified_quality_phred 20 \
    --html ${LOGDIR}/${FROGNAME}_report.html --json ${LOGDIR}/${FROGNAME}_report.json \
    1> ${LOGDIR}/${FROGNAME}.trimmed.stdout

echo 'running alignment to reference index...'
python -m bsbolt Align \
    -DB ${INDEXDIR} -F1 ${OUTDIR}/${FROGNAME}_R1_trimmed.fastq.gz -F2 ${OUTDIR}/${FROGNAME}_R2_trimmed.fastq.gz \
    -t ${THREADS} -OT ${THREADS} \
    -O ${OUTDIR}/${FROGNAME} \
    1> ${LOGDIR}/${FROGNAME}.bam.stdout
rm ${OUTDIR}/${FROGNAME}_R{1,2}_trimmed.fastq.gz

echo 'running fixmate bam...'
samtools fixmate -p -m --threads ${THREADS} \
    ${OUTDIR}/${FROGNAME}.bam ${OUTDIR}/${FROGNAME}.fixmate.bam \
    1> ${LOGDIR}/${FROGNAME}.fixmate.bam.stdout
rm ${OUTDIR}/${FROGNAME}.bam

echo 'running sort bam...'
samtools sort -@ ${THREADS} ${OUTDIR}/${FROGNAME}.fixmate.bam -o ${OUTDIR}/${FROGNAME}.sorted.fixmate.bam \
    1> ${LOGDIR}/${FROGNAME}.sorted.fixmate.bam.stdout
rm ${OUTDIR}/${FROGNAME}.fixmate.bam

echo 'running dedup bam...'
samtools markdup -r ${OUTDIR}/${FROGNAME}.sorted.fixmate.bam ${OUTDIR}/${FROGNAME}.dedup.sorted.fixmate.bam \
    1> ${LOGDIR}/${FROGNAME}.dedup.sorted.fixmate.bam.stdout
rm ${OUTDIR}/${FROGNAME}.sorted.fixmate.bam

echo 'running index bam...'  # creates .bam.bai
samtools index -@ ${THREADS} ${OUTDIR}/${FROGNAME}.dedup.sorted.fixmate.bam \
    1> ${LOGDIR}/${FROGNAME}.dedup.sorted.fixmate.bam.bai.stdout

echo 'running CallMethylation...'
python -m bsbolt CallMethylation -t ${THREADS} -BQ 10 -MQ 20 -IO -CG \
    -I ${OUTDIR}/${FROGNAME}.dedup.sorted.fixmate.bam -DB ${INDEXDIR} \
    -O ${OUTDIR}/${FROGNAME} \
    1> ${LOGDIR}/${FROGNAME}.meth.stdout

# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "

