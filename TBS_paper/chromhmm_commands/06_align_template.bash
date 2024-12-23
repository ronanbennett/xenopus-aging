#!/bin/bash
# $1 is filename .fastq.gz to align
# $2 is reference genome .fasta
# run this with qsub -N name   then the filename and the input param(s)
#$ -V
#$ -l h_rt=2:00:00,h_data=8G
#$ -pe shared 10

module load bowtie2
echo "aligning ${1}"
# arguments from QCBIO 2020 ChIP-seq
srr=$(basename ${1} .fastq.gz)
dir_path=$(dirname ${1})
bowtie2 -q -p 10 -k 1 --no-unal -x ${2} -U ${dir_path}/${srr}.fastq.gz -S ${dir_path}/${srr}.sam
echo "finished aligning"

