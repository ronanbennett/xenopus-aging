#!/bin/bash
# $1 is name of directory containing .sam files. No trailing /. Output will be in same directory
module load samtools
for samfile in ${1}/*.sam; do
    filename=$(basename "${samfile}" .sam)
    samtools view -@ 8 -S -b ${samfile} > ${1}/${filename}.bam
    samtools sort -@ 8 ${1}/${filename}.bam -o ${1}/${filename}.sorted.bam
    samtools markdup -@ 8 ${1}/${filename}.sorted.bam ${1}/${filename}.sorted.markedduplicates.bam
    samtools view -@ 8 -b -F 0x400 ${1}/${filename}.sorted.markedduplicates.bam > ${1}/${filename}.sorted.deduplicated.bam
    samtools index -@ 8 ${1}/${filename}.sorted.deduplicated.bam
    # Optional: Remove intermediate unsorted BAM if desired
    # rm "${1}/${filename}.bam"
done
