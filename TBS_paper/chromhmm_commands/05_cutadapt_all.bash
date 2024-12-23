#! /bin/bash
# $1 is directory of source .fastq.gz files (w/o trailing /) 
# $2 is dir where destination trimmed .fastq.gz will live
for file in ${1}/*.fastq.gz; do
    filename=$(basename ${file})
   bash 05_cutadapt.bash ${1}/${filename} ${2}/${filename}
done
