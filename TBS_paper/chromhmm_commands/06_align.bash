# $1 is directory storing the .fastq.gz files (no trailing /).
# $2 is reference genome .fasta (with indexing information in the same directory).
for file in ${1}/*.fastq.gz; do 
    qsub -N "align_$(basename ${file})" 06_align_template.bash ${file} ${2}
done
