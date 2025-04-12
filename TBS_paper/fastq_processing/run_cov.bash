for file in ${SCRATCH}/outputs/*.bam; do
    qsub -V cov.bash ${file} 
done
