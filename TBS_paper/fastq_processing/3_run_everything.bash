
filename="sample_names.txt"
for id in $(cat "$filename"); do
    qsub -N "methyl_$id" fastq_to_methyl.bash "$id"
done
