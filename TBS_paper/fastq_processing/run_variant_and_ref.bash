
filename="sample_names.txt"
for id in $(cat "$filename"); do
    qsub -N "variant_$id" variant_and_ref.bash "$id"
done
