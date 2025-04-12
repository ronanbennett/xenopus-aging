
gunzip $SCRATCH/outputs/*.CGmap.gz
ls $SCRATCH/outputs/*.CGmap > all_filenames.txt
python -m bsbolt AggregateMatrix -F all_filenames.txt -O $SCRATCH/outputs/methylation_matrix.txt -min-coverage 100 -min-sample 1.0 -CG 
