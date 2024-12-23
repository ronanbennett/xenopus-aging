# $1 is genome .fasta filename
# this will put index information into the same directory as the $1 fname
module load bowtie2
bowtie2-build XENTR_10.0_genome.fasta XENTR_10.0_genome.fasta
