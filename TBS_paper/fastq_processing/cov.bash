# Computes coverage

samtools depth ${1} | awk '$3 >= 10 {print $1"\t"$2-1"\t"$2}' > ${1}.10x.bed
sort-bed ${1}.10x.bed > ${1}.10x.sorted.bed
