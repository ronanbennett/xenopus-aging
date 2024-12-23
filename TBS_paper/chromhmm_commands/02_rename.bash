# $1 is the folder (w/o trailing slash) that stores all the files of the format fastq?acc=*
for file in ${1}/fastq?acc=*; do
    # Extract the part after '=' and before the space (if any)
    newname=$(echo ${file} | sed -e 's/fastq?acc=\(.*\)/\1.fastq.gz/')

    # Rename the file
    mv ${file} ${newname}
done
