# ${1} is filename of list of SRR accession codes, separated by \n
# {2} is location to download the files to 
cat ${1} | parallel -j 10 wget -P ${2} "https://www.be-md.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc={}"
