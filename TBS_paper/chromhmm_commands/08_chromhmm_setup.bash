#!/bin/bash

# NOTE the chromosome names in chromsizes and ncbiRefSeq had to be manually changed from chr to Chr

java -mx4000M -jar ChromHMM/ChromHMM.jar ConvertGeneTable -l xentro_10.0_chromsizes.tsv -noheader ncbiRefSeq.txt xentro10.0 xentro10.0
