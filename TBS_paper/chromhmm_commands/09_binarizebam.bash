#!/bin/bash

# $1 is dir containing input bams
# $2 is cellmarkfiletable
# $3 is output dir
# must make cellmarkfiletable before this
java -mx4000M -jar ChromHMM/ChromHMM.jar BinarizeBam \
    xentro_10.0_chromsizes.tsv \
    ${1} \  # dir containing input .bams
    ${2} \  # 
    ${3} # output dir for binarized data files
