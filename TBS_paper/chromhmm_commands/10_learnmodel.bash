#!/bin/bash

# $1 is directory containing binarized input files (assumed to contain '_binary')
# $2 is output directory
# $3 is number of hidden states to learn

java -mx4000M -jar ChromHMM.jar LearnModel \
    ${1} \
    ${2} \
    ${3} \
    xentro10.0
