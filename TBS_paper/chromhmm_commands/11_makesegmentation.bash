#!/bin/bash

java -mx4000M -jar ChromHMM/ChromHMM.jar MakeSegmentation \
    histones/CHROMHMM_OUT_10 \
    histones/ \
    histones/segmentation_out
