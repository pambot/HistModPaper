#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe threaded 20

module load java/1.6

java -jar ChromHMM.jar LearnModel \
-color 0,70,120 \
-p $NSLOTS \
-init load \
-m Hist-2_elim/elim_${1}_model_60.txt \
-nobrowser \
-nobed \
Hist-2_bin Hist-2_load/Model_$1 $1 hg19


