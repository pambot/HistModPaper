#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe threaded 25
#$ -M pamela.wu@nyumc.org
#$ -m e

module load java/1.6

java -jar ChromHMM.jar LearnModel \
-p $NSLOTS \
-color 0,70,120 \
-init information \
-nobed \
-nobrowser \
Hist-2_bin Hist-2_init 60 hg19
