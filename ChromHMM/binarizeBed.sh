#!/bin/bash
#$ -S /bin/bash
#$ -cwd

cd $HOME/gpredict/ChromHMM
module load java/1.6

java -jar ChromHMM.jar BinarizeBed \
CHROMSIZES/hg19.txt \
bedFiles-2 \
Hist-2.files \
Hist-2_bin
