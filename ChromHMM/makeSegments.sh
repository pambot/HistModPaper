#!/bin/bash
#$ -S /bin/bash
#$ -cwd

module load java/1.6
cd $HOME/gpredict/ChromHMM

n=$1

# segmentation files (required)
java -jar ChromHMM.jar MakeSegmentation \
-printposterior \
-printstatesbyline \
Hist-2_load/Model_$n/model_$n.txt \
Hist-2_bin \
Hist-2_load/Model_$n

