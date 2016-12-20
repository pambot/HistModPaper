#!/bin/bash
#$ -S /bin/bash
#$ -cwd

module load java/1.6
cd $HOME/gpredict/ChromHMM

n=$1
cell=$2
ttype=$3

java -jar ChromHMM.jar NeighborhoodEnrichment \
-a $cell \
-color 0,70,120 \
-nostrand \
-posterior \
Hist-2_load/Model_$n/POSTERIOR \

Hist-2_load/Model_$n/${cell}_${ttype}_neigh
