#!/bin/bash
#$ -S /bin/bash
#$ -cwd

cd $HOME/gpredict/ChromHMM
module load bedtools/2.17.0

f=$1
reps=$(ls bedFiles-2/*.bed | grep $f)
cat $reps > $f.bed
if [ -s $f.bed ] ; then
rm -f $reps
fi

bedtools sort -i $f.bed > $f.sort.bed
if [ -s $f.sort.bed ] ; then
rm -f $f.bed
fi

mv $f.sort.bed $f.bed

