#!/bin/bash
#$ -S /bin/bash
#$ -cwd

# make expression indices for bed files
for cell in Gm12878 Helas3 Hepg2 H1hesc Huvec K562 Nhek ; do
for rtype in Coding Lncrna ; do
python multicov_indices.py $cell $rtype
done
done

for cell in Gm12878 Helas3 Hepg2 H1hesc Huvec K562 Nhek ; do
for rtype in Coding Lncrna ; do
python multicov_unique.py $cell $rtype
done
done

# merge replicates
cd $HOME/gpredict/ChromHMM
files=$(ls bedFiles-2/*.bed | sed 's/Rep.*.bed//g' | uniq)

for f in $files ; do
qsub mergeBedReps.sh $f
done

# make the file list
qsub makeFiles.sh

# bin data for ChromHMM
qsub binarizeBed.sh

# train the initial model
qsub ChromHMM_init.sh

# eliminate one state at a time
qsub statePruning.sh

# load each model
for model in Hist-2_elim/*.txt ; do
n=$(echo $model | cut -d'_' -f 3)
qsub ChromHMM_load.sh $n
done

# transfer to Dropbox
python bic.py

# choose number of states
s=27

# get segmentation files
qsub makeSegments.sh $s

# make enrichments
for cell in Gm12878 Helas3 Hepg2 H1hesc Huvec K562 Nhek ; do 
for ttype in TSS TTS ; do 
qsub makeEnrich.sh $s $cell $ttype
done
done

# send to another directory
mkdir $HOME/gpredict/predict
cp Hist-2_load/Model_$s/* $HOME/gpredict/predict
mv Hist-2_load/Model_$s/POSTERIOR $HOME/gpredict/predict/POSTERIOR
mv Hist-2_load/Model_$s/STATEBYLINE $HOME/gpredict/predict/STATEBYLINE



