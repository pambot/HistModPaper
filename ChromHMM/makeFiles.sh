#!/bin/bash
#$ -S /bin/bash
#$ -cwd

cd $HOME/gpredict/ChromHMM
rm -f Hist-2.files

for cell in Gm12878 H1hesc Helas3 Hepg2 Huvec K562 Nhek ; do
for sig in "H2az" "H3k0\?9me3" "H3k27ac" "H3k27me3" "H3k36me3" "H3k0\?4me1" "H3k0\?4me2" "H3k0\?4me3" "H3k79me2" "H3k0\?9ac" "H4k20me1" ; do
rsig=${sig/0\\\?/}
bed=$(ls bedFiles-2/* | grep $cell | grep $sig | xargs basename)
cntl=$(ls bedFiles-2/* | grep $cell | grep Control | xargs basename)
echo -e "$cell\t$rsig\t$bed\t$cntl" >> Hist-2.files
done
done
