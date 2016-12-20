#!/bin/bash
#$ -S /bin/bash
#$ -cwd

module load java/1.6
cd $HOME/gpredict/ChromHMM

# for cell in Gm12878 Helas3 Hepg2 H1hesc Huvec K562 Nhek ; do for ttype in TSS TTS ; do qsub makeEnrich.sh 20 $cell $ttype ; done ; done

n=$1
cell=$2
ttype=$3

# enrichment
echo -e "RefSeqExon.bed.gz\nRefSeqGene.bed.gz\nCodingTSS.bed\nLncRNATSS.bed" > $cell.$n.$ttype.enrich
for f in $(ls COORDS/hg19/$cell*.bed | grep $ttype | grep Coding | grep Expr) $(ls COORDS/hg19/$cell*.bed | grep TSS | grep Coding | grep Zero) $(ls COORDS/hg19/$cell*.bed | grep $ttype | grep Lncrna | grep Expr) $(ls COORDS/hg19/$cell*.bed | grep TSS | grep Lncrna | grep Zero) ; do
echo $(basename $f) >> $cell.$n.$ttype.enrich
done
java -jar ChromHMM.jar OverlapEnrichment \
-a $cell \
-color 0,70,120 \
-f $cell.$n.$ttype.enrich \
-multicount \
-posterior \
Hist_load/Model_$n/POSTERIOR \
COORDS/hg19 \
Hist_load/Model_$n/${cell}_${ttype}_enrich

rm -f $cell.$n.$ttype.enrich
