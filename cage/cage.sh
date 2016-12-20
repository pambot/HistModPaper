#!/bin/bash
#$ -S /bin/bash
#$ -cwd

# wget bam and bai files
cd $HOME/gpredict/cage2
encode=ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRikenCage/
wget "${encode}files.txt"
awk 'BEGIN{FS="\t"}{print $1}' files.txt | grep '.*.ba[mi]' > bamfiles.txt
while read line ; do echo $encode$line >> temp ; done < bamfiles.txt
mv -f temp bamfiles.txt
rm -f temp

mkdir -p bam
qsub wget.sh -i bamfiles.txt -P bam

# get samtools reads for all bam
qsub samcounts.sh

# populate array of cell line names and get multicovs
cat CodingTSSpm50.bed LncrnaTSSpm50.bed > TSSpm50.bed

mkdir -p multiraw
for bam in bam/*.bam ; do
	qsub multicov.sh $bam TSSpm50.bed
done

# when multicovs are processed
mkdir -p multicov

cells=$(ls multiraw/*.multicov | sed -e "s/multiraw\///g;s/wg.*Cage//g;s/Cell.*//g;s/Nuc.*//g;s/Poly.*//g;s/Chrom.*//g;s/Cyto.*//g" | uniq)

# combine by cell line and log2 transform
for cell in $cells ; do
	qsub python.sh cell_multicov.py $cell
done

# clean up large files
rm -rf bam

# move multicov files
cp multicov ../predict2

# switch to predict2






