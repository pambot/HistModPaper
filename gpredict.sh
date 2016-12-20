#!/bin/bash
#$ -S /bin/bash
#$ -cwd


cd $HOME/gpredict
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
gunzip gencode.v19.annotation.gtf.gz

# segment gencode
awk -F '\t' '{if ($3 == "transcript") print $0}' gencode.v19.annotation.gtf | grep 'transcript_type "protein_coding"' | awk '$5-$4 > 200' | awk -F"\t" 'BEGIN{OFS="\t";}{print $1,$4,$5,$9,1,$7;}' > coding.bed

awk -F '\t' '{if ($3 == "transcript") print $0}' gencode.v19.long_noncoding_RNAs.gtf | awk '$5-$4 > 200' | awk -F"\t" 'BEGIN{OFS="\t";}{print $1,$4,$5,$9,1,$7;}'> lncRNA.bed

# get the TSS and TTS info
# eliminate overlapping TSS by +/- 4kbp
module load bedtools

function elimOverlapsMakeBed() {
	python get_intervals.py coding.bed $1
	python get_intervals.py lncRNA.bed 0
	bedtools intersect -v -a coding_tss.bed -b lncRNA_tss.bed > coding_tss_filt.bed
	bedtools intersect -v -a coding_tts.bed -b lncRNA_tts.bed > coding_tts_filt.bed
	echo >&2 'Done coding filtering.'
	python get_intervals.py coding.bed 0
	python get_intervals.py lncRNA.bed $1
	bedtools intersect -v -a lncRNA_tss.bed -b coding_tss.bed > lncRNA_tss_filt.bed
	bedtools intersect -v -a lncRNA_tts.bed -b coding_tts.bed > lncRNA_tts_filt.bed
	echo >&2 'Done lncRNA filtering.'
	python combine_filt.py coding
	python combine_filt.py lncRNA
	rm -f *_t{t,s}s_filt.bed
	sed -i '/^\s*$/d' *_t{t,s}s.bed
	echo >&2 'Done combining and cleaning. Run complete.'
	count=$(wc -l *_t{t,s}s.bed)
	echo >&2 ${count[*]} 
}

elimOverlapsMakeBed 4000
mv lncRNA_tss.bed lncRNA_g1.bed
mv coding_tss.bed coding_g1.bed
mv lncRNA_tts.bed lncRNA_g2.bed
mv coding_tts.bed coding_g2.bed

## get plus/minus of intervals
# if plus strand, TSS=start and TTS=end; if minus, vice versa
for bed in lncRNA_g1.bed coding_g1.bed lncRNA_g2.bed coding_g2.bed ; do
python recalibrate_intervals.py $bed 1000 pm1K
echo "Done recalibrating $bed"
done

cat CodingTSSpm1K.bed LncrnaTSSpm1K.bed > TSSpm1K.bed
cat CodingTTSpm1K.bed LncrnaTTSpm1K.bed > TTSpm1K.bed

# before this, download ChromHMM.zip and unzip
mv ChromHMM/COORDS/hg19 ChromHMM/COORDS/hg19_orig
mkdir ChromHMM/COORDS/hg19

for bed in *pm1K.bed ; do
mv $bed ChromHMM/COORDS/hg19/${bed%pm1K.bed}.bed
done

# do same but with plus/minus 50 for locus expression at TSS
for bed in lncRNA_g1.bed coding_g1.bed lncRNA_g2.bed coding_g2.bed ; do
python recalibrate_intervals.py $bed 50 pm50
echo "Done recalibrating $bed"
done

mkdir cage
for bed in *pm50.bed ; do
mv $bed cage/$bed
done

# set some stuff up
for cell in Gm12878 Helas3 Hepg2 H1hesc Huvec K562 Nhek ; do
mkdir -p preprocess/$cell
done


## process fastq files for all cell lines
# alignment of reads
#screen
#qlogin
cd $HOME/gpredict
ENCODE=http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeBroadHistone

rm -f files.txt
wget $ENCODE/files.txt

rm -f wgetFiles.txt
for cell in Gm12878 Helas3 Hepg2 H1hesc Huvec K562 Nhek ; do
for sig in "Dnase" "H2az" "H3k0\?9me3" "H3k27ac" "H3k27me3" "H3k36me3" "H3k0\?4me1" "H3k0\?4me2" "H3k0\?4me3" "H3k79me2" "H3k0\?9ac" "H4k20me1" "Control" ; do
awk 'BEGIN{OFS=FS="\t"}{print $1}' files.txt | grep "fastq" | grep $cell | grep $sig  >> wgetFiles.txt
done
done
sed -i 's%^%'$ENCODE'/%' wgetFiles.txt

wget -i wgetFiles.txt -P preprocess/
# cntl+a+d to get out of qlogin

# wait until wget is finished
for cell in Gm12878 Helas3 Hepg2 H1hesc Huvec K562 Nhek ; do
mv preprocess/*$cell*.fastq.gz preprocess/$cell
done

# run bwa for alignment
for cell in Gm12878 Helas3 Hepg2 H1hesc Huvec K562 Nhek ; do
for fqgz in preprocess/$cell/*.fastq.gz ; do
qsub bwa.sh $fqgz
done
done

# convert BAM to BED files
for cell in Gm12878 Helas3 Hepg2 H1hesc Huvec K562 Nhek ; do 
for bam in preprocess/$cell/*.bam ; do 
qsub bam2bed.sh $bam
done
done

# delete to save space
for cell in Gm12878 Helas3 Hepg2 H1hesc Huvec K562 Nhek ; do
rm -f preprocess/$cell/*.bam
done


## get histone bigwig info
rm -f wgetBigWig.txt
for cell in Gm12878 Helas3 Hepg2 H1hesc Huvec K562 Nhek ; do
for sig in "H2az" "H3k0\?9me3" "H3k27ac" "H3k27me3" "H3k36me3" "H3k0\?4me1" "H3k0\?4me2" "H3k0\?4me3" "H3k79me2" "H3k0\?9ac" "H4k20me1" ; do
awk 'BEGIN{OFS=FS="\t"}{print $1}' files.txt | grep "bigWig" | grep $cell | grep $sig  >> wgetBigWig.txt
done
done
sed -i 's%^%'$ENCODE'/%' wgetBigWig.txt

wget -i wgetBigWig.txt -P preprocess/
# cntl+a+d to get out of qlogin

# wait until wget is finished
for cell in Gm12878 Helas3 Hepg2 H1hesc Huvec K562 Nhek ; do
mv preprocess/*$cell*.bigWig preprocess/$cell
done

# run bigwigsummary
mkdir -p predict3/summary

for cell in Gm12878 Helas3 Hepg2 H1hesc Huvec K562 Nhek ; do
for histone in $(ls preprocess/$cell/*.bigWig) ; do
qsub bigWigSummary.sh $histone TSSpm1K.bed
qsub bigWigSummary.sh $histone TTSpm1K.bed
done
done

## switch to cage to calculate gene expression information
cd $HOME/gpredict/cage

# now switch to predict.sh in the same directory


