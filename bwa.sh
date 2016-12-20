#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe threaded 20

cd $HOME/gpredict
module load bwa/0.7.7
module load igenomes
module load samtools/0.1.19
REF=$IGENOMES_ROOT/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex
prefix=$(basename ${1%.fastq.gz})

cp $1 $TMPDIR/$prefix.fastq.gz
gunzip $TMPDIR/$prefix.fastq.gz

bwa aln -l 25 -t $NSLOTS $REF/genome.fa $TMPDIR/$prefix.fastq > $TMPDIR/$prefix.sai

bwa samse $REF/genome.fa $TMPDIR/$prefix.sai $TMPDIR/$prefix.fastq | \
samtools view -bS - > $TMPDIR/$prefix.bam

mv $TMPDIR/$prefix.bam ${1%.fastq.gz}.bam
rm -f $TMPDIR/$prefix.*

if [ -s ${1%.fastq.gz}.bam ]; then rm -f $1 ; fi

