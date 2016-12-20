#!/bin/bash
#$ -S /bin/bash
#$ -cwd

module load samtools/0.1.19

for bam in bam/*.bam ; do
	ncounts=$(samtools view -c -F 4 $bam)
	echo -e $(basename $bam)"\t"$ncounts >> samcounts.txt
done
