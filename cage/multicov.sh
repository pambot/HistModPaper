#!/bin/bash
#$ -S /bin/bash
#$ -cwd

module load bedtools/2.17.0

bname=$(basename ${1%.*})
bedtools multicov -s -bams $1 -bed $2 > multiraw/$bname.multicov

