#!/bin/bash
#$ -S /bin/bash
#$ -cwd

while read -a bed ; do
	ans=$(./bigWigSummary $1 ${bed[0]} ${bed[1]} ${bed[2]} 10)
	if [[ ! $ans =~ [0-9na] ]] ; then
		ans=$(awk 'BEGIN {while (c++<10) printf "n/a "}')
	fi
	echo ${bed[0]} ${bed[1]} ${bed[2]} $ans >> predict3/summary/$(basename ${1%.*}).$(basename ${2%.*}).summary
done < $2
