#!/bin/bash
#$ -S /bin/bash
#$ -cwd

module load python/2.7.3

python $*
