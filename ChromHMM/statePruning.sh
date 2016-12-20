#!/bin/bash
#$ -S /bin/bash
#$ -cwd

module load java/1.6

java -jar ChromHMM.jar StatePruning Hist-2_init Hist-2_elim
