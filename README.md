# HistModPaper
This is a collection of the scripts used for the publication, *Comparing models for predicting gene properties using histone modification signals and chromatin states*, that has recently been submitted.

NOTE: [ChromHMM](http://compbio.mit.edu/ChromHMM/), created by Dr. Jason Ernst, was used heavily to compute features, but it's not mine. To get the full experience (not included due to Github space restrictions), download ChromHMM.zip, unzip it, and put all the contents of ChromHMM from this repository into the unzipped folder.

## Abstract
Histone modification signals have been used extensively in models to predict transcriptional regulation. Recently, the idea emerged of using combinations of histone modification signals, called chromatin states, as indicators of regulatory regions. This study aims to examine under what circumstances histone signals values or chromatin states can better predict various gene properties, such as gene expression, coding vs. non-coding loci, and cell-type specificity by constructing features that either encoded the signal values or chromatin states surrounding a gene locus. We compared each feature type's prediction ability and found that expression is best predicted by histone modification signals, whereas cell-type specificity is best predicted by chromatin states. We looked at the biological generalizability of each feature type and how feature types fared under non-linear models, and we also found gene property-specific associations. We also found that chromatin states had the best rank correlation of histone modification importances between linear and non-linear models, and three clusters of histone modification importances could be clearly identified and that chromatin state features were disproportionally kept in all models after feature selection, suggesting that chromatin states have higher information fidelity than histone signal features. Overall, we found that gene properties were consistently better predicted by one feature type over the other, but that chromatin states features have unique properties that were consistent across all models studied.

## Caveat Emptor
The idea behind this repository is to make the code and the materials open for examination and, if it is found satisfactory, to be open for use by other people in the field. As with many scientific code repositories, the scripts have been run on multiple machines (with different compatibility requirements) and have been moved around and sometimes kept without being deleted. In producing this public repository, many things have once again been moved around so that their organization makes sense, and I've tried to modify the scripts so that things are still coherent, but unlike software engineering packages and programs, these scripts were never meant to act as one automatic entity or used in their entirety by another user. These scripts are meant be used as loosely connected discrete entities, as I used them to cobble together my results. It's possible that not all of the scripts work perfectly as is (some of them I wrote years ago), but I'm certain that with minor modifications only, everything from the paper can be reproduced in its entirety. If you experience any problems using these scripts, please email me using the email template: `<firstName>.<lastName>@nyumc.org`.

## How I Ran Scripts
I started with the main script, `gpredict.sh`. Generally, I would run one code block as delineated by a comment header and wait for it to run before going on to the next code block. Whenever you see a line beginning with `qsub`, that meant I ran it on the HPC at NYULMC in order to take advantage of parallel processing. As more subprojects were created, more subdirectories were made with a `<dirname>.sh` as the main script file to roughly preserve the order of execution of the other scripts in the directory.

### Dependencies UPDATE
* Python 2.7+
  * Numpy 1.10.1
  * Scipy 0.16.0
  * Pandas 0.13.1
  * Matplotlib 1.5.0
  * Pybedtools 0.6.8
  * Scikit-Learn 0.18.0
* BWA 0.7.7
* Samtools 0.1.19
* Bedtools 2.17.0
* Java 1.6
* ChromHMM 1.11

### Items Not Included
Luckily, all of the raw data here can be gotten with `wget` and that is all included in the scripts, since all of the data comes from public data hubs like ENCODE and GENCODE. However, biological data files are gigantic, and even including a subset of them would max out this account. That is why I made sure that everything can be either downloaded or derived. I tried to document every single step I used to generate all intermediate files and results files, along with comments to make the code interpretable. 


