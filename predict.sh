#!/bin/bash
#$ -S /bin/bash
#$ -cwd

module load bedtools/2.22.0

# assign list of all genes group categories
awk '{print $4}' predict2/coords/CodingTSS.bed | sed -e 's/^/C/g' > temp1
awk '{print $4}' predict2/coords/LncrnaTSS.bed | sed -e 's/^/L/g' > temp2
cat temp1 temp2 > loci.txt
rm temp*

# now format the data for learning
cd predict2
python process_multicov.py

for cell in Gm12878 Helas3 Hepg2 H1hesc Huvec K562 Nhek ; do
python train_labels.py $cell
done
cd ..

for folder in predict1/ predict3/ predict4/ ; do
cp predict2/train/*.labels $folder/train
done

# regression
for folder in predict1/ predict2/ predict3/ predict4/ ; do
for cell in Gm12878 Helas3 Hepg2 H1hesc Huvec K562 Nhek ; do
qsub python.sh regr_spec.py $folder $cell 0.0
qsub python.sh regr_expr.py $folder $cell 0.0
done
done

for folder in predict1/ predict2/ predict3/ predict4/ ; do
qsub python.sh regr_expr_cellvscell.py $folder 0
qsub python.sh regr_expr_cellvscell.py $folder 1
qsub python.sh regr_spec_cellvscell.py $folder 0
qsub python.sh regr_spec_cellvscell.py $folder 1
done

for folder in predict1/ predict2/ predict3/ predict4/ ; do
python plot_cellvscell.py $folder rgrExprValid_full.heatmap
python plot_cellvscell.py $folder rgrExprValid_res.heatmap
python plot_cellvscell.py $folder rgrSpecValid_full.heatmap
python plot_cellvscell.py $folder rgrSpecValid_res.heatmap
done

# predict classes by expression thesholds: S1=full, S2=expressed
for folder in predict1/ predict2/ predict3/ predict4/ ; do
for cell in Gm12878 Helas3 Hepg2 H1hesc Huvec K562 Nhek ; do
qsub python.sh predict.py $folder $cell Coding S2 0.0 0.1
qsub python.sh predict.py $folder $cell Coding S2 0.1 2.0
qsub python.sh predict.py $folder $cell Coding S2 2.0 4.0
qsub python.sh predict.py $folder $cell Coding S2 4.0 6.0
qsub python.sh predict.py $folder $cell Coding S2 6.0 20.0 # catch all leftover
done
done

# predict in general
for folder in predict1/ predict2/ predict3/ predict4/ ; do
for cell in Gm12878 Helas3 Hepg2 H1hesc Huvec K562 Nhek ; do
qsub python.sh predict.py $folder $cell Coding S2 0.0 20.0
done
done

# stability selection for features
for folder in predict1/ predict2/ predict3/ predict4/ ; do
for cell in Gm12878 Helas3 Hepg2 H1hesc Huvec K562 Nhek ; do
qsub python.sh feature_select_lin_coding.py $folder $cell
done
done

for folder in predict1/ predict2/ predict3/ predict4/ ; do
for cell in Gm12878 Helas3 Hepg2 H1hesc Huvec K562 Nhek ; do
qsub python.sh feature_select_lin_expr.py $folder $cell full
qsub python.sh feature_select_lin_expr.py $folder $cell res
done
done

for folder in predict1/ predict2/ predict3/ predict4/ ; do
for cell in Gm12878 Helas3 Hepg2 H1hesc Huvec K562 Nhek ; do
qsub python.sh feature_select_lin_spec.py $folder $cell full
qsub python.sh feature_select_lin_spec.py $folder $cell res
done
done

# feature selection for mlp
n_layers=1
for folder in predict1/ predict2/ predict3/ predict4/ ; do
for cell in Gm12878 Helas3 Hepg2 H1hesc Huvec K562 Nhek ; do
for seed in $(seq 0 20) ; do
qsub python.sh feature_select_mlp_coding.py $folder $cell $n_layers $seed
done
done
done

n_layers=1
for folder in predict1/ predict2/ predict3/ predict4/ ; do
for cell in Gm12878 Helas3 Hepg2 H1hesc Huvec K562 Nhek ; do
for seed in $(seq 0 20) ; do
qsub python.sh feature_select_mlp_expr.py $folder $cell full $n_layers $seed
qsub python.sh feature_select_mlp_expr.py $folder $cell res $n_layers $seed
done
done
done

n_layers=1
for folder in predict1/ predict2/ predict3/ predict4/ ; do
for cell in Gm12878 Helas3 Hepg2 H1hesc Huvec K562 Nhek ; do
for seed in $(seq 0 20) ; do
qsub python.sh feature_select_mlp_spec.py $folder $cell full $n_layers $seed
qsub python.sh feature_select_mlp_spec.py $folder $cell res $n_layers $seed
done
done
done

# get significant coefs
for folder in predict1/ predict2/ predict3/ predict4/ ; do
for cell in Gm12878 Helas3 Hepg2 H1hesc Huvec K562 Nhek ; do
for condition in Coding ; do
qsub python.sh process_mlp_perm.py $folder $cell $condition
done
done
done

for folder in predict1/ predict2/ predict3/ predict4/ ; do
for cell in Gm12878 Helas3 Hepg2 H1hesc Huvec K562 Nhek ; do
for condition in Multicov Specificity ; do
for ttype in full res ; do
qsub python.sh process_mlp_perm.py $folder $cell $condition $ttype ;
done
done
done
done

# predict for cells vs cells
python predict_cellvscell.py Coding S2 fscore

for folder in predict1/ predict2/ predict3/ predict4/ ; do
python plot_cellvscell.py $folder cellvsMeanCodingFscore.heatmap
python plot_cellvscell.py $folder rgrExprValid_full.heatmap
python plot_cellvscell.py $folder rgrExprValid_res.heatmap
python plot_cellvscell.py $folder rgrSpecValid_full.heatmap
python plot_cellvscell.py $folder rgrSpecValid_res.heatmap
done

# run statistical tests and plotting
python plot_byexpr.py Coding S2
python plot_cellvscell_compare.py > compare_cellvscell.txt
python plot_feature_select.py > compare_feature_select.txt
python plot_lin_compare.py > compare_lin.txt
python plot_mlp_compare.py > compare_mlp.txt

# mlp time
for folder in predict1/ predict2/ predict3/ predict4/ ; do
for cell in Gm12878 Helas3 Hepg2 H1hesc Huvec K562 Nhek ; do
for n_layers in $(seq 1 5); do
qsub python.sh predict_mlp.py $folder $cell $n_layers
done
done
done

for folder in predict1/ predict2/ predict3/ predict4/ ; do
for cell in Gm12878 Helas3 Hepg2 H1hesc Huvec K562 Nhek ; do
for n_layers in $(seq 1 5); do
qsub python.sh regr_expr_mlp.py $folder $cell full 0.0 $n_layers
qsub python.sh regr_expr_mlp.py $folder $cell res 0.0 $n_layers
done
done
done

for folder in predict1/ predict2/ predict3/ predict4/ ; do
for cell in Gm12878 Helas3 Hepg2 H1hesc Huvec K562 Nhek ; do
for n_layers in $(seq 1 5); do
qsub python.sh regr_spec_mlp.py $folder $cell full 0.0 $n_layers
qsub python.sh regr_spec_mlp.py $folder $cell res 0.0 $n_layers
done
done
done





