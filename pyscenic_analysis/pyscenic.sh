dir=/data/work 
tfs=$dir/hs_hgnc_tfs.txt
feather=$dir/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather
tbl=$dir/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
input_loom=/data/work/results/CS9/CS9_bin50.loom
ls $tfs  $feather  $tbl  

pyscenic grn \
--num_workers 16 \
--output adj.sample.tsv \
--method grnboost2 \
$input_loom  $tfs 


pyscenic ctx \
adj.sample.tsv $feather \
--annotations_fname $tbl \
--expression_mtx_fname $input_loom  \
--mode "dask_multiprocessing" \
--output reg.csv \
--num_workers 16  \
--mask_dropouts

pyscenic aucell \
$input_loom \
reg.csv \
--output out_SCENIC.loom \
--num_workers 16

