#3D Hotspot analysis
import os
import glob
import sys
#sys.path.append("/jdfssz1/ST_SUPERCELLS/P18Z10200N0350/USER/haoshijie/software/miniconda3/lib/python3.7/site-packages/")
#sys.path.append('/hwfssz1/ST_SUPERCELLS/P21Z10200N0090/haoshijie/software/local/lib/python3.7/site-packages/')
#sys.path.append("/jdfssz1/ST_SUPERCELLS/P18Z10200N0350/USER/haoshijie/software/miniconda3/lib/python3.7/")
import numpy as np
import pandas as pd
import hotspot
import matplotlib.pyplot as plt
import matplotlib.colors
import seaborn as sns
import pickle
import scanpy as sc
import mplscience

adata = sc.read_h5ad("../input/Files/CS9_CS10_3d_h5ad/CS10_align.h5ad")
metadata = pd.read_csv("CS10_SPATEO_align_metadata.csv",index_col=0)

tmp = adata.obs.index==metadata.index
tmp1 = tmp.tolist()
tmp1.count(True)

adata.obs['x_align'] = metadata['x']
adata.obs['y_align'] = metadata['y']
adata.obs['z_align'] = metadata['z']

adata.layers["counts"] = adata.X.copy()
adata.obsm["spatial"]=adata.obs[["x_align","y_align","z_align"]].to_numpy()
adata.layers["counts_csc"] = adata.layers["counts"].tocsc()

save_path = "./02.HotSpot_Spateo/"
HOTSPOT = ''.join([save_path, "/", "CS10", "_hotspot_spateo.v2.p"])
HS_RESULTS = ''.join([save_path,"/","CS10","_hs_results_spateo.v2.p"])

fl = "CS10"
hs = hotspot.Hotspot(
        adata,
        layer_key="counts_csc",
        model='bernoulli',
        latent_obsm_key="spatial",
        umi_counts_obs_key="nCount_RNA")
hs.create_knn_graph(
weighted_graph=False, n_neighbors=20,)
hs_results = hs.compute_autocorrelations(jobs=16)
#hs_results.head(15)
# Select the genes with significant lineage autocorrelation
hs_genes = hs_results.loc[hs_results.FDR < 0.05].sort_values('Z', ascending=False).head(500).index

# Compute pair-wise local correlations between these genes
lcz = hs.compute_local_correlations(hs_genes, jobs=16)
modules = hs.create_modules(
min_gene_threshold=12, core_only=True, fdr_threshold=0.05)
modules.value_counts()
module_scores = hs.calculate_module_scores()
module_scores.head()
module_cols = []
for c in module_scores.columns:
    key = f"Module {c}"
    adata.obs[key] = module_scores[c]
    module_cols.append(key)
adata.obs.to_csv(save_path+fl+"_bin50-HotSpot_spateo.v2.csv",index=False,sep=',')
with mplscience.style_context():
    sc.pl.spatial(adata, color=module_cols, frameon=False, vmin="p0", vmax="p99", spot_size=1,save = fl+"_bin50-HotSpot_spatial_spateo.v2.png")
plt.rcParams['figure.figsize'] = (15.0, 12.0)
hs.plot_local_correlations()
plt.savefig("".join([save_path,fl,"-regulon_module_number_spateo.v2.pdf"]), dpi = 600)
hs.local_correlation_z.to_csv(save_path+fl+"_bin50-local_correlation_z_spateo.v2.csv",index=True,sep=',')
with open(HOTSPOT, "wb") as f:
    pickle.dump(hs,f,protocol=4)
with open(HS_RESULTS, "wb") as f:
    pickle.dump(hs_results,f)
moduleresult = hs.results.join(hs.modules)
moduleresult.to_csv(save_path+fl+"_bin50-HotSpotmodule_spateo.v2.csv",index=True,sep=',')