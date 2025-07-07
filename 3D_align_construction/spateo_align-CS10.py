import os
os.environ['CUDA_VISIBLE_DEVICES'] = '1'

import torch
device = 'cuda' if torch.cuda.is_available() else 'cpu'
print("Running this notebook on: ", device)

import spateo as st
print("Last run with spateo version:", st.__version__)

# Other imports
import warnings
warnings.filterwarnings('ignore')
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
# Uncomment the following if running on the server
# import pyvista as pv
# pv.start_xvfb()

%config InlineBackend.print_figure_kwargs={'facecolor' : "w"}
%config InlineBackend.figure_format='retina'

# Load the slices
slice1 = st.read('h5ad/CS10/CS10_1_D03251F413.h5ad')
slice2 = st.read('h5ad/CS10/CS10_2_D03251G111.h5ad')
#slice1, slice2

slice1.obs['x'] = slice1.obs['x'].astype(int)
slice1.obs['y'] = slice1.obs['y'].astype(int)
slice2.obs['x'] = slice2.obs['x'].astype(int)
slice2.obs['y'] = slice2.obs['y'].astype(int)
slice1.obsm['spatial'] = slice1.obs.loc[:,['x','y']]
slice2.obsm['spatial'] = slice2.obs.loc[:,['x','y']]
slice1.obsm['spatial'] = slice1.obsm['spatial'].values
slice2.obsm['spatial'] = slice2.obsm['spatial'].values
#slice1, slice2

spatial_key = 'spatial'
cluster_key = 'celltype'
st.pl.slices_2d(
    slices = [slice1, slice2],
    label_key = cluster_key,
    spatial_key = spatial_key,
    height=2,
    center_coordinate=True,
    show_legend=True,
    legend_kwargs={'loc': 'upper center', 'bbox_to_anchor': (0.5, 0) ,'ncol': 6, 'borderaxespad': -2, 'frameon': False},
)

# preprocess slice1
sc.pp.filter_cells(slice1, min_genes=1)  # we use min_genes=10 as 100 is too large for ST data
sc.pp.filter_genes(slice1, min_cells=1)
# Saving count data
slice1.layers["counts"] = slice1.X.copy()
# Normalizing to median total counts
sc.pp.normalize_total(slice1)
# Logarithmize the data
sc.pp.log1p(slice1)
# annotates highly variable genes
sc.pp.highly_variable_genes(slice1, n_top_genes=2000)

# preprocess slice2
sc.pp.filter_cells(slice2, min_genes=1)
sc.pp.filter_genes(slice2, min_cells=1)
# Saving count data
slice2.layers["counts"] = slice2.X.copy()
# Normalizing to median total counts
sc.pp.normalize_total(slice2)
# Logarithmize the data
sc.pp.log1p(slice2)
# annotates highly variable genes
sc.pp.highly_variable_genes(slice2, n_top_genes=2000)
#slice1, slice2

key_added = 'align_spatial'
# spateo return aligned slices as well as the mapping matrix
aligned_slices, pis = st.align.morpho_align(
    models=[slice1, slice2],
    ## Uncomment this if use highly variable genes
    # models=[slice1[:, slice1.var.highly_variable], slice2[:, slice2.var.highly_variable]],
    ## Uncomment the following if use pca embeddings
    # rep_layer='X_pca',
    # rep_field='obsm',
    # dissimilarity='cos',
    verbose=False,
    spatial_key=spatial_key,
    key_added=key_added,
    device=device,
)

st.pl.overlay_slices_2d(slices = aligned_slices, spatial_key = key_added+'_nonrigid', height=5, overlay_type='backward')

df = pd.DataFrame(aligned_slices[0].obsm['align_spatial_nonrigid'])  
df.to_csv('11.3D/SPATEO/CS10/CS10_1_D03251F413_align.csv', index=False,header=False) 
df = pd.DataFrame(aligned_slices[1].obsm['align_spatial_nonrigid'])  
df.to_csv('11.3D/SPATEO/CS10/CS10_2_D03251G111_align.csv', index=False,header=False)  

allslices = ["CS10_2_D03251G111","CS10_3_D03252G213","CS10_4_D03254C314","CS10_5_D03255E612",
             "CS10_6_D03257C314","CS10_7_D03257C412","CS10_8_D03257D411","CS10_9_D03257E112","CS10_10_D03257G212"]
             

for i in range(0,9) : 
    slice1 = st.read('h5ad/CS10/'+allslices[i]+'.h5ad')
    slice2 = st.read('h5ad/CS10/'+allslices[i+1]+'.h5ad')
    
    df = pd.read_csv('11.3D/SPATEO/CS9/'+allslices[i]+'_align.csv',header=None)
    slice1.obs['x'] = df.iloc[:,0].values
    slice1.obs['y'] = df.iloc[:,1].values
    
    slice2.obs['x'] = slice2.obs['x'].astype(int)
    slice2.obs['y'] = slice2.obs['y'].astype(int)
    slice1.obsm['spatial'] = slice1.obs.loc[:,['x','y']]
    slice2.obsm['spatial'] = slice2.obs.loc[:,['x','y']]
    slice1.obsm['spatial'] = slice1.obsm['spatial'].values
    slice2.obsm['spatial'] = slice2.obsm['spatial'].values
    spatial_key = 'spatial'
    cluster_key = 'celltype'
    st.pl.slices_2d(
        slices = [slice1, slice2],
        label_key = cluster_key,
        spatial_key = spatial_key,
        height=2,
        center_coordinate=True,
        show_legend=True,
        legend_kwargs={'loc': 'upper center', 'bbox_to_anchor': (0.5, 0) ,'ncol': 6, 'borderaxespad': -2, 'frameon': False},
    )
    # preprocess slice1
    sc.pp.filter_cells(slice1, min_genes=1)  # we use min_genes=10 as 100 is too large for ST data
    sc.pp.filter_genes(slice1, min_cells=1)
    # Saving count data
    slice1.layers["counts"] = slice1.X.copy()
    # Normalizing to median total counts
    sc.pp.normalize_total(slice1)
    # Logarithmize the data
    sc.pp.log1p(slice1)
    # annotates highly variable genes
    sc.pp.highly_variable_genes(slice1, n_top_genes=2000)
    # preprocess slice2
    sc.pp.filter_cells(slice2, min_genes=1)
    sc.pp.filter_genes(slice2, min_cells=1)
    # Saving count data
    slice2.layers["counts"] = slice2.X.copy()
    # Normalizing to median total counts
    sc.pp.normalize_total(slice2)
    # Logarithmize the data
    sc.pp.log1p(slice2)
    # annotates highly variable genes
    sc.pp.highly_variable_genes(slice2, n_top_genes=2000)
    #slice1, slice2
    key_added = 'align_spatial'
    # spateo return aligned slices as well as the mapping matrix
    aligned_slices, pis = st.align.morpho_align(
        models=[slice1, slice2],
        ## Uncomment this if use highly variable genes
        # models=[slice1[:, slice1.var.highly_variable], slice2[:, slice2.var.highly_variable]],
        ## Uncomment the following if use pca embeddings
        # rep_layer='X_pca',
        # rep_field='obsm',
        # dissimilarity='cos',
        verbose=False,
        spatial_key=spatial_key,
        key_added=key_added,
        device=device,
    )
    st.pl.overlay_slices_2d(slices = aligned_slices, spatial_key = key_added+'_nonrigid', height=5, overlay_type='backward')
    df = pd.DataFrame(aligned_slices[1].obsm['align_spatial_nonrigid'])  
    df.to_csv('11.3D/SPATEO/CS10/'+allslices[i+1]+'_align.csv', index=False,header=False)  