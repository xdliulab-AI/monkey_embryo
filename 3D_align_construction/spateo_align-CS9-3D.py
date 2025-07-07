import warnings

import numpy as np
from scipy.spatial import KDTree
from sklearn.decomposition import PCA
import spateo as st
warnings.filterwarnings('ignore')

import os
import pyvista as pv

adata = st.read_h5ad('CS9.h5ad')  
import pandas as pd
df = pd.read_csv('11.3D/SPATEO/CS9_SPATEO_align_metadata.csv',index_col=0)
print(df)

adata.obsm['3d_align_spatial'] = df.loc[:,['x','y','z']]

cpo = [(553, 1098, 277), (1.967, -6.90, -2.21), (0, 0, 1)]

adata.uns["__type"] = "UMI"
adata
adata.obsm['3d_align_spatial']
embryo_pc, plot_cmap = st.tdr.construct_pc(adata=adata.copy(), spatial_key="3d_align_spatial", groupby="celltype", key_added="celltype", colormap="rainbow")
st.tdr.save_model(model=embryo_pc, filename="D11.3D/SPATEO/CS9_SPATEO_embryo_pc_model.vtk")

embryo_mesh, _, _ = st.tdr.construct_surface(pc=embryo_pc, key_added="celltype", alpha=0.6, cs_method="marching_cube", cs_args={"mc_scale_factor": 0.8}, smooth=5000, scale_factor=1.08)
st.pl.three_d_plot(model=st.tdr.collect_models([embryo_mesh, embryo_pc]), key="celltype", model_style=["surface", "points"], jupyter="static", cpo=cpo)
st.tdr.save_model(model=embryo_mesh, filename="11.3D/SPATEO/CS9_SPATEO_embryo_mesh_model.vtk")

embryo_voxel, _ = st.tdr.voxelize_mesh(mesh=embryo_mesh, voxel_pc=None, key_added="celltype", label="embryo_voxel", color="gainsboro", smooth=500)
#st.pl.three_d_plot(model=embryo_voxel, key="celltype", jupyter="static", cpo=cpo)
st.tdr.save_model(model=embryo_voxel, filename="11.3D/SPATEO/CS9_SPATEO_embryo_voxel_model.vtk")