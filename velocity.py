import os
import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad

os.getcwd()

adata = scv.read("CD8scv.h5ad")

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)

#Plot velocity onto UMAP with celltypes
scv.pl.velocity_embedding_stream(adata, basis="umap", color="celltype", save = "velo.celltype")
scv.pl.velocity_embedding(adata, arrow_length=3, arrow_size=2, dpi=120)
#plot velocity of selected gene
scv.pl.velocity(adata, ['GZMK',  'GZMB', 'GNLY', 'IL7R'], ncols=2, save = "selected.genes.png")

#Identify important genes
scv.tl.rank_velocity_genes(adata, groupby='celltype', min_corr=.3)

df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head()

scv.pl.scatter(adata, df['0'][:3], ylabel='0', frameon=False, color='celltype', size=10, linewidth=1.5)

#Plot velo length and confidence
scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95], save = "velo.length.conf.png")

#Plot velo graph
scv.pl.velocity_graph(adata, threshold=.1, color='celltype', save = "velo.graph.png")

x, y = scv.utils.get_cell_transitions(adata, basis='umap', starting_cell=70)
ax = scv.pl.velocity_graph(adata, c='lightgrey', edge_width=.05, show=False)
ax = scv.pl.scatter(adata, x=x, y=y, s=120, c='ascending', cmap='gnuplot', ax=ax)

#Pseudotime
scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot', save = "velo.pseudotime")

#Run Paga
adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']

scv.tl.paga(adata, groups='celltype')
pd.set_option("display.precision", 2)
df = scv.get_df(adata, 'paga/transitions_confidence').T #Fix bug

scv.pl.paga(adata, basis='umap', size=50, alpha=.1, min_edge_width=2, node_size_scale=1.5, save = "paga.png")

#Dynamical Modeling
scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)

df = adata.var
df = df[(df['fit_likelihood'] > .1) & df['velocity_genes'] == True]

kwargs = dict(xscale='log', fontsize=16)

scv.tl.latent_time(adata)
top_genes = adata.var["fit_likelihood"].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata, var_names=top_genes, sortby="latent_time", col_color="celltype", n_convolve=100, save = "heat.velo")

