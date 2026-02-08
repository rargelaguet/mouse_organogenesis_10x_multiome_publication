
import os
import dfply
from re import search

#####################
## Define settings ##
#####################

host = os.uname()[1]
if search("^(BAI-|SDI-)", host):
    exec(open('/Users/rargelaguet/altos_github/fibrosis/Peyser2019_fibrosis/settings.py').read())
    exec(open('/Users/rargelaguet/altos_github/fibrosis/Peyser2019_fibrosis/utils.py').read())
else:
    exec(open('/share/rargelaguet/altos_github/fibrosis/Peyser2019_fibrosis/settings.py').read())
    exec(open('/share/rargelaguet/altos_github/fibrosis/Peyser2019_fibrosis/utils.py').read())

io["metadata"] = io["basedir"] + "/processed/cell_metadata.txt.gz"
io["anndata"] = io["basedir"] + "/processed/anndata.h5ad"
io["outdir"] = io["basedir"] + "/results/dimensionality_reduction/cells"


opts["samples"] = [
  "", 
  "", 
  ""
]

########################
## Load cell metadata ##
########################

cell_metadata = pd.read_table(io["metadata"])
# cell_metadata = (pd.read_table(io["metadata"]) >>
#     mask(X["sample"].isin(opts["samples"]))
# ).set_index("cell", drop=False)

print(cell_metadata.shape)
print(cell_metadata["sample"].value_counts())

#########################
## Load anndata object ##
#########################

# adata = sc.read(io["anndata"])
adata = load_adata(io["anndata"], metadata_file = io["metadata"], filter_lowly_expressed_genes = True, normalise = True, cell_column = "cell", cells = cell_metadata["cell"].values)
adata

#######################
## Feature selection ##
#######################

sc.pp.highly_variable_genes(adata, n_top_genes=2500)

#########
## PCA ##
#########

sc.tl.pca(adata, svd_solver='arpack')

# sc.pl.pca(adata, components=['1,2'], color=["sample"], size=25)

###############
## kNN graph ##
###############

sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30)

##########
## UMAP ##
##########

sc.tl.umap(adata, min_dist=0.5, n_components=2)
sc.pl.umap(adata, color=["sample"], size=25, legend_loc="on data",  legend_fontsize=8)
sc.pl.umap(adata, color=["cell_type"], size=25, legend_loc="on data",  legend_fontsize=8)
sc.pl.umap(adata, color=["condition"], size=25, legend_fontsize=8)

###########################
## Force-directed layout ##
###########################

sc.tl.draw_graph(adata, layout='fa')
sc.pl.draw_graph(adata, color=["sample"], legend_loc='on data', size=20)

##########################
## Plot gene expression ##
##########################

"HERVK-int" in adata.var_names
sc.pl.draw_graph(adata, color="HERVK-int", legend_loc='on data')

sc.pl.umap(adata, color="ENSG00000118137")

sc.pl.violin(adata, ['ELN','ACTA2'], groupby='cell_type', jitter=False, inner='box')

###############################################
## Save dimensionality reduction coordinates ##
###############################################

# save scVI coordinates
# tmp = pd.DataFrame(np.around(adata.obsm["X_scVI"],3), index=adata.obs_names, columns=["scvi_%s" % i for i in range(adata.obsm["X_scVI"].shape[1])])
# tmp.to_csv(io["outdir"]+"/X_scVI.txt")

# save FA coordinates
tmp = pd.DataFrame(np.around(adata.obsm["X_draw_graph_fa"],3), index=adata.obs_names, columns=["fa_%s" % i for i in range(adata.obsm["X_draw_graph_fa"].shape[1])])
tmp.to_csv(io["outdir"]+"/X_draw_graph_fa.txt")

########################
## Save cell metadata ##
########################

cell_metadata = adata.obs.copy()

cell_metadata.drop(columns=['fastq_1','fastq_2','log1p_n_genes_by_counts','log1p_total_counts','log1p_total_counts_mt','log1p_total_counts_ribo','log1p_total_counts_hb'], inplace=True)

cell_metadata["cell"] = adata.obs.index.values

cell_metadata.to_csv(io["metadata"], sep="\t", header=True, index=False)
