
import os
import dfply
from re import search

#####################
## Define settings ##
#####################

host = os.uname()[1]
if search("^(BAI-|SDI-)", host):
    exec(open('/Users/rargelaguet/altos_github/reprogramming/Xing2020_bj_oskm/settings.py').read())
    exec(open('/Users/rargelaguet/altos_github/reprogramming/Xing2020_bj_oskm/utils.py').read())
else:
    exec(open('/share/rargelaguet/altos_github/reprogramming/Xing2020_bj_oskm/settings.py').read())
    exec(open('/share/rargelaguet/altos_github/reprogramming/Xing2020_bj_oskm/utils.py').read())

io["metadata"] = io["basedir"] + "/processed/scTE/cell_metadata.txt.gz"
io["anndata"] = io["basedir"] + "/processed/scTE/anndata.h5ad"
io["outdir"] = io["basedir"] + "/results/dimensionality_reduction/cells"


opts["samples"] = [
  "GM00731_D0_S2", 
  "GM00731_D3_S4", 
  "GM00731_D7_S6", 
  "GM00731_D10_S8" 
  # "GM23815_D0_S1", 
  # "GM23815_D3_S3", 
  # "GM23815_D7_S5", 
  # "GM23815_D10_S7"
]

########################
## Load cell metadata ##
########################

cell_metadata = (pd.read_table(io["metadata"]) >>
    mask(X["sample"].isin(opts["samples"]))
).set_index("cell", drop=False)

print(cell_metadata.shape)
print(cell_metadata["sample"].value_counts())

#########################
## Load anndata object ##
#########################

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

sc.pl.violin(adata, ['LMNB1', 'LMNA'], groupby='sample', jitter=False, inner='box')

###############################################
## Save dimensionality reduction coordinates ##
###############################################

# save scVI coordinates
# tmp = pd.DataFrame(np.around(adata.obsm["X_scVI"],3), index=adata.obs_names, columns=["scvi_%s" % i for i in range(adata.obsm["X_scVI"].shape[1])])
# tmp.to_csv(io["outdir"]+"/X_scVI.txt")

# save FA coordinates
tmp = pd.DataFrame(np.around(adata.obsm["X_draw_graph_fa"],3), index=adata.obs_names, columns=["fa_%s" % i for i in range(adata.obsm["X_draw_graph_fa"].shape[1])])
tmp.to_csv(io["outdir"]+"/X_draw_graph_fa.txt")
