######################
## Import libraries ##
######################

import os
from re import search

###########################
## Load default settings ##
###########################

exec(open('/share/rargelaguet/altos_github/reprogramming/Xing2020_bj_oskm/settings.py').read())
exec(open('/share/rargelaguet/altos_github/reprogramming/Xing2020_bj_oskm/utils.py').read())

sc.settings.autoshow = False

################################
## Initialise argument parser ##
################################

p = argparse.ArgumentParser( description='' )
p.add_argument( '--anndata',               type=str,                required=True,           help='Anndata file')
p.add_argument( '--metadata',               type=str,                required=True,           help='Cell metadata file')
p.add_argument( '--outdir',               type=str,                required=True,           help='Output directory')
p.add_argument( '--samples',            type=str, nargs="+",             default="all",             help='samples to use')
args = p.parse_args()

## START TEST ##
# args = {}
# args["anndata"] = io["basedir"] + "/processed/rna/anndata.h5ad"
# args["metadata"] = io["basedir"] + "/processed/rna/cell_metadata.txt.gz"
# args["outdir"] = io["basedir"] + "/results/rna/dimensionality_reduction/scRNA#GSM7763424_D10"
# args["samples"] = ["scRNA#GSM7763424_D10"]
## END TEST ##

# convert args to dictionary
args = vars(args)

print(args)

#####################
## Parse arguments ##
#####################

# I/O
args["outdir"] = Path(args["outdir"])
if not os.path.isdir(args["outdir"]): 
  os.makedirs(args["outdir"]/"pdf")

sc.settings.figdir = args["outdir"]/"pdf"

# if isinstance(args["samples"],list): 
#   if args["samples"][0]=="all":
#     args["samples"] = opts["samples"]
#   else:
#     assert set(args["samples"]).issubset(opts["samples"])

# else:
#   print('args["samples"] has to be a list')

###################
## Load metadata ##
###################

metadata = (pd.read_table(args["metadata"]) >>
    # mask(X["sample"].isin(args["samples"]))
    mask(X["group"].isin(args["samples"]))
).set_index("cell", drop=False)

print(metadata.shape)
print(metadata.head())

##################
## Load AnnData ##
##################

adata = load_adata(
    adata_file = args["anndata"], 
    metadata_file = args["metadata"], 
    cells = metadata.index.values, 
    normalise = True, 
    keep_counts = True,
    filter_lowly_expressed_genes = True, 
    set_colors = False
)

#######################
## Feature selection ##
#######################

sc.pp.highly_variable_genes(adata, n_top_genes=1500)

#########
## PCA ##
#########

sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca(adata, color=["sample"], size=25, legend_loc="on data",  legend_fontsize=8, save=True)

###############
## kNN graph ##
###############

sc.pp.neighbors(adata, n_neighbors=25, n_pcs=25)

##########
## UMAP ##
##########

sc.tl.umap(adata, min_dist=0.5, n_components=2)
sc.pl.umap(adata, color=["sample"], size=25, legend_loc="on data",  legend_fontsize=8, save=True)

###########################
## Force-directed layout ##
###########################

sc.tl.draw_graph(adata, layout='fa')
sc.pl.draw_graph(adata, color=["sample"], legend_loc='on data', size=20, save=True)

######################
## Save coordinates ##
######################

# save PCA coordinates
tmp = pd.DataFrame(np.around(adata.obsm["X_pca"],3), index=adata.obs_names, columns=["fa_%s" % i for i in range(adata.obsm["X_pca"].shape[1])])
tmp.to_csv(args["outdir"]/"X_pca.txt")

# save UMAP coordinates
tmp = pd.DataFrame(np.around(adata.obsm["X_umap"],3), index=adata.obs_names, columns=["fa_%s" % i for i in range(adata.obsm["X_umap"].shape[1])])
tmp.to_csv(args["outdir"]/"X_umap.txt")

# save FA coordinates
tmp = pd.DataFrame(np.around(adata.obsm["X_draw_graph_fa"],3), index=adata.obs_names, columns=["fa_%s" % i for i in range(adata.obsm["X_draw_graph_fa"].shape[1])])
tmp.to_csv(args["outdir"]/"X_fa.txt")
