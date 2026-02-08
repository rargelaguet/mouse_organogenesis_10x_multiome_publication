######################
## Import libraries ##
######################

import os
import SEACells
from re import search

###########################
## Load default settings ##
###########################

host = os.uname()[1]
if search("BAI-JTYCW2FJ91", host):
    exec(open('/Users/rargelaguet/altos_github/fibrosis/Peyser2019_fibrosis/settings.py').read())
    exec(open('/Users/rargelaguet/altos_github/fibrosis/Peyser2019_fibrosis/utils.py').read())
else:
    exec(open('/share/rargelaguet/altos_github/fibrosis/Peyser2019_fibrosis/settings.py').read())
    exec(open('/share/rargelaguet/altos_github/fibrosis/Peyser2019_fibrosis/utils.py').read())

# sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(8, 7), facecolor='white')
sc.settings.autoshow = False

################################
## Initialise argument parser ##
################################

p = argparse.ArgumentParser( description='' )
p.add_argument( '--anndata',               type=str,                required=True,           help='Anndata file')
p.add_argument( '--metadata',               type=str,                required=True,           help='Cell metadata file')
p.add_argument( '--outdir',               type=str,                required=True,           help='Output directory')
p.add_argument( '--samples',            type=str, nargs="+",             default="all",             help='samples to use')
p.add_argument( '--percent_metacells',            type=float,              default=0.05,             help='Number of metacells (as a fraction of the total number of cells)')
args = p.parse_args()

## START TEST ##
# args = {}
# args["anndata"] = io["basedir"] + "/processed/anndata.h5ad"
# args["metadata"] = io["basedir"] + "/processed/cell_metadata.txt.gz"
# args["outdir"] = io["basedir"] + "/results/metacells"
# args["percent_metacells"] = 0.05
# args["samples"] = ["Bone"]
## END TEST ##

# convert args to dictionary
args = vars(args)
print(args)

#####################
## Parse arguments ##
#####################

if not os.path.isdir(args["outdir"]): os.makedirs(args["outdir"])
args["outdir"] = Path(args["outdir"])
os.makedirs(args["outdir"]/"pdf", exist_ok=True)

sc.settings.figdir = args["outdir"]/"pdf"

# if isinstance(args["samples"],list): 
#   if args["samples"][0]=="all":
#     args["samples"] = opts["samples"]
#   else:
#     assert set(args["samples"]).issubset(opts["samples"])

# else:
#   print('args["samples"] has to be a list')

# print(args)

###################
## Load metadata ##
###################

# metadata = pd.read_table(args["metadata"]).set_index("cell", drop=False)

metadata = (pd.read_table(args["metadata"]) >>
    mask(X["condition"].isin(args["samples"]))
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


# Filter genes
# adata = adata[:,adata.layers["raw"].sum(axis=0)>=50]
sc.pp.filter_genes(adata, min_counts=50)

#######################
## Feature selection ##
#######################

sc.pp.highly_variable_genes(adata, n_top_genes=2500)

##############################
## Dimensionality reduction ##
##############################

# Run PCA
# sc.tl.pca(adata, svd_solver='arpack', n_comps=25)

# Plot PCA
# sc.pl.pca(adata, components=[1,2], color=["celltype","stage"], size=25, legend_loc=None)

# Build kNN graph
sc.pp.neighbors(adata, n_neighbors=25, use_rep='X_pca')
# sc.pp.neighbors(adata, n_neighbors=25, use_rep='pca')

# Run UMAP
# sc.tl.umap(adata, min_dist=0.5, n_components=2)

# Plot UMAP
# sc.pl.umap(adata, color=["alias"], size=25, legend_loc=None, save="_umap_cells.pdf")

########################
## Fit metacell model ##
########################

n_metacells = round(args["percent_metacells"] * adata.shape[0])

print("Fitting SEACells with %d metacells..." % (n_metacells))

model = SEACells.core.SEACells(adata, 
                  build_kernel_on = 'X_pca',
                  # build_kernel_on = 'pca', 
                  n_SEACells = n_metacells, 
                  n_waypoint_eigs=10,
                  convergence_epsilon = 1e-5)

# Construct kernel matrix
model.construct_kernel_matrix()

# Initialize archetypes
model.initialize_archetypes()

# Plot the initilization to ensure they are spread across phenotypic space
# SEACells.plot.plot_initialization(adata, model)

# Fit the model
model.fit(min_iter=10, max_iter=100)
# model.fit(min_iter=1, max_iter=1)

# adata.obs[['SEACell']].head()

#######################
## Plot model output ##
#######################

model.plot_convergence(save_as=args["outdir"] / "pdf/model_convergence.pdf", show=False)

SEACells.plot.plot_2D(adata, key='X_umap', colour_metacells=False, save_as=args["outdir"] / "pdf/umap_highlight_metacells.pdf", show=False)

#######################
## Aggregate counts  ##
#######################

adata_metacells = SEACells.core.summarize_by_SEACell(adata, SEACells_label='SEACell', summarize_layer='raw')

#########################################################
## Plot dimensionality reduction at the metacell level ##
#########################################################

# adata_metacells.uns = adata.uns
# sc.pp.normalize_total(adata_metacells)
# sc.pp.log1p(adata_metacells)
# sc.pp.highly_variable_genes(adata_metacells, n_top_genes=2500)
# sc.tl.pca(adata_metacells, n_comps=25)
# sc.pp.neighbors(adata_metacells, n_neighbors=25, use_rep='X_pca')
# sc.tl.umap(adata, min_dist=0.5, n_components=2)
# sc.pl.umap(adata, color=["alias"], size=25, legend_loc=None, save="_umap_metacells.pdf")

##########
## Save ##
##########

cell2metacell = model.get_hard_assignments()
cell2metacell.reset_index(inplace=True)
cell2metacell.columns = ["cell","metacell"]
cell2metacell.to_csv(args["outdir"]/"cell2metacell_assignment.txt.gz", sep="\t", header=True, index=False)
# adata_metacells.write_h5ad(args["outdir"] / "anndata_metacells.h5ad")

##########
## TEST ##
##########

# umap = pd.DataFrame(ad.obsm[key]).set_index(ad.obs_names).join(ad.obs['SEACell'])
# umap['SEACell'] = umap['SEACell'].astype("category")
# mcs = umap.groupby('SEACell').mean().reset_index()
