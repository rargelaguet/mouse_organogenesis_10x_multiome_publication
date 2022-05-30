######################
## Import libraries ##
######################

import os
from re import search
import SEACells

###########################
## Load default settings ##
###########################

if search("BI2404M", os.uname()[1]):
    exec(open('/Users/argelagr/gastrulation_multiome_10x/settings.py').read())
    exec(open('/Users/argelagr/gastrulation_multiome_10x/utils.py').read())
elif search("pebble|headstone", os.uname()[1]):
    exec(open('/bi/group/reik/ricard/scripts/gastrulation_multiome_10x/settings.py').read())
    exec(open('/bi/group/reik/ricard/scripts/gastrulation_multiome_10x/utils.py').read())
else:
    exit("Computer not recognised")

################################
## Initialise argument parser ##
################################

p = argparse.ArgumentParser( description='' )
p.add_argument( '--anndata',               type=str,                required=True,           help='Anndata file')
p.add_argument( '--metadata',               type=str,                required=True,           help='Cell metadata file')
p.add_argument( '--outdir',               type=str,                required=True,           help='Output directory')
p.add_argument( '--samples',            type=str, nargs="+",             default="all",             help='samples to use')
p.add_argument( '--percent_metacells',            type=float,              default=0.05,             help='Number of metacells (as a fraction of the total number of cells)')
p.add_argument( '--n_features',            type=int,              default=1500,             help='Number of features')
p.add_argument( '--n_pcs',            type=int,              default=25,             help='Number of PCs')
# p.add_argument( '--seed',                  type=int,                default=42,               help='Random seed')
# p.add_argument( '--n_iter',       type=int,              default=50,              help='Number of iterations')
args = p.parse_args()

## START TEST ##
# args = {}
# args["anndata"] = io["basedir"] + "/processed/rna/anndata.h5ad"
# args["metadata"] = io["basedir"] + "/results/rna/mapping/sample_metadata_after_mapping.txt.gz"
# args["outdir"] = io["basedir"] + "/results/rna/metacells/test"
# args["samples"] = ["E8.5_rep2"]
# args["percent_metacells"] = 0.05
# args["n_features"] = 1500
# args["n_pcs"] = 25
## END TEST ##

# convert args to dictionary
args = vars(args)

#####################
## Parse arguments ##
#####################

# I/O
# io["pca_rna"] = io["basedir"] + "/results/rna/dimensionality_reduction/all_cells/E7.5_rep1-E7.5_rep2-E8.0_rep1-E8.0_rep2-E8.5_rep1-E8.5_rep2_pca_features2500_pcs30_batchcorrectionbysample.txt.gz"
# io["pca_atac"] = io["basedir"] + "/results/atac/archR/dimensionality_reduction/PeakMatrix/all_cells/E7.5_rep1-E7.5_rep2-E8.0_rep1-E8.0_rep2-E8.5_rep1-E8.5_rep2_umap_nfeatures50000_ndims50_neigh45_dist0.45.txt.gz"

if not os.path.isdir(args["outdir"]): os.makedirs(args["outdir"])
args["outdir"] = Path(args["outdir"])

sc.settings.figdir = args["outdir"] / "pdf"

if isinstance(args["samples"],list): 
  if args["samples"][0]=="all":
    args["samples"] = opts["samples"]
  else:
    assert set(args["samples"]).issubset(opts["samples"])

else:
  print('args["samples"] has to be a list')

print(args)

###################
## Load metadata ##
###################

metadata = (pd.read_table(args["metadata"]) >>
    # mask(X["pass_rnaQC"]==True, X["pass_atacQC"]==True, X["doublet_call"]==False, X["celltype"].isin(opts["celltypes"])) >>
    mask(X["pass_rnaQC"]==True, X["doublet_call"]==False, X["celltype"].isin(opts["celltypes"])) >>
    mask(X["sample"].isin(args["samples"]))
).set_index("cell", drop=False)

# TO-DO: USE ONLY WT CELLS 

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

# Set colors
adata.obs = adata.obs.rename(columns={"celltype.mapped":"celltype"})
colPalette_celltypes = [opts["celltype_colors"][i.replace(" ","_")] for i in sorted(np.unique(adata.obs['celltype']))]
adata.uns['celltype_colors'] = colPalette_celltypes
#colPalette_stages = [opts["stages_colors"][i.replace(" ","_")] for i in sorted(np.unique(adata.obs['stage']))]
#adata.uns['stage_colors'] = colPalette_stages


#######################
## Feature selection ##
#######################

sc.pp.highly_variable_genes(adata, n_top_genes=args["n_features"])

##############################
## Dimensionality reduction ##
##############################

# Load precomputed PCA coordinates
# pca_mtx = pd.read_csv(io["pca_rna"]).set_index("cell", drop=True).loc[adata.obs.index].to_numpy()
# pca_mtx = pd.read_csv(io["pca_atac"]).set_index("cell", drop=True).loc[adata.obs.index].to_numpy()
# adata.obsm["X_pca"] = pca_mtx

# Run PCA
sc.tl.pca(adata, svd_solver='arpack', n_comps=args["n_pcs"])

# Plot PCA
# sc.pl.pca(adata, components=[1,2], color=["celltype","stage"], size=25, legend_loc=None)

# Build kNN graph
sc.pp.neighbors(adata, n_neighbors=25, use_rep='X_pca')

# Run UMAP
sc.tl.umap(adata, min_dist=0.5, n_components=2)

# Plot UMAP
sc.pl.umap(adata, color=["celltype"], size=25, legend_loc=None, save="_umap_cells.pdf")

########################
## Fit metacell model ##
########################

n_metacells = round(args["percent_metacells"] * adata.shape[0])

print("Fitting SEACells with %d metacells..." % (n_metacells))

model = SEACells.core.SEACells(adata, 
                  build_kernel_on = 'X_pca', 
                  n_SEACells = n_metacells, 
                  n_waypoint_eigs=10,
                  waypt_proportion=1,
                  convergence_epsilon = 1e-6)

model.fit()

adata.obs[['SEACell']].head()

#######################
## Plot model output ##
#######################

model.plot_convergence(save_as=args["outdir"] / "pdf/model_convergence.pdf")

SEACells.plot.plot_2D(adata, key='X_umap', colour_metacells=False, save_as=args["outdir"] / "pdf/umap_highlight_metacells.pdf")

################################################################
## Aggregate counts and plot trajectory at the metacell level ##
################################################################

adata_metacells = SEACells.core.summarize_by_SEACell(adata, SEACells_label='SEACell', summarize_layer='raw')
adata_metacells.uns = adata.uns
adata_metacells.obs = (adata.obs.loc[adata_metacells.obs.index] >> 
    select(["sample","celltype","genotype"])
)
sc.pp.normalize_total(adata_metacells)
sc.pp.log1p(adata_metacells)
sc.pp.highly_variable_genes(adata_metacells, n_top_genes=2500)
sc.tl.pca(adata_metacells, n_comps=25)
sc.pp.neighbors(adata_metacells, n_neighbors=25, use_rep='X_pca')
sc.tl.umap(adata, min_dist=0.5, n_components=2)
sc.pl.umap(adata, color=["celltype"], size=25, legend_loc=None, save="_umap_metacells.pdf")

##########
## Save ##
##########

to_save = adata.obs[['SEACell']].reset_index()
to_save.columns = ["cell","metacell"]

outfile = args["outdir"] / "cell2metacell_assignment.txt.gz"
to_save.to_csv(outfile, sep="\t", header=True, index=False)

adata_metacells.write_h5ad(args["outdir"] / "anndata_metacells.h5ad")
