exec(open('../../../settings.py').read())
exec(open('../../../utils.py').read())

######################
## Define arguments ##
######################

p = argparse.ArgumentParser( description='' )
p.add_argument('--anndata',               type=str,                       help='Anndata object' )
p.add_argument('--metadata',               type=str,                       help='Metadata file' )
p.add_argument('--outdir',               type=str,                       help='Output directory' )
p.add_argument('--nfeatures',        type=int,    default=1000,                help='Number of features')
p.add_argument('--n_pcs',            type=int,    default=30,                  help='Number of PCs')
p.add_argument('--n_neighbors',     type=int,    default=30,     help='(UMAP) Number of neighbours')
p.add_argument('--min_dist',        type=float,     default=0.3,     help='(UMAP) Minimum distance')
p.add_argument('--colour_by',       type=str,  default="celltype.mapped",  nargs='+',  help='Metadata columns to colour the UMAP by')
p.add_argument('--seed',            type=int,    default=42,                  help='Random seed')
args = p.parse_args()

## START TEST ##
args.anndata = io["anndata"]
args.metadata = io["metadata"] # io["basedir"]+"/results/rna/mapping/sample_metadata_after_mapping.txt.gz"
args.samples = ["E8.5_rep1"]
args.celltype_label = "celltype.mapped"
args.nfeatures = 1500
args.n_pcs = 15
args.min_dist = 0.3
args.colour_by = "celltype.mapped"
args.outdir = io["basedir"]+"/results/rna/dimensionality_reduction/scanpy"
## END TEST ##

#####################
## Define settings ##
#####################

# Define I/O
outdir = Path(args.outdir)
if not outdir.exists: outdir.mkdir()

# Define options
sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(8,7), facecolor='white', format="png")
# sc.settings.figdir(str(outdir))
# sc.settings.autosave = True

########################
## Load cell metadata ##
########################

metadata = (pd.read_table(args.metadata) >>
    mask(X["pass_rnaQC"]==True, X["doublet_call"]==False) >>
    mask(X["sample"].isin(args.samples), X[args.celltype_label].isin(opts["celltypes"]))
)

assert args.colour_by in metadata.columns.values

###########################
## Load anndata object ##
###########################

adata = load_adata(
    adata_file = args.anndata, 
    cells = metadata.cell.values, 
    normalise = True, 
    filter_lowly_expressed_genes = True,
    set_colors = False
)
print(adata)


# Set colour palettes
# colPalette_celltypes = [opts["celltype_colors"][i.replace(" ","_")] for i in sorted(np.unique(adata.obs[args.celltype_label]))]
# adata.uns[args.celltype_label+"_colors"] = colPalette_celltypes
# colPalette_stages = [opts["stages_colors"][i.replace(" ","_")] for i in sorted(np.unique(adata.obs['stage']))]
# adata.uns['stage_colors'] = colPalette_stages

#######################
## Feature selection ##
#######################

sc.pp.highly_variable_genes(adata, n_top_genes=args.nfeatures)

#####################################
## Regress out technical variables ##
#####################################

# sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])

#########
## PCA ##
#########

sc.tl.pca(adata, svd_solver='arpack')
# sc.pl.pca(adata, components=[1,2], color=["celltype.mapped"], size=25, legend_loc=None)

#############################
## Batch effect correction ##
#############################


#####################
## Build kNN graph ##
#####################

sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, n_pcs=args.n_pcs)

##########
## UMAP ##
##########

# outfile ="umap_features%d_pcs%d_neigh%d_dist%s_%s.pdf" % (args.nfeatures, args.n_pcs, args.n_neighbors, args.min_dist)
sc.tl.umap(adata, min_dist=args.min_dist, n_components=2)
sc.pl.umap(adata, color=args.colour_by, size=25, legend_loc="on data", save=str(outdir/"umap.pdf"))

