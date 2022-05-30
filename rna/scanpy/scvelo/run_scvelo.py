######################
## Import libraries ##
######################

import os
from re import search
import scvelo as scv

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
p.add_argument( '--samples',            type=str, nargs="+",      required=True,       default="all",             help='Samples')
p.add_argument( '--celltypes',            type=str, nargs="+",      required=True,       default="all",             help='Celltypes')
p.add_argument( '--ncores',            type=int,              default=1,             help='Number of cores')
p.add_argument( '--outdir',               type=str,                required=True,           help='Output directory')
args = p.parse_args()

## START TEST ##
args = {}
args["anndata"] = io["basedir"] + "/processed/rna/velocyto/anndata_velocyto.h5ad"
args["metadata"] = io["basedir"] + "/results/rna/mapping/sample_metadata_after_mapping.txt.gz"
args["samples"] = "all"
args["celltypes"] = "all"
args["ncores"] = 1
args["outdir"] = io["basedir"] + "/processed/rna/velocyto"
## END TEST ##

# convert args to dictionary
args = vars(args)

#####################
## Parse arguments ##
#####################

if not os.path.isdir(args["outdir"]): os.makedirs(args["outdir"])
# if not os.path.isdir(os.path.dirname(args["outfile"])): 
#     os.makedirs(os.path.dirname(args["outfile"]))

if args["samples"]=="all":
    args["samples"] = opts["samples"]

if args["celltypes"]=="all":
    args["celltypes"] = opts["celltypes"]

print(args)

###################
## Load metadata ##
###################

print("Loading metadata...")

metadata = (pd.read_table(args["metadata"]) >>
    mask(X.pass_rnaQC==True, X.doublet_call==False, X["sample"].isin(args["samples"]), X["celltype"].isin(args["celltypes"]))
).set_index("cell", drop=False)

print(metadata.head())

##################
## Load anndata ##
##################

print("Loading anndata...")

adata = load_adata(
    adata_file = args["anndata"], 
    metadata_file = args["metadata"], 
    normalise = True, 
    cells = metadata.index.values
)
adata

assert "spliced" in list(adata.layers.keys())
assert "unspliced" in list(adata.layers.keys())

############
## scVelo ##
############

# Plot exonic vs intronic read proportions
scv.pl.proportions(adata, save = args["outdir"]+"/proportion_reads.pdf")

# Gene filter
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)

# Calculate means variances among nearest neighbors in PCA space
scv.pp.moments(adata, n_pcs=50, n_neighbors=25)

# Recover the full splicing kinetics of specified genes. 
# For each gene the model infers transcription rates, splicing rates, degradation rates, as well as cell-specific latent time and transcriptional states
scv.tl.recover_dynamics(adata, n_jobs=args["ncores"])

# Estimate velocities per gene
scv.tl.velocity(adata, mode="dynamical")

##########
## Save ##
##########

print("Saving output...")

# Rename index
adata.obs.index.name = "cells"

# delete unnecesary layers
del adata.layers["spliced"]

# cast some layers to float32 to reduce disk
adata.layers["fit_t"] = adata.layers["fit_t"].astype(np.float32)
adata.layers["fit_tau"] = adata.layers["fit_tau"].astype(np.float32)
adata.layers["fit_tau_"] = adata.layers["fit_tau_"].astype(np.float32)
adata.layers["velocity"] = adata.layers["velocity"].astype(np.float32)
adata.layers["velocity_u"] = adata.layers["velocity_u"].astype(np.float32)

# Logical columns need to be stored as str to avoid hdf5 complaining
adata.obs["pass_rnaQC"] = adata.obs["pass_rnaQC"].astype(str)
adata.obs["doublet_call"] = adata.obs["doublet_call"].astype(str)
adata.obs["genotype"] = adata.obs["genotype"].astype(str)

adata.write_h5ad(args["outdir"]+"/anndata_scvelo.h5ad", compression="gzip")
