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
p.add_argument( '--outfile',               type=str,                required=True,           help='Output file (anndata)')
args = p.parse_args()

## START TEST ##
# args = {}
# # args["samples"] = ["E7.5_rep1", "E7.5_rep2", "E7.75_rep1", "E8.0_rep1", "E8.0_rep2", "E8.5_rep1", "E8.5_rep2", "E8.75_rep1", "E8.75_rep2", "E8.5_CRISPR_T_KO", "E8.5_CRISPR_T_WT"]
# args["samples"] = ["E7.5_rep1"]
# # args["loom_files"] = [io["basedir"]+"/original/%s/velocyto"%i for i in args["samples"]]
# args["anndata"] = io["basedir"] + "/processed/rna/anndata.h5ad"
# args["metadata"] = io["basedir"] + "/results/rna/mapping/sample_metadata_after_mapping.txt.gz"
# args["outfile"] = io["basedir"] + "/processed/rna/velocyto/anndata_velocyto.h5ad"
## END TEST ##

# convert args to dictionary
args = vars(args)


#####################
## Parse arguments ##
#####################

if not os.path.isdir(os.path.dirname(args["outfile"])): 
    os.makedirs(os.path.dirname(args["outfile"]))

if args["samples"]=="all":
    args["samples"] = opts["samples"]

###################
## Load metadata ##
###################

print("Loading metadata...")

metadata = (pd.read_table(args["metadata"]) >>
    mask(X.pass_rnaQC==True, X.doublet_call==False, X["sample"].isin(args["samples"]))
).set_index("cell", drop=False)

print(metadata.head())

##################
## Load anndata ##
##################

print("Loading anndata...")

adata = load_adata(
    adata_file = args["anndata"], 
    metadata_file = args["metadata"], 
    normalise = False, 
    cells = metadata.index.values
)
adata

##########################
## Load velocyto output ##
##########################

print("Loading velocyto output files...")

looms = [None for i in range(len(args["samples"]))]
for i in range(len(args["samples"])):
    print(args["samples"][i])
    io["loom_velocyto"] = "%s/%s/velocyto/%s.loom" % (io["cellranger_output"],args["samples"][i],args["samples"][i])
    looms[i] = sc.read_loom(io["loom_velocyto"], sparse=True, X_name='spliced', obs_names='CellID', obsm_names=None, var_names='Gene')
    looms[i].var_names_make_unique()
    looms[i].obs.index = looms[i].obs.index.str.replace(":","#").str.replace("x","-1")
    print(looms[i].shape)

###########################
## Parse velocyto output ##
###########################

print("Parsing velocyto output files...")

adata_loom = anndata.AnnData.concatenate(*looms, join='inner', batch_key=None, index_unique=None)
adata_loom
del looms

# Remove non-used layers to save memory
del adata_loom.layers["ambiguous"]
del adata_loom.layers["matrix"]

# Convert to sparse matrices
# from scipy.sparse import csr_matrix
# adata.X = csr_matrix(adata.X)

###########################
## Merge anndata objects ##
###########################

print("Merging original anndata with the velocyto anndata...")

adata_final = scv.utils.merge(adata, adata_loom)
adata_final

adata_final.obs.index.name = None

# TO-DO: TRANSFE .UNS AND .OBSM

# del adata_loom
# del adata

##########
## Save ##
##########

print("Saving output...")

# Logical columns need to be stored as str to avoid hdf5 complaining
# adata_final.obs["pass_rnaQC"] = adata_final.obs["pass_rnaQC"].astype("category")
adata_final.obs["pass_rnaQC"] = adata_final.obs["pass_rnaQC"].astype(str)
adata_final.obs["doublet_call"] = adata_final.obs["doublet_call"].astype(str)
adata_final.obs["genotype"] = adata_final.obs["genotype"].astype(str)
# adata_final.obs["pass_atacQC"] = adata_final.obs["pass_atacQC"].astype(str)
# adata_final.obs = adata_final.obs.drop(["pass_rnaQC","pass_atacQC","barcode","doublet_call"], axis=1)

adata_final.write_h5ad(args["outfile"], compression="gzip")
