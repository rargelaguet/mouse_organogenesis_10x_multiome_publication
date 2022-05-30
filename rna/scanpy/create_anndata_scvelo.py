import argparse
import scvelo as scv

######################
## Define arguments ##
######################

p = argparse.ArgumentParser( description='' )
p.add_argument( '--loom_directory',               type=str,                       help='Input directory for the loom files (after velocyto)' )
p.add_argument( '--anndata',               type=str,                       help='Anndata object' )
p.add_argument( '--outfile',               type=str,                       help='Output file (anndata)' )
p.add_argument( '--metadata',               type=str,                       help='Metadata file' )
p.add_argument( '--samples',               type=str, nargs="+",                       help='Samples' )
args = p.parse_args()

#####################
## Define settings ##
#####################

exec(open('../../settings.py').read())
exec(open('../../utils.py').read())

## START TEST ##
args.outfile = io["basedir"]+"/processed/rna/anndata_scvelo.h5ad"
args.anndata = io["basedir"]+"/processed/rna/anndata.h5ad"
args.loom_directory = io["basedir"]+"/processed/rna/loom"
args.metadata = io["basedir"]+"/results_new2/rna/mapping/sample_metadata_after_mapping.txt.gz"
args.samples = ["E7.5_rep1", "E7.5_rep2"]
## END TEST ##


###################
## Load metadata ##
###################

print("Loading metadata...")

metadata = (pd.read_table(args.metadata) >>
    mask(X.pass_rnaQC==True, X.doublet_call==False) >>
    mask(X["sample"].isin(args.samples))
).set_index("cell", drop=False)
print(metadata.shape)

#########################
## Load anndata object ##
#########################

print("Loading anndata...")

adata = load_adata(adata_file = args.anndata, metadata_file = args.metadata, normalise = False, cells = metadata.index.values)

###############################################################
## Load spliced and unspliced count matrices from loom files ##
###############################################################

print("Loading loom files...")

looms = [None for i in range(len(args.samples))]

for i in range(len(args.samples)):
    loom_file = args.loom_directory + "/" + args.samples[i] + ".loom"
    looms[i] = sc.read_loom(loom_file, sparse=True, X_name='spliced', obs_names='CellID', obsm_names=None, var_names='Gene')
    # looms[i].var_names_make_unique()
    # looms[i].obs.index = looms[i].obs.index.str.replace(rename_dict[args.samples[i]]+":",args.samples[i]+"_").str.replace("x","-1")
    # print(looms[i].shape)
    # print(looms[i].obs.head())

####################
## Create anndata ##
####################

print("Creating anndata file...")

# Concatenate
adata_loom = anndata.AnnData.concatenate(*looms, join='inner', batch_key=None, index_unique=None)
del looms

# Remove non-used layers to save memory
del adata_loom.layers["ambiguous"]
del adata_loom.layers["matrix"]

# Merge anndata objects
adata_final = scv.utils.merge(adata, adata_loom)
del adata_loom
del adata
adata_final

adata_final.obs.index.name = None

##########
## Save ##
##########

print("Saving anndata object...")

adata.write_h5ad(args.outfile)