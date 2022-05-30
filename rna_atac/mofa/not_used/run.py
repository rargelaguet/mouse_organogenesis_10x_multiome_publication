from mofapy2.run.entry_point import entry_point
import pandas as pd
import numpy as np
import argparse
from scipy.io import mmread, mmwrite
from pathlib import Path

################################
## Initialise argument parser ##
################################

p = argparse.ArgumentParser( description='' )
# p.add_argument( '--input_folder',          type=str,              required=True,          help='Input data file (matrix format)' )
p.add_argument( '--rna_matrix',          type=str,              required=True,          help='' )
p.add_argument( '--atac_matrix',          type=str,              required=True,          help='' )
p.add_argument( '--rna_features',          type=str,              required=True,          help='' )
p.add_argument( '--atac_features',          type=str,              required=True,          help='' )
p.add_argument( '--cells',          type=str,              required=True,          help='' )
p.add_argument( '--outfile',               type=str,              required=True,          help='Output file to store the model (.hdf5)' )
p.add_argument( '--factors',               type=int,              default=25,             help='Number of factors' )
p.add_argument( '--seed',                  type=int,              default=42,             help='Random seed' )
p.add_argument( '--convergence_mode',      type=str,              default="fast",       help='Convergence mode')
p.add_argument( '--test',               action="store_true",                           help='Do stochastic inference?' )
args = p.parse_args()

## START TEST ##
# /bi/group/reik/ricard/data/gastrulation_multiome_10x
# /Users/argelagr/data/gastrulation_multiome_10x
# args = {}
# args["input_folder"] = "/bi/group/reik/ricard/data/gastrulation_multiome_10x/test/results/rna_atac/mofa/all_cells"
# args["rna_matrix"] = "/bi/group/reik/ricard/data/gastrulation_multiome_10x/test/results/rna_atac/mofa/all_cells/rna.mtx"
# args["atac_matrix"] = "/bi/group/reik/ricard/data/gastrulation_multiome_10x/test/results/rna_atac/mofa/all_cells/atac_tfidf.mtx"
# args["rna_features"] = "/bi/group/reik/ricard/data/gastrulation_multiome_10x/test/results/rna_atac/mofa/all_cells/rna_features.txt"
# args["atac_features"] = "/bi/group/reik/ricard/data/gastrulation_multiome_10x/test/results/rna_atac/mofa/all_cells/atac_features.txt"
# args["cells"] = "/bi/group/reik/ricard/data/gastrulation_multiome_10x/test/results/rna_atac/mofa/all_cells/cells.txt"
# args["outfile"] = "/bi/group/reik/ricard/data/gastrulation_multiome_10x/test/results/rna_atac/mofa/all_cells/test.hdf5"
# args["factors"] = 25
# args["seed"] = 42
# args["convergence_mode"] = "fast"
# args["test"] = True
## END TEST ##

# convert args to dictionary
args = vars(args)

###############
## Load data ##
###############

rna_mtx = mmread(args["rna_matrix"]).todense().T
atac_mtx = mmread(args["atac_matrix"]).todense().T

rna_features = pd.read_csv(args["rna_features"], header=None)[0].tolist()
atac_features = pd.read_csv(args["atac_features"], header=None)[0].tolist()
cells = pd.read_csv(args["cells"], header=None)[0].tolist()

# sample_metadata = pd.read_csv(input_folder/"sample_metadata.txt.gz")
# assert sample_metadata.cell.tolist() == cells

########################
## Create MOFA object ##
########################

# initialise entry point    
ent = entry_point()

# Set data
ent.set_data_matrix(
	data = [[rna_mtx], [atac_mtx]], 
	views_names = ["RNA","ATAC"], 
	samples_names = [ cells ], 
	features_names = [ rna_features, atac_features ]
)

# Set data options
ent.set_data_options(use_float32 = True)

# Set model options
ent.set_model_options(factors=args["factors"], spikeslab_factors=False, spikeslab_weights=False)

# Set training options
if args["test"]:
	ent.set_train_options(iter=3)
else:
	ent.set_train_options(convergence_mode=args["convergence_mode"], seed=args["seed"])

###############################
## Build and train the model ##
###############################

# Build the model
ent.build()

# Train the model
ent.run()

####################
## Save the model ##
####################

ent.save(args["outfile"], save_data=False)

##########
## TEST ##
##########

# Check multithreading
# import os
# os.environ["OMP_NUM_THREADS"] = "4"
# os.environ["OPENBLAS_NUM_THREADS"] = "4"
# os.environ["MKL_NUM_THREADS"] = "6"
# os.environ["VECLIB_MAXIMUM_THREADS"] = "4"
# os.environ["NUMEXPR_NUM_THREADS"] = "6"
