import io
import numpy as np
import anndata as anndata
import scipy as s
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn
import pandas as pd
import csv
import os
from dfply import *
from re import search
from pathlib import Path
import argparse

#########
## I/O ##
#########

io = {}
host = os.uname()[1]
if search("^(BAI-|SDI-|CI-)", host):
	io["basedir"] = "/Users/rargelaguet/data/Peyser2019_fibrosis"
	io["gene_metadata"] = "/Users/rargelaguet/data/ensembl/mouse/v102/BioMart/all_genes/Mmusculus_genes_BioMart.102.txt"
# elif search("rargelaguet|dockyard", host):
else:
	io["basedir"] = "/share/rargelaguet/data/Peyser2019_fibrosis"
	io["gene_metadata"] = "/share/rargelaguet/data/ensembl/mouse/v102/BioMart/all_genes/Mmusculus_genes_BioMart.102.txt"
# else:
# 	print("Computer not recognised"); exit()


io["metadata"] = io["basedir"] + "/processed/cell_metadata.txt.gz"
io["anndata"] = io["basedir"] + "/processed/anndata.h5ad"

#############
## Options ##
#############

opts = {}
opts["condition"] = ["Saline Control", "Bleomycin Fibrosis Model", "Bleomycin Fibrosis plus Nintedanib"]
