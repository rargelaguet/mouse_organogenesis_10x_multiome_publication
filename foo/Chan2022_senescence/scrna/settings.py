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
if search("^(BAI-|SDI-|CAM-)", host):
	io["basedir"] = "/Users/rargelaguet/data/Chan2022_senescence/scrna"
	io["gene_metadata"] = "/Users/rargelaguet/data/ensembl/human/v99/bioMart/all_genes/Hsapiens_genes_BioMart.99.txt.gz"
else:
	io["basedir"] = "/home/rargelaguet/data/Chan2022_senescence/scrna"
	io["gene_metadata"] = "/home/rargelaguet/data/ensembl/human/v99/bioMart/all_genes/Hsapiens_genes_BioMart.99.txt.gz"


io["metadata"] = io["basedir"] + "/processed/cellranger/cell_metadata.txt.gz"
io["anndata"] = io["basedir"] + "/processed/cellranger/anndata.h5ad"

#############
## Options ##
#############

opts = {}

opts["samples"] = [
  "sample-SRR14646407_GSM5341090_RS_sc_hTert_t2",
  "sample-SRR14646408_GSM5341091_RS_sc_hTert_t4",
  "sample-SRR14646409_GSM5341092_RS_sc_hTert_t5",
  "sample-SRR14646410_GSM5341093_RS_sc_hTert_t6",
  "sample-SRR14646411_GSM5341094_RS_sc_hTert_t7",
  "sample-SRR14646412_GSM5341095_RS_sc_PDL_25_1",
  "sample-SRR14646413_GSM5341096_RS_sc_PDL_29_1",
  "sample-SRR14646414_GSM5341097_RS_sc_PDL_29_2",
  "sample-SRR14646415_GSM5341098_RS_sc_PDL_33_1",
  "sample-SRR14646416_GSM5341099_RS_sc_PDL_33_2",
  "sample-SRR14646417_GSM5341100_RS_sc_PDL_37_1",
  "sample-SRR14646418_GSM5341101_RS_sc_PDL_37_2",
  "sample-SRR14646419_GSM5341102_RS_sc_PDL_46_1",
  "sample-SRR14646420_GSM5341103_RS_sc_PDL_46_2",
  "sample-SRR14646421_GSM5341104_RS_sc_PDL_50_1",
  "sample-SRR14646422_GSM5341105_RS_sc_PDL_50_2"
]

opts["aliases"] = [
]

