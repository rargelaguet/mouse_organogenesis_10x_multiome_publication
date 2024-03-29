# Decoding gene regulation in the mouse embryo using single-cell multi-omics

This repository contains the scripts to reproduce the results of the manuscript [Decoding gene regulation in the mouse embryo using single-cell multi-omics](https://www.biorxiv.org/content/10.1101/2022.06.15.496239v1). 


Abstract
--------
Following gastrulation, the three primary germ layers develop into the major organs in a process known as organogenesis. Single-cell RNA sequencing has enabled the profiling of the gene expression dynamics of these cell fate decisions, yet a comprehensive map of the interplay between transcription factors and cis-regulatory elements is lacking, as are the underlying gene regulatory networks. Here we generate a multi-omics atlas of mouse early organogenesis by simultaneously profiling gene expression and chromatin accessibility from tens of thousands of single cells. We develop a computational method to leverage the multi-modal readouts to predict transcription factor binding events in cis-regulatory elements, which we then use to infer gene regulatory networks that underpin lineage commitment events. Finally we show that these models can be used to generate in silico predictions of the effect of transcription factor perturbations. We validate this experimentally by showing that Brachyury is essential for the differentiation of neuromesodermal progenitors to somitic mesoderm fate by priming cis-regulatory elements.

<p align="center"> 
<img src="images/overview_github.png" width="900" height="400"/>
</p>


Content
-------
* `/acc/`: analysis of chromatin accessibility data
* `/rna/`: analysis of RNA expression data
* `/accrna/`: simultaneous analysis of RNA expression and chromatin accessibility data

Snakemake pipeline
-------
We provide snakemake pipelines that can be used to reproduce many results. 
* `/rna/snakemake`: snakemake pipeline for RNA expression
* `/atac/ArchR/snakemake`: snakemake pipeline for chromatin accessibility using ArchR
* `/rna_atac/snakemake`: snakemake pipeline to integrate RNA expression and chromatin accessibility results (MOFA, in silico ChIP-seq, etc.)

Please note that the snakemake pipeline is rather complex and needs to be simplified and polished. It is currently useful to get an idea of the pipeline, but bare in mind that it  won't work straight away.

IGV Genome browser session
-------
We provide a precomputed IGV Genome Browser Session that can be used to interactively explore the ATAC-seq profiles, as shown in the screenshot below:

<p align="center"> 
<img src="images/igv_screenshot_github.png" width="650" height="350"/>
</p>

It can be downloaded running the following command line:
```
wget ftp://ftpusr92:5FqIACU9@ftp1.babraham.ac.uk/igv_session_celltype.tar.gz
```

Then load the file `igv_session.xml` using `File -> Open Session`.

<!-- The following [videotutorial](XXX) shows how to download and load the IGV session -->

R Shiny app
-------
The R shiny app for interactive data analaysis is available [here](https://www.bioinformatics.babraham.ac.uk/shiny/shiny_multiome_organogenesis/)

<!-- Pre-recorded talk
-------
This precorded talk by Ricard Argelaguet presents an overview of the study. -->

Directories
-------
* Mapping to the reference atlas: `/rna/mapping`
* MOFA dimensionality reduction: `/rna_atac/mofa`
* Analysis of gene markers: `/rna_atac/rna_vs_acc/pseudobulk/gene_markers_rna_vs_acc`
* in silico ChIP-seq: `/rna_atac/virtual_chipseq_library`
* Metacell inference: `/rna/metacells/run`
* Catalogue of TF activities per cell type (Figure 3): `/rna_atac/rna_vs_chromvar_chip/pseudobulk/per_celltype`
* Gene regulatory network of NMP differentiation (Figure 4): `/rna_atac/gene_regulatory_networks/metacells/trajectories`

Data
----
<!-- The raw data is accessible at GEO ([XXXX](XXXX)).  -->
The data can be downloaded from the following FTP server: 
```
Hostname 	ftp1.babraham.ac.uk
Username 	ftpusr92
Password 	5FqIACU9
FTP URL 	ftp://ftpusr92:5FqIACU9@ftp1.babraham.ac.uk
```

Directory structure:

- `sample_metadata.txt.gz`: cell metadata file
- `results`: results folder
	- `rna`: results based on RNA expression alone
	- `atac`: results based on chromatin accessibility alone
	- `rna_atac`: results based on both RNA expression and chromatin accessibility
- `data`: data folder
	- `original`: CellRanger output files
	- `processed`: processed data objects
		- `rna`: Seurat, anndata and SingleCellExperiment objects.
		- `atac/archR`: ArchR objects
- `igv_session_celltype.tar.gz`: IGV session of celltype-specific ATAC profiles
- `igv_session_brachyury_ko.tar.gz`: IGV session of celltype-specific ATAC profiles for the Brachyury KO study

To download a specific file:
```
wget ftp://ftpusr92:5FqIACU9@ftp1.babraham.ac.uk/data/processed/rna/SingleCellExperiment.rds .
```

To download everything (~60GB):
```
wget -r ftp://ftpusr92:5FqIACU9@ftp1.babraham.ac.uk/ .
```

If your download from the FTP server is slow, we also provide a temporary download via Dropbox [here](https://www.dropbox.com/sh/y4drtqi82vwl8vf/AAAsLcrye8jUTm1XPv7VNYhFa?dl=0)


Twitter thread
--------
We all know this is the most important day in the road towards a scientific publication:  

[![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/bukotsunikki.svg?style=social)](https://twitter.com/RArgelaguet/status/1537146799772815366)

Contact
-------
We have created a Slack group to discuss results, questions, collaborations, etc. related to the study. Feel free to drop by [using this link](https://join.slack.com/t/mouseembryo10-waq1273/shared_invite/zt-1dxn064kk-garRxOLAhLOUFNZBwqzfqQ). Alternatively, feel free to reach me via email at rargelaguet@altoslabs.com

