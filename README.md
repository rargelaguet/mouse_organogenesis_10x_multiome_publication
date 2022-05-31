# Decoding gene regulation in the mouse embryo using single-cell multi-omics

This repository contains the scripts to reproduce the results of the manuscript [Decoding gene regulation in the mouse embryo using single-cell multi-omics](XXXX). 


Abstract
--------
Following gastrulation, the three primary germ layers develop into the major organs in a process known as organogenesis. Single-cell RNA sequencing has enabled the profiling of the gene expression dynamics of these cell fate decisions, yet a comprehensive map of the interplay between transcription factors and cis-regulatory elements is lacking, as are the underlying gene regulatory networks. Here we generate a multi-omics atlas of mouse early organogenesis by simultaneously profiling gene expression and chromatin accessibility from tens of thousands of single cells. We develop a computational method to leverage the multi-modal readouts to predict transcription factor binding events in cis-regulatory elements, which we then use to infer gene regulatory networks that underpin lineage commitment events. Finally we show that these models can be used to generate in silico predictions of the effect of transcription factor perturbations. We validate this experimentally by showing that Brachyury is essential for the differentiation of neuromesodermal progenitors to somitic mesoderm fate by priming cis-regulatory elements.

<!-- <p align="center"> 
<img src="images/figure1.png" width="650" height="350"/>
</p> -->


Twitter thread
--------
XXX

Content
-------
* `/acc/`: analysis of chromatin accessibility data
* `/rna/`: analysis of RNA expression data
* `/accrna/`: simultaneous analysis of RNA expression and chromatin accessibility data

Snakemake pipeline
-------
We provide snakemake pipelines that can be used to reproduce many results. 
* `/rna/snakemake`: snakemake pipeline for RNA expression
* `/atac/ArchR`: snakemake pipeline for chromatin accessibility using ArchR
* `/rna_atac`: snakemake pipeline to integrate RNA expression and chromatin accessibility results (MOFA, in silico ChIP-seq, etc.)

IGV Genome browser session
-------
We provide a precomputed IGV Genome Browser Session [here](XXX), which you can use to interactively explore the ATAC-seq profiles, as shown in the screenshot below:

<!-- <p align="center"> 
<img src="images/igv_screenshot.png" width="650" height="350"/>
</p> -->

Shiny app
-------
The shiny app is in preparation...

Shortcuts
-------
* MOFA dimensionality reduction (Figure 1c): XXX
* Analysis of gene markers (Figure 1d-f): XXX
* in silico ChIP-seq (Figure 2): XXX
* Metacell inference: XXX
* Catalogue of TF activities per cell type (Figure 3): XXX
* Gene regulatory network of NMP differentiation (Figure 4): XXX


Data
----
The raw data is accessible at GEO ([XXXX](XXXX)). 
The parsed data can be downloaded [XXX](XXXXX)

Contact
-------
* Ricard Argelaguet (ricard.argelaguet@gmail.com), [![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/bukotsunikki.svg?style=social)](https://twitter.com/RArgelaguet)

