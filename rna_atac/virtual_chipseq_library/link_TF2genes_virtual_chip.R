here::i_am("rna_atac/virtual_chipseq_library/link_TF2genes_virtual_chip.R")

# load default setings
source(here::here("settings.R"))
source(here::here("utils.R"))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--motif_annotation',  type="character",              help='Motif annotation') 
p$add_argument('--virtual_chip_matrix',       type="character",                help='TF2peak virtual chip matrix')
p$add_argument('--peak2gene',       type="character",                help='Links between peaks and genes based on genomic distance')
p$add_argument('--distance',  type="integer",            default=1e5,      help='Maximum distance for a linkage between a peak and a gene')
p$add_argument('--min_chip_score',  type="double",            default=0.25,      help='Minimum in silico ChIP-seq score')
p$add_argument('--outfile',       type="character",                help='Output file')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST
# args$motif_annotation <- "CISBP"
# args$virtual_chip_matrix <- file.path(io$basedir,sprintf("results/rna_atac/virtual_chipseq/pseudobulk/%s/virtual_chip_matrix.rds",args$motif_annotation))
# args$peak2gene <- file.path(io$basedir,"results/atac/archR/peak_calling/peaks2genes/peaks2genes_all.txt.gz")
# args$distance <- 5e4
# args$min_chip_score <- 0.25
# args$outfile <- file.path(io$basedir,sprintf("results/rna_atac/virtual_chipseq/pseudobulk/%s/TF2gene_after_virtual_chip.txt.gz",args$motif_annotation))
## END TEST

## START TEST
# io$basedir <- file.path(io$basedir,"test")
# args$motif_annotation <- "CISBP"
# args$virtual_chip_matrix <- file.path(io$basedir,sprintf("results/rna_atac/virtual_chipseq/metacells/trajectories/nmp/%s/virtual_chip_matrix.rds",args$motif_annotation))
# args$peak2gene <- file.path(io$basedir,"results/atac/archR/peak_calling/peaks2genes/peaks2genes_all.txt.gz")
# args$distance <- 5e4
# args$min_chip_score <- 0.15
# args$outfile <- file.path(io$basedir,sprintf("results/rna_atac/virtual_chipseq/metacells/trajectories/nmp/%s/TF2gene_after_virtual_chip.txt.gz",args$motif_annotation))
## END TEST

## START TEST
# io$basedir <- file.path(io$basedir,"test")
# args$motif_annotation <- "CISBP"
# args$virtual_chip_matrix <- file.path(io$basedir,sprintf("results/rna_atac/virtual_chipseq/metacells/%s/virtual_chip_matrix.rds",args$motif_annotation))
# args$peak2gene <- file.path(io$basedir,"results/atac/archR/peak_calling/peaks2genes/peaks2genes_all.txt.gz")
# args$distance <- 5e4
# args$min_chip_score <- 0.15
# args$outfile <- file.path(io$basedir,sprintf("results/rna_atac/virtual_chipseq/metacells/%s/TF2gene_after_virtual_chip.txt.gz",args$motif_annotation))
## END TEST

###########################
## Load virtual ChIP-seq ##
###########################

virtual_chip.mtx <- readRDS(args$virtual_chip_matrix)

# unique(as.numeric(virtual_chip.mtx[1:1000,1:100]))

#########################################################
## Load peak2gene linkages using only genomic distance ##
#########################################################

peak2gene.dt <- fread(args$peak2gene) %>%
  # peak2gene.dt <- fread(io$archR.peak2gene.nearest) %>%
  .[dist<=args$distance] %>%
  .[,peak:=sprintf("chr%s:%s-%s",chr,peak.start,peak.end)]

# Sanity checks
stopifnot(length(intersect(rownames(virtual_chip.mtx),unique(peak2gene.dt$peak)))>1e5)

#########################################################
## Link TFs to target genes using the virtual ChIP-seq ##
#########################################################

tf2gene_chip.dt <- colnames(virtual_chip.mtx) %>% map(function(i) {
  print(i)
  # Select target peaks (note that we only take positive correlations into account)
  target_peaks_i <- names(which(virtual_chip.mtx[,i]>=args$min_chip_score))
  
  if (length(target_peaks_i)>=1) {
    
    tmp <- data.table(
      tf = i,
      peak = target_peaks_i,
      chip_score = virtual_chip.mtx[target_peaks_i,i]
    ) %>% merge(peak2gene.dt[peak%in%target_peaks_i,c("peak","gene","dist")], by="peak")
    
    # print(sprintf("%s: %s target peaks & %s target genes",i,length(target_peaks_i),length(unique(tmp$gene))))
    return(tmp)
  }
}) %>% rbindlist

# Save
fwrite(tf2gene_chip.dt, args$outfile, sep="\t")
