here::i_am("rna_atac/virtual_chipseq_library/metacells/create_virtual_chipseq_library_metacells.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(furrr))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--tf2peak_cor',  type="character",              help='Correlations between TF RNA expression and peak accessibility') 
p$add_argument('--atac_peak_matrix',  type="character",              help='ATAC Peak matrix (pseudobulk)') 
p$add_argument('--peak_metadata',             type="character",        help='Peak metadata file')
p$add_argument('--motif2gene',  type="character",              help='Motif annotation')
p$add_argument('--motifmatcher',  type="character",              help='Motif annotation') 
p$add_argument('--motif_annotation',  type="character",              help='Motif annotation') 
# p$add_argument('--max_acc_score',     type="double",    default=1.5,    help='Maximum peak accessibility')
p$add_argument('--min_number_peaks',     type="integer",    default=50,    help='Minimum number of peaks per TF')
# p$add_argument('--min_peak_score',     type="integer",    default=15,    help='Minimum peak score')
p$add_argument('--outdir',          type="character",                help='Output directory')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')
p$add_argument('--test',  action="store_true",  help='Test mode?')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args <- list()
# args$motif_annotation <- "CISBP"
# args$motifmatcher <- file.path(io$basedir,sprintf("processed/atac/archR/Annotations/%s-Scores.rds",args$motif_annotation))
# # args$atac_peak_matrix <- file.path(io$basedir,"results/atac/archR/metacells/all_cells/PeakMatrix/PeakMatrix_summarized_experiment_metacells.rds")
# args$atac_peak_matrix <- file.path(io$basedir,"results/atac/archR/metacells/trajectories/nmp/PeakMatrix/PeakMatrix_summarized_experiment_metacells.rds")
# args$peak_metadata <- file.path(io$basedir,"processed/atac/archR/PeakCalls/peak_metadata.tsv.gz")
# args$motif2gene <- file.path(io$basedir,sprintf("processed/atac/archR/Annotations/%s_motif2gene.txt.gz",args$motif_annotation))
# # args$tf2peak_cor <- file.path(io$basedir,sprintf("results/rna_atac/rna_vs_acc/metacells/TFexpr_vs_peakAcc/%s_cor_TFexpr_vs_peakAcc.rds",args$motif_annotation))
# args$tf2peak_cor <- file.path(io$basedir,sprintf("results/rna_atac/rna_vs_acc/metacells/trajectories/nmp/TFexpr_vs_peakAcc/%s_cor_TFexpr_vs_peakAcc.rds",args$motif_annotation))
# # args$max_acc_score <- 4.5
# # args$min_peak_score <- 15
# args$min_number_peaks <- 50
# # args$outdir <- file.path(io$basedir,"results/rna_atac/virtual_chipseq/metacells/JASPAR/test")
# args$outdir <- file.path(io$basedir,sprintf("results/rna_atac/virtual_chipseq/metacells/trajectories/nmp/%s/test",args$motif_annotation))
# args$threads <- 1
# args$test <- FALSE
## END TEST ##

#####################
## Define settings ##
#####################

# I/O
dir.create(args$outdir, showWarnings = F, recursive = T)

# Parallel processing
# plan(sequential)
plan(multisession, workers=args$threads)

###############################
## Load pseudobulk ATAC data ##
###############################

# Load ATAC SummarizedExperiment
atac_peak_matrix_metacells.se <- readRDS(args$atac_peak_matrix)

# Normalise
assayNames(atac_peak_matrix_metacells.se) <- "counts"
assay(atac_peak_matrix_metacells.se,"logcounts") <- log2(1e6*(sweep(assay(atac_peak_matrix_metacells.se),2,colSums(assay(atac_peak_matrix_metacells.se)),"/"))+0.5)

# Load peak metadata
peak_metadata.dt <- fread(args$peak_metadata) %>% 
  # .[score>=args$min_peak_score] %>%
  .[,idx:=sprintf("%s:%s-%s",chr,start,end)]

######################
## Load motifmatchR ##
######################

motifmatcher.se <- readRDS(args$motifmatcher)

# Subset peaks
stopifnot(sort(rownames(motifmatcher.se))==sort(rownames(atac_peak_matrix_metacells.se)))
motifmatcher.se <- motifmatcher.se[rownames(atac_peak_matrix_metacells.se),]

###############################
## Load TF2peak correlations ##
###############################

tf2peak_cor.se <- readRDS(args$tf2peak_cor)

stopifnot(rownames(tf2peak_cor.se)%in%rownames(motifmatcher.se))
stopifnot(rownames(tf2peak_cor.se)%in%rownames(atac_peak_matrix_metacells.se))

##################
## Subset peaks ##
##################

peaks <- Reduce(intersect,list(rownames(motifmatcher.se),rownames(tf2peak_cor.se),peak_metadata.dt$idx))

print(sprintf("Number of peaks: %s",length(peaks)))

# motifmatcher.se <- motifmatcher.se[peaks,]
tf2peak_cor.se <- tf2peak_cor.se[peaks,]
atac_peak_matrix_metacells.se <- atac_peak_matrix_metacells.se[peaks,]
peak_metadata.dt <- peak_metadata.dt[idx%in%peaks]

################################
## Load motif2gene annotation ##
################################

motif2gene.dt <- fread(args$motif2gene)

# Manually remove some motifs
motif2gene.dt <- motif2gene.dt[!motif%in%c("T_789")]

# opts$motif_annotation <- args$motif_annotation
# source(here::here("atac/archR/load_motif_annotation.R"))

# (TO FIX) Make sure that there is a one-to-one match between TFs and genes
# tmp <- motif2gene.dt[,.N,by=c("gene")] %>% .[N>1]
# motif2gene.dt[gene=="TFAP2B"]
# motif2gene.dt[gene=="TFAP2A"]
# foo <- motif2gene.dt[,N:=.N,by=c("gene")] %>% .[N==1] %>% .[,N:=NULL] 
# bar <- motif2gene.dt[,N:=.N,by=c("gene")] %>% .[N>1] %>% .[,N:=NULL] %>% .[,head(.SD,n=1), by="gene"]
# motif2gene.dt <- rbind(foo,bar)

################
## Rename TFs ##
################

motifs <- intersect(colnames(motifmatcher.se),motif2gene.dt$motif)
TFs <- intersect(colnames(tf2peak_cor.se),motif2gene.dt$gene)

motif2gene_filt.dt <- motif2gene.dt[motif%in%motifs & gene%in%TFs]
motifs <- motif2gene_filt.dt$motif
TFs <- motif2gene_filt.dt$gene

tmp <- TFs; names(tmp) <- motifs

stopifnot(motif2gene_filt.dt$motif%in%colnames(motifmatcher.se))
stopifnot(motif2gene_filt.dt$gene%in%colnames(tf2peak_cor.se))

motifmatcher.se <- motifmatcher.se[,motifs]
colnames(motifmatcher.se) <- tmp[colnames(motifmatcher.se)]
tf2peak_cor.se <- tf2peak_cor.se[,TFs]
stopifnot(colnames(motifmatcher.se)==colnames(tf2peak_cor.se))

##################
## Prepare data ##
##################

tf2peak_cor.mtx <- assay(tf2peak_cor.se,"cor")
motifmatcher.mtx <- assay(motifmatcher.se,"motifScores")
atac.mtx <- assay(atac_peak_matrix_metacells.se[peaks,],"logcounts") %>% round(3)

######################################
## Create virtual chip-seq library ##
######################################

print("Predicting TF binding sites...")

# opts$TFs <- c("T", "ZIC2", "SOX2", "TAL1", "GATA1", "FOXA2", "GATA4", "NKX2-5")
opts$TFs <- TFs
if (args$test) {
  print("Testing mode activated...")
  opts$TFs <- opts$TFs %>% head(n=3)
} 
stopifnot(!duplicated(opts$TFs))
print(sprintf("Number of TFs: %s",length(opts$TFs)))

# i <- "GATA1"
# virtual_chip.dt <- opts$TFs %>% future_map(function(i) {
virtual_chip.dt <- opts$TFs %>% map(function(i) {
  print(sprintf("%s (%s/%s)",i,match(i,opts$TFs),length(opts$TFs)))
  
  # motif <- motif2gene.dt[gene==i,motif]
  # stopifnot(length(motif)==1)

  # peaks <- names(which(abs(assay(tf2peak_cor.se[peaks,i],"cor")[,1])>=0))
  peaks <- names(which(abs(tf2peak_cor.mtx[,i])>0)) # we only consider chromatin activators
  
  if (length(peaks)>=args$min_number_peaks) {
    
    # calculate accessibility score
    # max_accessibility_score <- apply(assay(atac_peak_matrix_metacells.se[peaks,]),1,max) %>% round(2)
    max_accessibility_score <- apply(atac.mtx[peaks,],1,max) %>% round(2)
    # max_accessibility_score[max_accessibility_score>=args$max_acc_score] <- args$max_acc_score
    
    # calculate correlation score
    # correlation_score <- assay(tf2peak_cor.se[,i],"cor")[peaks,1] %>% round(2)
    correlation_score <- tf2peak_cor.mtx[peaks,i] %>% round(2)
    correlation_score[correlation_score==0] <- NA
    
    # calculate motif score
    # motif_score <- assay(motifmatcher.se[peaks,motif2gene.dt[gene==i,motif]],"motifScores")[,1]
    motif_score <- motifmatcher.mtx[peaks,i]
    motif_score <- round(motif_score/max(motif_score),2)
    
    # calculate motif counts
    # motif_counts <- assay(motifmatcher.se[peaks,i],"motifCounts")[,1] %>% round(2)
    
    # predicted_score <- minmax.normalisation(max_accessibility_score * correlation_score * motif_score * motif_counts) 
    predicted_score <- correlation_score * minmax.normalisation(max_accessibility_score * motif_score)
    # predicted_score <- predicted_score/max(predicted_score,na.rm=T) %>% round(2)
    
    tmp <- data.table(
      peak = peaks, 
      correlation_score = correlation_score,
      max_accessibility_score = max_accessibility_score,
      motif_score = motif_score,
      # motif_counts = motif_counts,
      score = round(predicted_score,2)
    ) %>% sort.abs("score") %>% 
      # .[motif_score==0,motif_score:=NA] %>%
      .[,c("peak","score","correlation_score","max_accessibility_score","motif_score")]
    
    bed.dt <- tmp %>% copy %>%
      .[,peak:=str_replace(peak,":","-")] %>%
      .[,chr:=strsplit(peak,"-") %>% map_chr(1)] %>%
      .[,start:=strsplit(peak,"-") %>% map_chr(2)] %>%
      .[,end:=strsplit(peak,"-") %>% map_chr(3)] %>%
      .[,c("chr","start","end","score")]
    
    # Save 
    fwrite(tmp, sprintf("%s/%s.txt.gz",args$outdir,i), sep="\t", quote=F, col.names = T)
    fwrite(bed.dt, sprintf("%s/%s.bed.gz",args$outdir,i), sep="\t", quote=F, col.names = F)
    
    to_return.dt <- tmp[!is.na(score),c("peak","score")] %>% .[,tf:=i]
    return(to_return.dt)
  }
}) %>% rbindlist

# Load precomputed
# virtual_chip.dt <- opts$TFs %>% map(function(i) {
#     print(sprintf("%s (%s/%s)",i,match(i,opts$TFs),length(opts$TFs)))
#     file <- sprintf("%s/%s.txt.gz",args$outdir,i)
#     if (file.exists(file)) {
#       fread(file) %>% .[!is.na(score),c("peak","score")] %>% .[,tf:=i] %>% return
#     }
#   }) %>% rbindlist

########################
## Save sparse matrix ##
########################

virtual_chip.mtx <- virtual_chip.dt %>% 
  .[,peak:=factor(peak,levels=rownames(motifmatcher.se))] %>%
  data.table::dcast(peak~tf, value.var="score", fill=0, drop=F) %>%
  matrix.please %>% Matrix::Matrix(.)

saveRDS(virtual_chip.mtx, file.path(args$outdir,"virtual_chip_matrix.rds"))

###################################################################
## Update motifmatchr results using the virtual ChIP-seq library ##
###################################################################

print("Updating motifmatchr results using the virtual ChIP-seq library...")

motifmatcher_chip.se <- motifmatcher.se[,colnames(virtual_chip.mtx)]
# assays(motifmatcher_chip.se) <- assays(motifmatcher_chip.se)[c("motifScores")]

# reset motif matches
stopifnot(rownames(virtual_chip.mtx)==rownames(motifmatcher_chip.se))
stopifnot(colnames(virtual_chip.mtx)==colnames(motifmatcher_chip.se))
assay(motifmatcher_chip.se,"VirtualChipScores") <- virtual_chip.mtx
saveRDS(motifmatcher_chip.se, file.path(args$outdir,"motifmatchr_virtual_chip.rds"))
