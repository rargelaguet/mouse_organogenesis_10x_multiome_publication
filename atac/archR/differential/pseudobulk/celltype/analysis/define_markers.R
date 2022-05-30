here::i_am("atac/archR/differential/pseudobulk/celltype/analysis/define_markers.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--differential_results',        type="character",     help='')
p$add_argument('--matrix',        type="character",     help='')
p$add_argument('--min_score',       type="double",       default=0.75,      help='Minimum marker score')
p$add_argument('--min_fold_change',       type="double",       default=1.0,      help='Minimum marker score')
p$add_argument('--fdr',       type="double",       default=0.01,      help='False discovery rate')
p$add_argument('--outdir',          type="character",     help='Output directory')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# io$basedir <- file.path(io$basedir,"test")
# args <- list()
# args$matrix <- "PeakMatrix"
# args$differential_results <- file.path(io$basedir,sprintf("results/atac/archR/differential/pseudobulk/celltype/%s/parsed/diff_results.txt.gz",args$matrix))
# args$min_fold_change <- 1.0
# args$min_score <- 0.85
# args$fdr <- 0.01
# args$outdir <- file.path(io$basedir,sprintf("results/atac/archR/differential/pseudobulk/celltype/%s/parsed",args$matrix))
## END TEST ##

dir.create(args$outdir, showWarnings=F, recursive=T)

##################
## Load results ##
##################

# i <- opts$celltypes[1]; j <- opts$celltypes[2]
# diff_results.dt <- opts$celltypes %>% map(function(i) { opts$celltypes %>% map(function(j) {
#   file <- sprintf("%s/%s_%s_vs_%s.txt.gz", args$differential_results,args$matrix,i,j)
#   if (file.exists(file)) {
#     fread(file, select=c(1,2,3)) %>% 
#       .[,c("celltypeA","celltypeB"):=list(as.factor(i),as.factor(j))] %>% 
#       setnames("logFC","diff") %>%
#       return
#   }
# }) %>% rbindlist }) %>% rbindlist %>% 
#   .[is.na(diff),c("diff","padj_fdr"):=list(0,1)] %>%
#   .[,sig:=ifelse(abs(diff)>=args$min_fold_change & padj_fdr<=args$fdr,TRUE,FALSE)] %>%
#   .[,direction:=c("up","down")[as.numeric(diff>0)+1]]  # up = higher accessibility in celltype A

diff_results.dt <- fread(args$differential_results) %>% 
  .[abs(logFC)>=args$min_fold_change & padj_fdr<=args$fdr] %>%
  .[,celltypeA:=factor(celltypeA,levels=opts$celltypes)] %>%
  .[,celltypeB:=factor(celltypeB,levels=opts$celltypes)] %>%
  .[,direction:=as.factor(c("down","up"))[as.numeric(logFC<0)+1]]

ncelltypes <- unique(c(as.character(unique(diff_results.dt$celltypeA)),as.character(unique(diff_results.dt$celltypeB)))) %>% length

#########################
## Define marker genes ##
#########################

foo <- diff_results.dt[,.(N=.N, logFC=mean(logFC), score=sum(direction=="up")), by=c("celltypeA","feature")] %>% setnames("celltypeA","celltype")
bar <- diff_results.dt[,.(N=.N, logFC=mean(logFC), score=sum(direction=="down")), by=c("celltypeB","feature")] %>% setnames("celltypeB","celltype")

markers_peaks.dt <- merge(foo[,c("celltype","feature","score")],bar[,c("celltype","feature","score")],by=c("celltype","feature"), all=TRUE) %>% 
  .[is.na(score.x),score.x:=0] %>% .[is.na(score.y),score.y:=0] %>%
  .[,score:=score.x+score.y] %>%
  .[,c("score.x","score.y"):=NULL] %>%
  .[,score:=round(score/(ncelltypes-1),2)] %>%
  setorder(celltype,-score)
# rm(foo,bar)

stopifnot(max(markers_peaks.dt$score,na.rm=T)==1)

#########################
## Add MeanDiff values ##
#########################

tmp <- merge(foo[,c("celltype","feature","N","logFC")], bar[,c("celltype","feature","N","logFC")], by=c("celltype","feature"), all=TRUE) %>% 
  .[is.na(N.x),N.x:=0] %>% .[is.na(N.y),N.y:=0] %>%
  .[,logFC:=weighted.mean(c(-logFC.x,logFC.y), c(N.x,N.y)) %>% round(2),by=c("celltype","feature")]

markers_peaks.dt <- markers_peaks.dt %>% 
  merge(tmp[,c("celltype","feature","logFC")], by=c("celltype","feature"))

##########
## Save ##
##########

# Save marker score for all combination of genes and cell types
length(unique(markers_peaks.dt$feature))
length(unique(markers_peaks.dt$celltype))
fwrite(markers_peaks.dt, file.path(args$outdir,"markers_all.txt.gz"), sep="\t")

# Save marker score for strong markers
markers_peaks_filt.dt <- markers_peaks.dt %>% .[score>=args$min_score & logFC>=args$min_fold_change]
length(unique(markers_peaks_filt.dt$feature))
length(unique(markers_peaks_filt.dt$celltype))
table(markers_peaks_filt.dt$celltype)
fwrite(markers_peaks_filt.dt, file.path(args$outdir,"markers_filt.txt.gz"), sep="\t")
