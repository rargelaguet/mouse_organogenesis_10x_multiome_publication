here::i_am("atac/archR/differential/cells/celltype/analysis/define_markers.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--differential_results',        type="character",     help='')
p$add_argument('--differential_stats',        type="character",     help='')
p$add_argument('--matrix',        type="character",     help='')
p$add_argument('--min_score',       type="double",       default=0.75,      help='Minimum marker score')
p$add_argument('--min_diff',       type="double",       default=1.0,      help='Minimum marker score')
p$add_argument('--fdr',       type="double",       default=0.01,      help='False discovery rate')
p$add_argument('--outdir',          type="character",     help='Output directory')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# io$basedir <- file.path(io$basedir,"test")
# args <- list()
# args$matrix <- "GeneScoreMatrix_TSS" # "PeakMatrix"
# args$differential_results <- file.path(io$basedir,sprintf("results/atac/archR/differential/cells/celltype/%s",args$matrix))
# args$differential_stats <- file.path(args$differential_results,"diff_stats.txt")
# args$min_diff <- 0.15
# args$min_score <- 0.85
# args$fdr <- 0.01
# args$outdir <- file.path(io$basedir,sprintf("results/atac/archR/differential/cells/celltype/%s/markers",args$matrix)); dir.create(args$outdir, showWarnings = F)
## END TEST ##

print(args)

dir.create(args$outdir, showWarnings=F, recursive=T)

#####################
## Define settings ##
#####################

# Annoying stuff
if (grepl("Peak",args$matrix)) {
  select.cols <- c(1,2,3)
} else if (grepl("Gene",args$matrix)) {
  select.cols <- c(1,3,4)
} else {
  stop()
}

##################
## Load results ##
##################

# i <- opts$celltypes[2]; j <- opts$celltypes[1]
diff.dt <- opts$celltypes %>% map(function(i) { opts$celltypes %>% map(function(j) {
  file <- sprintf("%s/%s_%s_vs_%s.txt.gz", args$differential_results,args$matrix,i,j)
  if (file.exists(file)) {
    fread(file, select=select.cols) %>% 
      setnames(c("feature","diff","padj_fdr")) %>%
      .[,c("celltypeA","celltypeB"):=list(as.factor(i),as.factor(j))] %>% 
      return
  }
}) %>% rbindlist }) %>% rbindlist %>% 
  .[is.na(diff),c("diff","padj_fdr"):=list(0,1)] %>%
  .[,sig:=ifelse(abs(diff)>=args$min_diff & padj_fdr<=args$fdr,TRUE,FALSE)] %>%
  .[,direction:=c("up","down")[as.numeric(diff>0)+1]]  # up = higher accessibility in celltype A

###################
## Sanity checks ##
###################

stopifnot(sum(is.na(diff.dt$sig))==0)

# Load stats
diff_stats.dt <- fread(args$differential_stats) %>% 
  setnames(c("celltypeA","celltypeB","N_groupA","N_groupB"))

# check if some DA comparison is missing
tmp <- diff_stats.dt %>% 
  merge(diff.dt[,c("celltypeA","celltypeB")] %>% unique %>% .[,done:=TRUE], all.x=TRUE, by=c("celltypeA","celltypeB")) %>%
  .[is.na(done),done:=FALSE]
stopifnot(tmp$done==TRUE)

#########################
## Define marker genes ##
#########################

ncelltypes <- unique(c(as.character(unique(diff.dt$celltypeA)),as.character(unique(diff.dt$celltypeB)))) %>% length

foo <- diff.dt[,.(N=.N, diff=mean(diff), score=sum(sig==T & direction=="up")), by=c("celltypeA","feature")] %>% setnames("celltypeA","celltype")
bar <- diff.dt[,.(N=.N, diff=mean(diff), score=sum(sig==T & direction=="down")), by=c("celltypeB","feature")] %>% setnames("celltypeB","celltype")
  
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

tmp <- merge(foo[,c("celltype","feature","N","diff")], bar[,c("celltype","feature","N","diff")], by=c("celltype","feature"), all=TRUE) %>% 
  .[is.na(N.x),N.x:=0] %>% .[is.na(N.y),N.y:=0] %>%
  .[,diff:=weighted.mean(c(-diff.x,diff.y), c(N.x,N.y)),by=c("celltype","feature")]

markers_peaks.dt <- markers_peaks.dt %>% 
  merge(tmp[,c("celltype","feature","diff")], by=c("celltype","feature"))

##########
## Save ##
##########

# Save marker score for all combination of genes and cell types
length(unique(markers_peaks.dt$feature))
length(unique(markers_peaks.dt$celltype))
fwrite(markers_peaks.dt, file.path(args$outdir,sprintf("%s_markers_all.txt.gz",args$matrix)), sep="\t")

# Save marker score for strong markers
markers_peaks_filt.dt <- markers_peaks.dt %>% .[score>=args$min_score & diff>=args$min_diff]
length(unique(markers_peaks_filt.dt$feature))
length(unique(markers_peaks_filt.dt$celltype))
table(markers_peaks_filt.dt$celltype)
fwrite(markers_peaks_filt.dt, file.path(args$outdir,sprintf("%s_markers_filt.txt.gz",args$matrix)), sep="\t")

##########
## TEST ##
##########

