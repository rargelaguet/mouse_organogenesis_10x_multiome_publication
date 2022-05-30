here::i_am("atac/archR/differential/metacells/celltype/analysis/define_markers.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--differential_results',        type="character",     help='')
p$add_argument('--differential_stats',        type="character",     help='')
p$add_argument('--motif_annotation',        type="character",     help='')
p$add_argument('--min_score',       type="double",       default=0.75,      help='Minimum marker score')
p$add_argument('--min_diff',       type="double",       default=1.0,      help='Minimum marker score')
p$add_argument('--fdr',       type="double",       default=0.01,      help='False discovery rate')
p$add_argument('--outdir',          type="character",     help='Output directory')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# io$basedir <- file.path(io$basedir,"test")
# args <- list()
# args$motif_annotation <- "CISBP"
# args$differential_results <- file.path(io$basedir,sprintf("results/atac/archR/chromvar_chip/metacells/differential/celltypes/%s",args$motif_annotation))
# args$differential_stats <- file.path(args$differential_results,"diff_stats.txt")
# args$min_diff <- 10
# args$min_score <- 0.85
# args$fdr <- 0.01
# args$outdir <- file.path(io$basedir,sprintf("results/atac/archR/chromvar_chip/metacells/differential/celltypes/%s/markers",args$motif_annotation)); dir.create(args$outdir, showWarnings = F)
## END TEST ##

dir.create(args$outdir, showWarnings=F, recursive=T)

##################
## Load results ##
##################

# i <- opts$celltypes[1]; j <- opts$celltypes[2]
diff.dt <- opts$celltypes %>% map(function(i) { opts$celltypes %>% map(function(j) {
  file <- sprintf("%s/chromVAR_%s_vs_%s.txt.gz", args$differential_results,i,j)
  if (file.exists(file)) {
    fread(file, select=c(1,2,3)) %>% 
      .[,c("celltypeA","celltypeB"):=list(as.factor(i),as.factor(j))] %>% 
      # setnames("logFC","diff") %>%
      return
  }
}) %>% rbindlist }) %>% rbindlist %>% 
  .[is.na(diff),c("diff","padj"):=list(0,1)] %>%
  .[,sig:=ifelse(abs(diff)>=args$min_diff & padj<=args$fdr,TRUE,FALSE)] %>%
  .[,direction:=c("up","down")[as.numeric(diff>0)+1]]  # up = higher accessibility in celltype A

diff.dt %>% setnames("gene","tf")

# hist(diff.dt$diff)

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
## Define marker tfs ##
#########################

ncelltypes <- unique(c(as.character(unique(diff.dt$celltypeA)),as.character(unique(diff.dt$celltypeB)))) %>% length

foo <- diff.dt[,.(N=.N, diff=mean(diff), score=sum(sig==T & direction=="up")), by=c("celltypeA","tf")] %>% setnames("celltypeA","celltype")
bar <- diff.dt[,.(N=.N, diff=mean(diff), score=sum(sig==T & direction=="down")), by=c("celltypeB","tf")] %>% setnames("celltypeB","celltype")
  
markers.dt <- merge(foo[,c("celltype","tf","score")],bar[,c("celltype","tf","score")],by=c("celltype","tf"), all=TRUE) %>% 
  .[is.na(score.x),score.x:=0] %>% .[is.na(score.y),score.y:=0] %>%
  .[,score:=score.x+score.y] %>%
  .[,c("score.x","score.y"):=NULL] %>%
  .[,score:=round(score/(ncelltypes-1),2)] %>%
  setorder(celltype,-score)
# rm(foo,bar)

stopifnot(max(markers.dt$score,na.rm=T)==1)

#########################
## Add MeanDiff values ##
#########################

tmp <- merge(foo[,c("celltype","tf","N","diff")], bar[,c("celltype","tf","N","diff")], by=c("celltype","tf"), all=TRUE) %>% 
  .[is.na(N.x),N.x:=0] %>% .[is.na(N.y),N.y:=0] %>%
  .[,diff:=weighted.mean(c(-diff.x,diff.y), c(N.x,N.y)),by=c("celltype","tf")]

markers.dt <- markers.dt %>% 
  merge(tmp[,c("celltype","tf","diff")], by=c("celltype","tf"))

##########
## Save ##
##########

# Save marker score for all combination of tfs and cell types
length(unique(markers.dt$tf))
length(unique(markers.dt$celltype))
fwrite(markers.dt, file.path(args$outdir,"chromVAR_markers_all.txt.gz"), sep="\t")

# Save marker score for strong markers
markers_filt.dt <- markers.dt %>% .[score>=args$min_score & diff>=args$min_diff]
length(unique(markers_filt.dt$tf))
length(unique(markers_filt.dt$celltype))
table(markers_filt.dt$celltype)
fwrite(markers_filt.dt, file.path(args$outdir,"chromVAR_markers_filt.txt.gz"), sep="\t")

##########
## TEST ##
##########

# diff.dt[tf=="DMBX1" & celltypeA=="Epiblast"]
# diff.dt[tf=="DMBX1" & celltypeB=="Epiblast"]
