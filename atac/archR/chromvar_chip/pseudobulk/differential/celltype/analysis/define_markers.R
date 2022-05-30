here::i_am("atac/archR/chromvar_chip/pseudobulk/differential/celltype/analysis/define_markers.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--differential_results',        type="character",     help='')
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
# args$differential_results <- file.path(io$basedir,sprintf("results/atac/archR/chromvar_chip/pseudobulk_with_replicates/differential/celltypes/%s/parsed/diff_results.txt.gz",args$motif_annotation))
# args$min_diff <- 15
# args$min_score <- 0.85
# args$fdr <- 0.01
# args$outdir <- file.path(io$basedir,sprintf("results/atac/archR/chromvar_chip/pseudobulk_with_replicates/differential/celltypes/%s/parsed",args$motif_annotation)); dir.create(args$outdir, showWarnings = F)
## END TEST ##

dir.create(args$outdir, showWarnings=F, recursive=T)

##################
## Load results ##
##################

diff.dt <- fread(args$differential_results) %>% 
  .[abs(diff)>=args$min_diff & padj<=args$fdr] %>%
  .[,celltypeA:=factor(celltypeA,levels=opts$celltypes)] %>%
  .[,celltypeB:=factor(celltypeB,levels=opts$celltypes)] %>%
  .[,direction:=as.factor(c("down","up"))[as.numeric(diff<0)+1]]

ncelltypes <- unique(c(as.character(unique(diff.dt$celltypeA)),as.character(unique(diff.dt$celltypeB)))) %>% length

####################
## Define markers ##
####################

ncelltypes <- unique(c(as.character(unique(diff.dt$celltypeA)),as.character(unique(diff.dt$celltypeB)))) %>% length

foo <- diff.dt[,.(N=.N, diff=mean(diff), score=sum(direction=="up")), by=c("celltypeA","gene")] %>% setnames("celltypeA","celltype")
bar <- diff.dt[,.(N=.N, diff=mean(diff), score=sum(direction=="down")), by=c("celltypeB","gene")] %>% setnames("celltypeB","celltype")
  
markers.dt <- merge(foo[,c("celltype","gene","score")], bar[,c("celltype","gene","score")], by=c("celltype","gene"), all=TRUE) %>% 
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

tmp <- merge(foo[,c("celltype","gene","N","diff")], bar[,c("celltype","gene","N","diff")], by=c("celltype","gene"), all=TRUE) %>% 
  .[is.na(N.x),N.x:=0] %>% .[is.na(N.y),N.y:=0] %>%
  .[,diff:=weighted.mean(c(-diff.x,diff.y), c(N.x,N.y)) %>% round(2),by=c("celltype","gene")]

markers.dt <- markers.dt %>% 
  merge(tmp[,c("celltype","gene","diff")], by=c("celltype","gene"))

##########
## Save ##
##########

# Save marker score for all combination of tfs and cell types
length(unique(markers.dt$gene))
length(unique(markers.dt$celltype))
fwrite(markers.dt, file.path(args$outdir,"markers_all.txt.gz"), sep="\t")

# Save marker score for strong markers
markers_filt.dt <- markers.dt %>% .[score>=args$min_score & diff>=args$min_diff]
length(unique(markers_filt.dt$gene))
length(unique(markers_filt.dt$celltype))
table(markers_filt.dt$celltype)
fwrite(markers_filt.dt, file.path(args$outdir,"markers_filt.txt.gz"), sep="\t")
