here::i_am("rna/differential/metacells/celltype/analysis/define_marker_genes.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--differential_results',        type="character",     help='')
p$add_argument('--differential_stats',        type="character",     help='')
p$add_argument('--min_score',       type="double",       default=0.75,      help='Minimum marker score')
p$add_argument('--min_fold_change',       type="double",       default=1,      help='Minimum marker score')
p$add_argument('--fdr',       type="double",       default=0.01,      help='False discovery rate')
p$add_argument('--outdir',          type="character",     help='Output directory')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# io$basedir <- file.path(io$basedir,"test")
# args <- list()
# args$differential_results <- file.path(io$basedir,"results/rna/differential/metacells/celltype/parsed/diff_expr_results.txt.gz")
# args$differential_stats <- file.path(io$basedir,"results/rna/differential/metacells/celltype/parsed/diff_expr_stats.txt.gz")
# args$min_score <- 0.75
# args$min_fold_change <- 1.0
# args$fdr <- 0.01
# args$outdir <- file.path(io$basedir,"results/rna/differential/metacells/celltype/parsed"); dir.create(args$outdir, showWarnings = F)
## END TEST ##

dir.create(args$outdir, showWarnings = F)

###############################
## Load differential results ##
###############################

# i <- "Gut"; j <- "NMP"
# diff_results.dt <- opts$celltypes %>% map(function(i) { opts$celltypes %>% map(function(j) {
#   file <- file.path(args$differential_results_folder,sprintf("%s_vs_%s.txt.gz",i,j))
#   if (file.exists(file)) {
#     fread(file, select = c(1,2,4)) %>% 
#       .[abs(logFC)>=args$min_fold_change & padj_fdr<=args$fdr] %>%
#       .[,c("celltypeA","celltypeB"):=list(i,j)] %>%
#       return
#   } }) %>% rbindlist }) %>% rbindlist %>%
#   .[,celltypeA:=factor(celltypeA,levels=opts$celltypes)] %>%
#   .[,celltypeB:=factor(celltypeB,levels=opts$celltypes)] %>%
#   .[,direction:=as.factor(c("down","up"))[as.numeric(logFC<0)+1]]

diff_results.dt <- fread(args$differential_results) %>% 
  .[abs(logFC)>=args$min_fold_change & padj_fdr<=args$fdr] %>%
  .[,celltypeA:=factor(celltypeA,levels=opts$celltypes)] %>%
  .[,celltypeB:=factor(celltypeB,levels=opts$celltypes)] %>%
  .[,direction:=as.factor(c("down","up"))[as.numeric(logFC<0)+1]]

ncelltypes <- unique(c(as.character(unique(diff_results.dt$celltypeA)),as.character(unique(diff_results.dt$celltypeB)))) %>% length

# diff_stats.dt <- unique(diff_results.dt[,c("celltypeA","celltypeB")])

#########################
## Define marker genes ##
#########################

foo <- diff_results.dt[,.(score=sum(direction=="up")), by=c("celltypeA","gene")] %>% setnames("celltypeA","celltype")
bar <- diff_results.dt[,.(score=sum(direction=="down")), by=c("celltypeB","gene")] %>% setnames("celltypeB","celltype")

markers_genes.dt <- merge(foo,bar,by=c("celltype","gene"), all=TRUE) %>% 
  .[is.na(score.x),score.x:=0] %>% .[is.na(score.y),score.y:=0] %>%
  .[,score:=score.x+score.y] %>%
  .[,c("score.x","score.y"):=NULL] %>%
  .[,score:=round(score/(ncelltypes-1),2)] %>%
  setorder(celltype,-score)
rm(foo,bar)

stopifnot(max(markers_genes.dt$score,na.rm=T)==1)

##########
## Save ##
##########

# Save marker score for all combination of genes and cell types
length(unique(markers_genes.dt$gene))
length(unique(markers_genes.dt$celltype))
fwrite(markers_genes.dt, file.path(args$outdir,"marker_genes_all.txt.gz"), sep="\t")

# Save marker score for strong markers
markers_genes_filt.dt <- markers_genes.dt %>% .[score>=args$min_score]
length(unique(markers_genes_filt.dt$gene))
length(unique(markers_genes_filt.dt$celltype))
fwrite(markers_genes_filt.dt, file.path(args$outdir,"marker_genes_filtered.txt.gz"), sep="\t")

