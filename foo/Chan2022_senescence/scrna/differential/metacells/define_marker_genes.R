# here::i_am("rna/differential/metacells/celltype/analysis/define_marker_genes.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

io$differential_results <- file.path(io$basedir,"results/rna/differential/metacells/alias/parsed/diff_expr_results.txt.gz")
io$differential_stats <- file.path(io$basedir,"results/rna/differential/metacells/alias/parsed/diff_expr_stats.txt.gz")
opts$min_score <- 0.75
opts$min_fold_change <- 1.0
opts$fdr <- 0.01
io$outdir <- file.path(io$basedir,"results/rna/differential/metacells/alias/parsed"); dir.create(io$outdir, showWarnings = F)

# opts$clusters <- paste0("c",1:14)
# opts$groups <- opts$clusters
opts$groups <- opts$aliases

dir.create(io$outdir, showWarnings = F)

###############################
## Load differential results ##
###############################

diff_results.dt <- fread(io$differential_results) %>% 
  .[abs(logFC)>=opts$min_fold_change & padj_fdr<=opts$fdr] %>%
  .[,groupA:=factor(groupA,levels=opts$groups)] %>%
  .[,groupB:=factor(groupB,levels=opts$groups)] %>%
  .[,direction:=as.factor(c("down","up"))[as.numeric(logFC<0)+1]]

ngroups <- unique(c(as.character(unique(diff_results.dt$groupA)),as.character(unique(diff_results.dt$groupB)))) %>% length

#########################
## Define marker genes ##
#########################

markers_genes.dt <- diff_results.dt %>% 
  .[,.(score=sum(direction=="up")), by=c("groupA","gene")] %>% setnames("groupA","group") %>%
  .[,score:=round(score/(ngroups-1),2)] %>%
  setorder(group,-score)
# bar <- diff_results.dt[,.(score=sum(direction=="down")), by=c("groupB","gene")] %>% setnames("groupB","group")

stopifnot(max(markers_genes.dt$score,na.rm=T)==1)

##########
## Save ##
##########

# Save marker score for all combination of genes and cell types
length(unique(markers_genes.dt$gene))
length(unique(markers_genes.dt$group))
fwrite(markers_genes.dt, file.path(io$outdir,"marker_genes_all.txt.gz"), sep="\t")

# Save marker score for strong markers
markers_genes_filt.dt <- markers_genes.dt %>% .[score>=opts$min_score]
length(unique(markers_genes_filt.dt$gene))
length(unique(markers_genes_filt.dt$group))
fwrite(markers_genes_filt.dt, file.path(io$outdir,"marker_genes_filtered.txt.gz"), sep="\t")

