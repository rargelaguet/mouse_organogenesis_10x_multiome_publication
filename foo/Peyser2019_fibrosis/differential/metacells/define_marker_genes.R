source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# Options
opts$min_score <- 0.75
opts$min_fold_change <- 1.0
opts$fdr <- 0.01
opts$group_variable <- "cluster"
opts$groups <- opts$clusters

# I/O
io$sce.pseudobulk <- file.path(io$basedir,sprintf("results/rna/pseudobulk/%s/SingleCellExperiment_pseudobulk.rds",opts$group_variable))
io$differential_results <- file.path(io$basedir,sprintf("results/rna/differential/metacells/%s/parsed/diff_expr_results.txt.gz",opts$group_variable))
io$outdir <- file.path(io$basedir,sprintf("results/rna/differential/metacells/%s/parsed",opts$group_variable)); dir.create(io$outdir, showWarnings = F)
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
  .[,score:=round(score/(ngroups-1),2)]

stopifnot(max(markers_genes.dt$score,na.rm=T)==1)

#######################
## Add average logFC ##
#######################

tmp <- diff_results.dt %>%
  .[,.(logFC=round(mean(logFC),2)), by=c("groupB","gene")] %>%
  setnames("groupB","group")

markers_genes.dt <- markers_genes.dt %>% 
  merge(tmp, by=c("group","gene"))

######################################
## Add mean expression (pseudobulk) ##
######################################

pseudobulk.sce <- readRDS(io$sce.pseudobulk)[,opts$groups]

genes.to.use <- intersect(unique(markers_genes.dt$gene),rownames(pseudobulk.sce))
pseudobulk.sce <- pseudobulk.sce[genes.to.use,]

expr.dt <- assay(pseudobulk.sce,"logcounts") %>% round(2) %>%
  as.data.table(keep.rownames=T) %>%
  setnames("rn","gene") %>%
  melt(id.vars="gene", variable.name="group", value.name="expr_pseudobulk")

markers_genes.dt <- markers_genes.dt %>% 
  merge(expr.dt, by=c("group","gene")) %>%
  setorder(group,-score)


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
table(markers_genes_filt.dt$group)
fwrite(markers_genes_filt.dt, file.path(io$outdir,"marker_genes_filtered.txt.gz"), sep="\t")

