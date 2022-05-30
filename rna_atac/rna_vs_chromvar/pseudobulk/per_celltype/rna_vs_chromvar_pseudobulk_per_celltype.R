here::i_am("rna_atac/rna_vs_chromvar/pseudobulk/per_celltype/rna_vs_chromvar_pseudobulk_per_celltype.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--sce',                type="character",              help='RNA SingleCellExperiment (pseudobulk)') 
p$add_argument('--motif_annotation',  type="character",              help='Motif annotation') 
p$add_argument('--atac_chromvar_matrix',  type="character",              help='ATAC chromVAR matrix (pseudobulk)') 
p$add_argument('--motif2gene',   type="character",              help='') 
p$add_argument('--cor_rna_vs_chromvar_per_gene',  type="character",              help='Motif annotation') 
p$add_argument('--outdir',          type="character",                help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
# io$basedir <- file.path(io$basedir,"test")
# args <- list()
# args$motif_annotation <- "CISBP"
# args$sce <- file.path(io$basedir,"results/rna/pseudobulk/celltype/SingleCellExperiment_pseudobulk.rds")
# args$atac_chromvar_matrix <- file.path(io$basedir,sprintf("results/atac/archR/chromvar/pseudobulk/celltype/chromVAR_deviations_%s_pseudobulk_archr.rds",args$motif_annotation))
# args$motif2gene <- file.path(io$basedir,sprintf("processed/atac/archR/Annotations/%s_motif2gene.txt.gz",args$motif_annotation))
# args$cor_rna_vs_chromvar_per_gene <- file.path(io$basedir,sprintf("results/rna_atac/rna_vs_chromvar/pseudobulk/per_gene/cor_rna_vs_chromvar_%s_per_gene_pseudobulk.txt.gz",args$motif_annotation))
# args$outdir <- file.path(io$basedir,sprintf("results/rna_atac/rna_vs_chromvar/pseudobulk/per_celltype/%s",args$motif_annotation))
## END TEST ##

# I/O
dir.create(args$outdir, showWarnings=F)

#######################################
## Load pseudobulk RNA and ATAC data ##
#######################################

# Load SingleCellExperiment
rna_pseudobulk.sce <- readRDS(args$sce)

# Load ATAC SummarizedExperiment
atac_chromvar_pseudobulk.se <- readRDS(args$atac_chromvar_matrix)
# rownames(atac_chromvar_pseudobulk.se)[rownames(atac_chromvar_pseudobulk.se)=="MA0009.1"] <- "TBXT"

################
## Parse data ##
################

# Load motif2gene annotation
motif2gene.dt <- fread(args$motif2gene) %>%
  .[motif%in%rownames(atac_chromvar_pseudobulk.se) & gene%in%toupper(rownames(rna_pseudobulk.sce))]# %>%
  # .[,N:=length(unique(motif)),by="gene"] %>% .[N==1] %>% .[,N:=NULL]

# Duplicated motifs
tmp <- motif2gene.dt[,.N,"gene"][N>1]
if (nrow(tmp)>0) {
  print(sprintf("%s genes with duplicated motifs:",nrow(tmp)))
  print(motif2gene.dt[gene%in%tmp$gene] %>% setorder(-gene))
}

# Subset TFs
rna_tf_pseudobulk.sce <- rna_pseudobulk.sce[rownames(rna_pseudobulk.sce)%in%str_to_title(motif2gene.dt$gene),]
rownames(rna_tf_pseudobulk.sce) <- toupper(rownames(rna_tf_pseudobulk.sce))
atac_chromvar_pseudobulk.se <- atac_chromvar_pseudobulk.se[rownames(atac_chromvar_pseudobulk.se)%in%motif2gene.dt$motif,]

###########################################
## Convert to long data.tables and merge ##
###########################################

rna_tf_pseudobulk.dt <- logcounts(rna_tf_pseudobulk.sce) %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","gene") %>%
  melt(id.vars="gene", variable.name="celltype", value.name="expr")

atac_chromvar_pseudobulk.dt <- as.matrix(assay(atac_chromvar_pseudobulk.se,"z")) %>% t %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","celltype") %>%
  melt(id.vars=c("celltype"), variable.name="motif", value.name="chromvar_zscore") %>%
  merge(motif2gene.dt[,c("motif","gene")], by="motif")

# Merge
rna_atac.dt <- merge(rna_tf_pseudobulk.dt, atac_chromvar_pseudobulk.dt, by = c("gene","celltype"))

# Print stats
print(sprintf("Number of TFs: %s",length(unique(rna_atac.dt$gene))))
print(sprintf("Number of celltypes: %s",length(unique(rna_atac.dt$celltype))))

#############################################
## Load pre-computed correlation estimates ##
#############################################

cor.dt <- fread(args$cor_rna_vs_chromvar_per_gene) %>%
  .[,cor_sign:=as.factor(c("Repressor","Activator")[(r>0)+1])]

###############################
## Scatterplot per cell type ##
###############################

opts$max.chromvar <- 12
opts$max.expr <- 12

celltypes.to.plot <- unique(rna_atac.dt$celltype)

# i <- "Neural_crest"
for (i in celltypes.to.plot) {
  
  to.plot <- rna_atac.dt[celltype==i]  %>%
    # merge(cor.dt[,c("gene","cor_sign")], by="gene") %>%
    .[chromvar_zscore<0,chromvar_zscore:=0] %>%
    .[chromvar_zscore>=opts$max.chromvar,chromvar_zscore:=opts$max.chromvar] %>%
    .[expr>=opts$max.expr,expr:=opts$max.expr]
  
  to.plot.text <- to.plot[expr>5 & chromvar_zscore>3.5] 
  to.plot.dots <- to.plot[!gene%in%to.plot.text$gene] 
  
  to.plot %>% .[,dot_size:=minmax.normalisation(expr)*minmax.normalisation(chromvar_zscore)]
  
  p <- ggplot(to.plot, aes(x=chromvar_zscore, y=expr)) +
    geom_point(aes(size=dot_size), shape=21) + 
    # geom_point(size=0.5, data=to.plot.text) +
    ggrepel::geom_text_repel(data=to.plot.text, aes(label=gene), size=4) +
    # geom_text(aes(label=gene), size=3, data=to.plot.text) +
    scale_size_continuous(range = c(0.1,7)) +
    # scale_fill_gradient(low = "gray80", high = "darkgreen") +
    coord_cartesian(
      xlim = c(0,opts$max.chromvar+0.1),
      ylim = c(min(rna_atac.dt$expr),opts$max.expr+0.1)
      ) +
    labs(x="Motif accessibility (chromVAR, z-score)", y="Gene expression") +
    guides(size="none", fill = guide_legend(override.aes = list(size=3))) +
    scale_fill_brewer(palette="Dark2") +
    theme_classic() +
    theme(
      legend.position = "none",
      legend.title = element_blank(),
      axis.text = element_text(size=rel(0.75), color="black")
    )
  
  outfile <- file.path(args$outdir,sprintf("%s_rna_vs_chromvar_%s_pseudobulk.pdf",i,args$motif_annotation))
  pdf(outfile, width = 6, height = 5)
  print(p)
  dev.off()
} 

# Create a completion token
# file.create(file.path(args$outdir,sprintf("%s_completed.txt",args$motif_annotation)))
file.create(file.path(args$outdir,"completed.txt"))


#############
## Explore ##
#############

# # Quantify number of active TFs per cell type

# to.plot <- rna_atac.dt %>%
#   .[,sum(expr>5 & chromvar_zscore>15),by="celltype"] %>%
#   # .[,sum(chromvar_zscore>25),by="celltype"] %>%
#   setorder(-V1) %>% .[,cellype:=factor(celltype,levels=celltype)]


# p <- ggbarplot(to.plot, x="celltype", y="V1", fill="celltype", sort.val = "asc") +
#   scale_fill_manual(values=opts$celltype.colors) +
#   labs(x="", y="Number of active TFs") +
#   guides(x = guide_axis(angle = 90)) +
#   theme_classic() +
#   theme(
#     legend.position = "none",
#     axis.text.x = element_text(color="black", size=rel(0.8)),
#     axis.title.x = element_blank(),
#     axis.ticks.x = element_blank()
# )

# pdf(sprintf("%s/number_active_TFs.pdf",io$outdir), width = 9, height = 5)
# print(p)
# dev.off()


# to.plot <- rna_atac.dt %>%
#   merge(cor_rna_vs_chromvar_per_gene.dt[,c("gene","cor_sign")], by="gene") %>%
#   .[,sum(expr>5 & chromvar_zscore>15),by=c("celltype","cor_sign")]# %>%
#   # .[,celltype:=factor(celltype,levels=opts$celltypes)]# %>%
#   # dcast(celltype~cor_sign, value.var="V1") %>%
#   # .[,ratio_activators:=Activator/(Activator+Repressor)]

# celltype.order <- to.plot %>%
#   dcast(celltype~cor_sign, value.var="V1") %>%
#   .[,ratio_activators:=Activator/(Activator+Repressor)] %>%
#   setorder(ratio_activators) %>% .$celltype
# to.plot[,celltype:=factor(celltype,levels=celltype.order)]

# p <- ggbarplot(to.plot, x="celltype", y="V1", fill="cor_sign") +
#   scale_fill_brewer(palette="Dark2") +
#   # scale_fill_manual(values=opts$celltype.colors) +
#   labs(x="", y="Number of active TFs") +
#   guides(x = guide_axis(angle = 90)) +
#   theme_classic() +
#   theme(
#     legend.position = "top",
#     axis.text.x = element_text(color="black", size=rel(0.8)),
#     axis.title.x = element_blank(),
#     axis.ticks.x = element_blank()
#   )
