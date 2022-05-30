here::i_am("rna_atac/rna_vs_acc/pseudobulk/gene_expr_vs_promoter_acc/plot_gene_expr_vs_promoter_acc_pseudobulk.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

library(patchwork)

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--sce',  type="character",              help='RNA SingleCellExperiment (pseudobulk)') 
p$add_argument('--gene_score_matrix',  type="character",              help='ATAC Gene score matrix (pseudobulk)') 
p$add_argument('--cor_results',  type="character",              help='RNA expression vs promoter accessibility correlation results') 
p$add_argument('--outdir',          type="character",                help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
args <- list()
args$sce <- file.path(io$basedir,"results/rna/pseudobulk/SingleCellExperiment_pseudobulk_celltype.mapped.rds")
args$gene_score_matrix <- file.path(io$basedir,"results/atac/archR/pseudobulk/celltype.mapped/pseudobulk_GeneScoreMatrix_TSS_summarized_experiment.rds") # io$archR.pseudobulk.GeneMatrix.se
args$cor_results <- file.path(io$basedir,"results/rna_atac/gene_expr_vs_promoter_acc/pseudobulk/cor_gene_expr_vs_promoter_acc_pseudobulk.txt.gz")
args$outdir <- file.path(io$basedir,"results/rna_atac/gene_expr_vs_promoter_acc/pseudobulk")
## END TEST ##

#####################
## Define settings ##
#####################

# I/O
dir.create(file.path(args$outdir,"individual_genes"), showWarnings=FALSE, recursive=TRUE)

# Options
opts$threshold_fdr <- 0.10

#######################################
## Load pseudobulk RNA and ATAC data ##
#######################################

# Load SingleCellExperiment
rna_pseudobulk.sce <- readRDS(args$sce)

# Load ATAC SummarizedExperiment
atac_pseudobulk_GeneScoreMatrix.se <- readRDS(args$gene_score_matrix)

# Rename genes
# if (any(grepl("f+",rownames(atac_pseudobulk_GeneScoreMatrix.se)))) {
#   rownames(atac_pseudobulk_GeneScoreMatrix.se) <- rowData(atac_pseudobulk_GeneScoreMatrix.se)$name
# }

###########################################
## Convert to long data.tables and merge ##
###########################################

rna_pseudobulk.dt <- logcounts(rna_pseudobulk.sce) %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","gene") %>%
  melt(id.vars="gene", variable.name="celltype", value.name="expr")

atac_gene_scores_pseudobulk.dt <- as.matrix(assay(atac_pseudobulk_GeneScoreMatrix.se)) %>% t %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","celltype") %>%
  melt(id.vars=c("celltype"), variable.name="gene", value.name="acc")

# Merge
rna_atac.dt <- merge(rna_pseudobulk.dt, atac_gene_scores_pseudobulk.dt, by = c("gene","celltype"))

#############################################
## Load pre-computed correlation estimates ##
#############################################

cor.dt <- fread(args$cor_results) %>%
  .[,cor_sign:=as.factor(c("Repressor","Activator")[(r>0)+1])] %>% 
  .[, sig:=p<=opts$threshold_fdr]

##########
## Plot ##
##########

to.plot <- cor.dt %>%
  .[padj_fdr<=1e-14,padj_fdr:=1e-14] %>%
  .[,log_pval:=-log10(padj_fdr+1e-100)] %>%
  .[,dot_size:=minmax.normalisation(log_pval)]

negative_hits <- to.plot[sig==TRUE & r<0,gene]
positive_hits <- to.plot[sig==TRUE & r>0,gene]
all <- nrow(to.plot)

xlim <- max(abs(to.plot$r), na.rm=T)
ylim <- max(-log10(to.plot$padj_fdr+1e-100), na.rm=T)

p <- ggplot(to.plot, aes(x=r, y=log_pval)) +
  geom_segment(aes(x=0, xend=0, y=0, yend=ylim-1), color="orange", size=0.5) +
  geom_jitter(aes(fill=sig, size=dot_size, alpha=dot_size), width=0.05, shape=21) + 
  ggrepel::geom_text_repel(data=head(to.plot[sig==T & r<(-0.30)],n=25), aes(x=r, y=log_pval, label=gene), size=3,  max.overlaps=100, segment.color = NA) +
  ggrepel::geom_text_repel(data=head(to.plot[sig==T & r>0.30],n=25), aes(x=r, y=log_pval, label=gene), size=3,  max.overlaps=100, segment.color = NA) +
  # ggrepel::geom_text_repel(data=to.plot[sig==T & r<(-0.50)][sample(.N,12)], aes(x=r, y=log_pval, label=gene), size=4,  max.overlaps=100, segment.color = NA) +
  # ggrepel::geom_text_repel(data=to.plot[sig==T & r>0.75][sample(.N,12)], aes(x=r, y=log_pval, label=gene), size=4,  max.overlaps=100, segment.color = NA) +
  scale_fill_manual(values=c("black","red")) +
  # scale_size_manual(values=c(0.5,1)) +
  scale_size_continuous(range = c(0.2,2)) + 
  scale_alpha_continuous(range=c(0.25,1)) +
  scale_x_continuous(limits=c(-xlim-0.15,xlim+0.15)) +
  scale_y_continuous(limits=c(0,ylim+1)) +
  annotate("text", x=0, y=ylim+1, size=4, label=sprintf("(%d)", all)) +
  annotate("text", x=-xlim-0.15, y=ylim+1, size=4, label=sprintf("%d (-)",length(negative_hits))) +
  annotate("text", x=xlim+0.15, y=ylim+1, size=4, label=sprintf("%d (+)",length(positive_hits))) +
  labs(x="Pearson correlation (RNA expression vs gene accessibility)", y=expression(paste("-log"[10],"(p.value)"))) +
  theme_classic() +
  theme(
    axis.text = element_text(size=rel(0.75), color='black'),
    axis.title = element_text(size=rel(1.0), color='black'),
    legend.position="none"
  )

pdf(file.path(args$outdir,"volcano_cor_gene_expr_vs_promoter_acc_pseudobulk.pdf"), width = 7, height = 5)
print(p)
dev.off()

###############
## Load PAGA ##
###############

opts$celltypes <- opts$celltypes[opts$celltypes%in%unique(rna_atac.dt$celltype)]

source(here::here("load_paga_graph.R"))

# Define colors
rna.col.seq <- atac.col.seq <- round(seq(0,1,0.1), 2)
rna.colors <- colorRampPalette(c("gray92", "darkgreen"))(length(rna.col.seq))
atac.colors <- colorRampPalette(c("gray92", "purple"))(length(atac.col.seq)) 

# Create base PAGA plot
p.paga <- ggnet2(
  net = net.paga,
  mode = c("x", "y"),
  node.size = 0,
  edge.size = 0.15,
  edge.color = "grey",
  label = FALSE,
  label.size = 2.3
)

##########################
## Scatterplot per gene ##
##########################

facet.labels <- c(expr = "RNA expression", acc = "Gene accessibility")

# genes.to.plot <- rna_atac.dt %>%
#   .[,.(var_expr=var(expr), var_acc=var(acc)),by="gene"] %>% 
#   setorder(-var_expr) %>% 
#   head(n=100) %>% .$gene

# genes.to.plot <- c("Ubb", "Top2a", "Actg1", "Gapdh")
genes.to.plot <- fread(io$rna.atlas.marker_genes.up)[["gene"]] %>% unique
genes.to.plot <- genes.to.plot[genes.to.plot%in%unique(rna_atac.dt$gene)]

# i <- "Ubb"
for (i in genes.to.plot) {
  outfile <- file.path(args$outdir,sprintf("individual_genes/%s_gene_expr_vs_promoter_acc_pseudobulk.png",i))
  if (!file.exists(outfile)) {
    to.plot <- rna_atac.dt[gene==i] %>% .[,celltype:=factor(celltype,levels=opts$celltypes)]
    to.plot2 <- to.plot %>% melt(id.vars=c("celltype","gene"))
    
    ## Scatterplot ##

    p.scatter <- ggscatter(to.plot, x="acc", y="expr", fill="celltype", size=4, shape=21, 
                    add="reg.line", add.params = list(color="black", fill="lightgray"), conf.int=TRUE) +
      stat_cor(method = "pearson", label.x.npc = "middle", label.y.npc = "bottom") +
      scale_fill_manual(values=opts$celltype.colors) +
      labs(y=sprintf("%s expression",i), x=sprintf("%s Gene accessibility",i)) +
      guides(fill="none") +
      theme(
        axis.text = element_text(size=rel(0.8))
      )

    ## Barplot ##
    
    p.barplot <- ggbarplot(to.plot2, x="celltype", y="value", fill="celltype") +
      facet_wrap(~variable, nrow=2, scales="free_y",  labeller = as_labeller(facet.labels)) +
      scale_fill_manual(values=opts$celltype.colors) +
      geom_hline(yintercept=0, linetype="dashed") +
      labs(x="", y="") +
      theme_classic() +
      guides(x = guide_axis(angle = 90)) +
      theme(
        axis.text.y = element_text(color="black", size=rel(0.9)),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()
      )
      

    ## PAGA ##
    expr.values <- rna_atac.dt[gene==i,c("celltype","expr")] %>% matrix.please %>% .[opts$celltypes,] %>% minmax.normalisation
    expr.colors <- round(expr.values,1) %>% map(~ rna.colors[which(rna.col.seq == .)]) %>% unlist
  
    p.paga.expr <- p.paga + geom_text(label = "\u25D0", aes(x=x, y=y), color=expr.colors, size=15, family = "Arial Unicode MS",
                  data = p.paga$data[,c("x","y")] %>% dplyr::mutate(expr=expr.colors)) +
      scale_colour_manual(values=expr.colors) + 
      labs(title="RNA expression") +
      theme(
        plot.margin = margin(t = 15, r = 15, b = 15, l = 15, unit = "pt"),
        plot.title = element_text(hjust = 0.5)
      )
  
    # motif accessibility
    acc.values <- rna_atac.dt[gene==i,c("celltype","acc")] %>% matrix.please %>% .[opts$celltypes,] %>% minmax.normalisation
    if (all(is.na(acc.values))) { acc.values <- rep(0,length(acc.values)) }
    acc.colors <- round(acc.values,1) %>% map(~ atac.colors[which(atac.col.seq == .)]) %>% unlist
    
    p.paga.acc <- p.paga + geom_text(label = "\u25D1", aes(x=x, y=y), color=acc.colors, size=15, family = "Arial Unicode MS",
                  data = p.paga$data[,c("x","y")] %>% dplyr::mutate(acc=acc.colors)) +
      scale_fill_manual(values=acc.colors) + 
      labs(title="Gene accessibility") +
      theme(
        plot.margin = margin(t = 15, r = 15, b = 15, l = 15, unit = "pt"),
        plot.title = element_text(hjust = 0.5)
      )
      
    
    p <- (p.scatter+p.barplot+p.paga.expr+p.paga.acc) + plot_layout(nrow=1, widths=c(1,1,0.6,0.6))
    
    png(outfile, width = 1600, height = 500)
    # pdf(outfile, height=8, width=11)
    print(p)
    dev.off()

  }
}

# Completion token
file.create(file.path(args$outdir,"completed.txt"))