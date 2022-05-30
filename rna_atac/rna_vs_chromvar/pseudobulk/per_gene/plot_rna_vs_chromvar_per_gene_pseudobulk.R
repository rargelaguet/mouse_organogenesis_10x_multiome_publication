here::i_am("rna_atac/rna_vs_chromvar/pseudobulk/per_gene/plot_rna_vs_chromvar_per_gene_pseudobulk.R")

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
p$add_argument('--cor_results',  type="character",              help='') 
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
# args$cor_results <- file.path(io$basedir,sprintf("results/rna_atac/rna_vs_chromvar/pseudobulk/per_gene/%s/cor_rna_vs_chromvar_%s_per_gene_pseudobulk.txt.gz",args$motif_annotation,args$motif_annotation))
# args$outdir <- file.path(io$basedir,sprintf("results/rna_atac/rna_vs_chromvar/pseudobulk/per_gene/%s/pdf",args$motif_annotation))
## END TEST ##

# I/O
dir.create(args$outdir, showWarnings=F, recursive = T)

# Options
# opts$motif_annotation <- args$motif_annotation
opts$threshold_fdr <- 0.10

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

cor.dt <- fread(args$cor_results) %>%
  .[,cor_sign:=as.factor(c("Repressor","Activator")[(r>0)+1])] %>% 
  .[,sig:=p<=opts$threshold_fdr]

##################
## Volcano plot ##
##################

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
  labs(x="Pearson correlation (TF expr vs motif acc. (chromVAR z-score)", y=expression(paste("-log"[10],"(p.value)"))) +
  theme_classic() +
  theme(
    axis.text = element_text(size=rel(0.75), color='black'),
    axis.title = element_text(size=rel(1.0), color='black'),
    legend.position="none"
  )

pdf(file.path(args$outdir,sprintf("volcano_rna_vs_chromvar_%s_per_gene_pseudobulk.pdf",args$motif_annotation)), width = 7, height = 5)
print(p)
dev.off()


######################################
## Scatter plot of individual genes ##
######################################

facet.labels <- c(expr = "RNA expression", chromvar_zscore = "Motif accessibility (z-score)")

genes.to.plot <- unique(rna_atac.dt$gene)

# i <- "T"; j <- "T_789"
for (i in genes.to.plot) {
  
  motifs.to.plot <- unique(rna_atac.dt[gene==i,motif])
  for (j in motifs.to.plot) {
    to.plot <- rna_atac.dt[gene==i & motif==j]
    to.plot2 <- to.plot %>% melt(id.vars=c("celltype","gene","motif"))
    
    p1 <- ggscatter(to.plot, x="expr", y="chromvar_zscore", fill="celltype", size=4, shape=21, 
                    add="reg.line", add.params = list(color="black", fill="lightgray"), conf.int=TRUE) +
      stat_cor(method = "pearson", label.x.npc = "middle", label.y.npc = "bottom") +
      scale_fill_manual(values=opts$celltype.colors) +
      labs(x=sprintf("%s expression",i), y=sprintf("%s Motif acc. (chromVAR z-score)",i)) +
      guides(fill="none") +
      theme(
        axis.text = element_text(size=rel(0.7))
      )
    
    p2 <- ggbarplot(to.plot2, x="celltype", y="value", fill="celltype") +
      facet_wrap(~variable, nrow=2, scales="free_y",  labeller = as_labeller(facet.labels)) +
      scale_fill_manual(values=opts$celltype.colors) +
      geom_hline(yintercept=0, linetype="dashed") +
      labs(x="", y="") +
      theme_classic() +
      guides(x = guide_axis(angle = 90)) +
      theme(
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()
      )
    
    p <- cowplot::plot_grid(plotlist=list(p1,p2), nrow = 1, rel_widths = c(1/2,1/2))
    
    # pdf(sprintf("%s/individual_genes/%s_rna_vs_chromvar.pdf",args$outdir,i), width=10, height=4)
    png(file.path(args$outdir,sprintf("%s(%s_%s)_rna_vs_chromvar_pseudobulk.png",i,j,args$motif_annotation)), width = 1000, height = 500)
    print(p)
    dev.off()
  }
}

# Completion token
file.create(file.path(args$outdir,"completed.txt"))