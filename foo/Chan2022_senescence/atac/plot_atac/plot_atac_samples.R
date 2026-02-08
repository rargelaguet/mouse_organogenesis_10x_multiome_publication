source(here::here("settings.R"))
source(here::here("utils.R"))

library(pheatmap)

#####################
## Define settings ##
#####################

# I/O
io$metadata <- file.path(io$basedir,"processed/atac/sample_metadata.txt")
io$sce <- paste0(io$basedir,"/processed/atac/SingleCellExperiment.rds")
io$outdir <- file.path(io$basedir,"results/atac"); dir.create(io$outdir, showWarnings = F, recursive = T)

# Options
opts$plot_logcounts <- FALSE

###################
## Load metadata ##
###################

sample_metadata.dt <- fread(io$metadata)

###############
## Load data ##
###############

sce <- readRDS(io$sce)[,sample_metadata.dt$sample]
# colnames(sce) <- sample_metadata.dt$sample
# rownames(sce) <- rowData(sce)$symbol

colData(sce) <- sample_metadata.dt %>% tibble::column_to_rownames("sample") %>% DataFrame

###########################
## Load feature metadata ##
###########################

# peak_metadata.dt <- fread(io$gene_metadata)

io$limma_results <- file.path(io$basedir,"results/atac/limma/control_samples/limma_results.rds")
limma.list <- readRDS(io$limma_results)
diff_acc.dt <- limma.list[["Passage"]] %>% as.data.table 

###################
## Normalisation ##
###################

if (opts$plot_logcounts) {
  assay(sce,"logcounts") <- log2(1e6*(sweep(assay(sce,"counts"),2,colSums(counts(sce)),"/"))+1) # log cpm
} else { 
  assay(sce,"logcounts") <- 1e6*(sweep(assay(sce,"counts"),2,colSums(counts(sce)),"/")) # cpm
}

################
## Line plots ##
###############

# i <- "Peak_38866"
for (i in genes.to.plot) {
  
  to.plot <- data.table(
    sample = colnames(sce),
    expr = assay(sce,"logcounts")[i,]
  ) %>% merge(sample_metadata.dt,by="sample")
  
  # to.plot[,Passage:=factor(Passage,levels=opts$passages)]
  to.plot.means <- to.plot[,.(expr=mean(expr), sd=sd(expr)), by=c("Passage","Treatment")]
  
  ylim.max <- max(5,max(to.plot$expr))
  
  p <- ggplot(to.plot, aes(x=Passage, y=expr)) +
    geom_line(aes(color=Treatment), alpha=0.75, data=to.plot.means) +
    geom_point(aes(fill=Treatment), size=4, shape=21, stroke=0.2, color="black") +
    geom_errorbar(aes(ymin=expr-sd, ymax=expr+sd), color="black", width=0.4, alpha=0.9, linewidth=0.5, data=to.plot.means) +
    # scale_fill_manual(values=c(brewer.pal(n=9,"OrRd"),"gray70")) +
    # facet_wrap(~Treatment, nrow=1, scales="fixed") +
    # facet_wrap(~sex) +
    coord_cartesian(ylim=c(0,ylim.max)) +
    labs(x="Passage number", y=sprintf("RNA expression (%s)",ifelse(opts$plot_logcounts, "log cpm", "cpm")), title=i) +
    theme_classic() +
    theme(
      strip.text = element_text(size=rel(1.15)),
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(colour="black",size=rel(1.10)),
      axis.text.y = element_text(colour="black",size=rel(1.10)),
      legend.position = "top",
      legend.title = element_blank(),
    )
  
  pdf(file.path(io$outdir,sprintf("%s_expr_lineplots_Chan2022_senescence.pdf",i)), width=6, height=5)
  print(p)
  dev.off()
}


###############
## Bar plots ##
###############

genes.to.plot <- gene_metadata.dt[chr=="chrMT",symbol]
genes.to.plot <- TEs
genes.to.plot <- gene_metadata_champ.dt[abs(tg_logFC)>=2 & tg_padj<0.01,symbol]
genes.to.plot <- gene_metadata_champ.dt[abs(passage_logFC)>=1 & passage_padj<0.01,symbol]

genes.to.plot <- genes.to.plot[genes.to.plot%in%rownames(sce) ]

# i <- "ST6GAL1"
for (i in genes.to.plot) {
  
  to.plot <- data.table(
    sample = colnames(sce),
    expr = assay(sce,"logcounts")[i,]
  ) %>% merge(sample_metadata.dt,by="sample")
  
  to.plot[,Passage:=factor(Passage,levels=opts$passages)]
  to.plot.means <- to.plot[,.(expr=mean(expr), sd=sd(expr)), by=c("Passage","Treatment")]
  
  ylim.max <- max(5,max(to.plot$expr))
  p <- ggplot(to.plot, aes(x=Passage, y=expr)) +
    geom_bar(stat="identity", alpha=0.75, fill="gray75", color="black", data=to.plot.means) +
    geom_point(size=4, shape=21, stroke=0.2, color="black", fill="gray70") +
    geom_errorbar(aes(ymin=expr-sd, ymax=expr+sd), color="black", width=0.4, alpha=0.9, linewidth=0.5, data=to.plot.means) +
    # scale_fill_manual(values=c(brewer.pal(n=9,"OrRd"),"gray70")) +
    facet_wrap(~Treatment, nrow=1, scales="fixed") +
    # facet_wrap(~sex) +
    coord_cartesian(ylim=c(0,ylim.max)) +
    labs(x="Passage number", y=sprintf("RNA expression (%s)",ifelse(opts$plot_logcounts, "log cpm", "cpm")), title=i) +
    theme_classic() +
    theme(
      strip.text = element_text(size=rel(1.15)),
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(colour="black",size=rel(1.10)),
      axis.text.y = element_text(colour="black",size=rel(1.10)),
      legend.position = "none"
    )
  
  pdf(file.path(io$outdir,sprintf("%s_expr_barplots_Chan2022_senescence.pdf",i)), width=8, height=5)
  print(p)
  dev.off()
}

#############
## Heatmap ##
#############

samples.order <- sample_metadata.dt[disease=="normal"] %>% setorder(mean_age) %>% .$sample

# tmp <- apply(assay(sce[,samples.order],"logcounts"),1,max) %>% sort
# features.max <- tmp[tmp>=7.5] %>% names
# tmp <- apply(assay(sce[,samples.order],"logcounts"),1,var) %>% sort
# features.var <- tmp[tmp>=0.10] %>% names
# genes.to.plot <- intersect(features.var,features.max)
# genes.to.plot <- fread("/Users/rargelaguet/data/DellaValle2022_line1_progeria/results/differential/genes/HGPS_diff.txt.gz")[abs(logFC)>=2 & padj_fdr<0.01,gene]
genes.to.plot <- fread("/Users/rargelaguet/data/Tsitsipatis2023_human_fibroblast_aging/results/differential/genes/pc1_group/diff.txt.gz")[abs(logFC)>=1 & padj_fdr<0.01,gene]
# genes.to.plot <- c("EZH1","EZH2","SUZ12","EED","RING1")
# genes.to.plot <- grep("HOX",genes.to.plot,value=T)

to.plot.mtx <- assay(sce,"logcounts")[genes.to.plot,samples.order] %>% as.matrix 
# to.plot.mtx[to.plot.mtx>=16] <- 16

breaks <- seq(from=min(abs(to.plot.mtx)), to=max(abs(to.plot.mtx)), by=1)

pheatmap(
  mat = to.plot.mtx,
  breaks = breaks,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaks)),
  cluster_cols = F, cluster_rows = T,
  show_colnames = F, show_rownames = T,
  legend = T, annotation_legend = F,
  fontsize_row = 4, fontsize_col = 8,
  treeheight_row = 10, treeheight_col = 10,
  scale = "none",
  filename = file.path(io$outdir,"heatmaps/gene_expr_heatmap_Hox_Tsitsipatis2023.pdf")
)


# scaled
to.plot.mtx <- assay(sce,"logcounts")[genes.to.plot,samples.order] %>% as.matrix 

pheatmap(
  mat = to.plot.mtx,
  # breaks = breaks,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaks)),
  cluster_cols = F, cluster_rows = T,
  show_colnames = F, show_rownames = T,
  legend = T, annotation_legend = F,
  fontsize_row = 4, fontsize_col = 8,
  treeheight_row = 10, treeheight_col = 10,
  scale = "row"
  # filename = file.path(io$outdir,"heatmaps/TE_expr_heatmap_scaled_Hox_Tsitsipatis2023.pdf")
)

##########
## TEST ##
##########

to.plot <- data.table(
  sample = colnames(sce),
  ngenes = colSums(counts(sce)>=10)
) %>% merge(sample_metadata.dt,by="sample")

ggboxplot(to.plot, x="age_group", y="ngenes")
