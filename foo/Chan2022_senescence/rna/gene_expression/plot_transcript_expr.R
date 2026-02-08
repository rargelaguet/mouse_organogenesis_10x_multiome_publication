source(here::here("settings.R"))
source(here::here("utils.R"))

library(pheatmap)

#####################
## Define settings ##
#####################

# I/O
io$metadata <- file.path(io$basedir,"processed/rna/sample_metadata.txt")
io$sce <- file.path(io$basedir,"processed/rna/salmon/salmon_transcripts_SingleCellExperiment.rds")
io$outdir <- file.path(io$basedir,"results/rna/transcript_expression"); dir.create(io$outdir, showWarnings = F, recursive = T)

# Options
# opts$plot_logcounts <- TRUE

###################
## Load metadata ##
###################

sample_metadata.dt <- fread(io$metadata) %>%
  .[,sample:=NULL] %>% setnames("alias","sample") %>%
  .[,passage_string:=sprintf("Passage %s",Passage)] %>%
  .[,passage_string:=factor(passage_string,levels=sprintf("Passage %s",opts$passages))]

#########################
## Load RNA expression ##
#########################

sce <- readRDS(io$sce)[,sample_metadata.dt$sample]

colData(sce) <- sample_metadata.dt %>% tibble::column_to_rownames("sample") %>% DataFrame

########################
## Load gene metadata ##
########################

transcript_metadata.dt <- rowData(sce) %>% as.data.table(keep.rownames=T) %>% setnames("rn","transcript_id")

###################
## Normalisation ##
###################

# if (opts$plot_logcounts) {
#   assay(sce,"logcounts") <- log2(1e6*(sweep(assay(sce,"counts"),2,colSums(counts(sce)),"/"))+1)
# } else { 
#   assay(sce,"logcounts") <- 1e6*(sweep(assay(sce,"counts"),2,colSums(counts(sce)),"/"))
# }

################
## Line plots ##
###############

genes.to.plot <- transcript_metadata.dt[symbol=="JARID2",transcript_id]

genes.to.plot <- genes.to.plot[genes.to.plot%in%rownames(sce)]

# i <- "ENST00000582401"
for (i in genes.to.plot) {
  
  to.plot <- data.table(
    sample = colnames(sce),
    expr = assay(sce,"tpm")[i,]
  ) %>% merge(sample_metadata.dt,by="sample")
  
  gene_name_i <- transcript_metadata.dt[transcript_id==i,ens_id]
  ens_id_i <- transcript_metadata.dt[transcript_id==i,symbol]
  
  # to.plot[,Passage:=factor(Passage,levels=opts$passages)]
  to.plot.means <- to.plot[,.(expr=mean(expr), sd=sd(expr)), by=c("Passage","Treatment")]
  
  ylim.max <- max(5,max(to.plot$expr))
  
  p <- ggplot(to.plot, aes(x=Passage, y=expr)) +
    geom_line(aes(color=Treatment), alpha=0.75, data=to.plot.means) +
    geom_point(aes(fill=Treatment), size=4, shape=21, stroke=0.2, color="black") +
    geom_errorbar(aes(ymin=expr-sd, ymax=expr+sd), color="black", width=0.4, alpha=0.9, linewidth=0.5, data=to.plot.means) +
    # scale_fill_manual(values=c(brewer.pal(n=9,"OrRd"),"gray70")) +
    # facet_wrap(~Treatment, nrow=1, scales="fixed") +
    # coord_cartesian(ylim=c(0,ylim.max)) +
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
  
  pdf(file.path(io$outdir,sprintf("%s_%s_%s_expr_lineplots_Chan2022_senescence.pdf",i,ens_id_i,gene_name_i)), width=6, height=5)
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


###############
## Box plots ##
###############

genes.to.plot <- genes.to.plot[genes.to.plot%in%rownames(sce) ]

samples.to.plot <- c("p8_Control","p8_TG","p11_Control","p11_TG","p24_Control","p24_TG","p29_Control","p29_TG")

# i <- "ST6GAL1"
for (i in genes.to.plot) {
  
  to.plot <- data.table(
    sample = colnames(sce),
    expr = assay(sce,"logcounts")[i,]
  ) %>% merge(sample_metadata.dt,by="sample") %>% .[Group%in%samples.to.plot]
  
  to.plot.means <- to.plot[,.(expr=mean(expr), sd=sd(expr)), by=c("Passage","Treatment")]
  
  p <- ggplot(to.plot, aes(x=Treatment, y=expr, fill=Treatment)) +
    # geom_bar(stat="identity", color="black", data=to.plot.means) +
    # geom_errorbar(aes(ymin=signature-sd, ymax=signature+sd), color="black", width=0.4, alpha=0.9, linewidth=0.5, data=to.plot.means) +
    geom_boxplot() +
    geom_point(size=3, shape=21, stroke=0.2, color="black") +
    facet_wrap(~passage_string, nrow=1) +
    labs(x="", y=sprintf("RNA expression (%s)",ifelse(opts$plot_logcounts, "log cpm", "cpm")), title=i) +
    theme_classic() +
    theme(
      strip.text = element_text(size=rel(1)),
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(colour="black",size=rel(1.20)),
      axis.text.y = element_text(colour="black",size=rel(1.0)),
      axis.title.y = element_text(colour="black",size=rel(1.15)),
      legend.position = "none",
      legend.title = element_blank(),
    )
  
  pdf(file.path(io$outdir,sprintf("%s_expr_boxplots_Chan2022_senescence.pdf",i)), width=7, height=5)
  print(p)
  dev.off()
}

################################################
## Heatmap of control samples: passage effect ##
################################################

genes.to.plot <- grep("HOX",genes.to.plot,value=T)

samples.order <- sample_metadata.dt[disease=="normal"] %>% setorder(mean_age) %>% .$sample

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
  filename = file.path(io$outdir,"heatmaps/gene_expr_heatmap_Hox_Chan2022_senescence.pdf")
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
  # filename = file.path(io$outdir,"heatmaps/TE_expr_heatmap_scaled_Hox_Chan2022_senescence.pdf")
)
