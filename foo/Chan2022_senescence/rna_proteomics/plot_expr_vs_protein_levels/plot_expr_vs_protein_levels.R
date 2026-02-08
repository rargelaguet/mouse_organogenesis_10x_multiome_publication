source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$sample_metadata.expr <- file.path(io$basedir,"processed/rna/sample_metadata.txt")
io$sample_metadata.proteomics <- file.path(io$basedir,"processed/proteomics/sample_metadata.txt")
io$feature_metadata.proteomics <- file.path(io$basedir,"processed/proteomics/feature_metadata.txt")
io$sce.expr <- paste0(io$basedir,"/processed/rna/SingleCellExperiment.rds")
io$sce.proteomics <- paste0(io$basedir,"/processed/proteomics/SingleCellExperiment.rds")
io$outdir <- file.path(io$basedir,"results/rna_proteomics/rna_vs_protein_levels"); dir.create(io$outdir, showWarnings = F, recursive = T)

# Options
opts$plot_logcounts_rna <- TRUE

###################
## Load metadata ##
###################

sample_metadata_expr.dt <- fread(io$sample_metadata.expr)
sample_metadata_proteomics.dt <- fread(io$sample_metadata.proteomics)

#########################
## Load RNA expression ##
#########################

expr.sce <- readRDS(io$sce.expr)[,sample_metadata_expr.dt$sample]
colnames(expr.sce) <- sample_metadata_expr.dt$alias
# colData(expr.sce) <- sample_metadata_expr.dt %>% tibble::column_to_rownames("sample") %>% DataFrame

proteomics.sce <- readRDS(io$sce.proteomics)[,sample_metadata_proteomics.dt$alias]
# colnames(sce) <- sample_metadata.dt$sample
rownames(proteomics.sce) <- rowData(proteomics.sce)$protein_entry
# colData(proteomics.sce) <- sample_metadata_proteomics.dt %>% tibble::column_to_rownames("sample") %>% DataFrame

assayNames(proteomics.sce) <- c("normalised","unnormalised")

###########################
## Load feature metadata ##
###########################

feature_metadata_proteomics.dt <- fread(io$feature_metadata.proteomics)

###################
## Normalisation ##
###################

# gene expression
if (opts$plot_logcounts_rna) {
  assay(expr.sce,"logcounts") <- log2(1e6*(sweep(assay(expr.sce,"counts"),2,colSums(counts(expr.sce)),"/"))+1)
} else {
  assay(expr.sce,"logcounts") <- 1e6*(sweep(assay(expr.sce,"counts"),2,colSums(counts(expr.sce)),"/"))
}

# proteomics
assay(proteomics.sce,"zscores") <- (assay(proteomics.sce,"normalised")-rowMeans(assay(proteomics.sce,"normalised"),na.rm=T)) / sqrt(rowVars(assay(proteomics.sce,"normalised"),na.rm=T))

################
## Line plots ##
###############

genes.intersect <- intersect(rownames(proteomics.sce), rownames(expr.sce))
"NSD1"%in%genes.intersect
grep("^JARID2",rownames(proteomics.sce),value=T)
grep("^JARID2",rownames(expr.sce),value=T)

# i <- "NSD1"
for (i in genes.to.plot) {
  
  stopifnot(i%in%rownames(expr.sce))
  stopifnot(i%in%rownames(proteomics.sce))
  
  # Expression
  
  to.plot <- data.table(
    alias = colnames(expr.sce),
    expr = assay(expr.sce,"logcounts")[i,]
  ) %>% merge(sample_metadata_expr.dt,by="alias")
  
  # to.plot[,Passage:=factor(Passage,levels=opts$passages)]
  to.plot.means <- to.plot[,.(expr=mean(expr), sd=sd(expr)), by=c("Passage","Treatment")]
  
  ylim.max <- max(5,max(to.plot$expr))
  
  p1 <- ggplot(to.plot, aes(x=Passage, y=expr)) +
    geom_line(aes(color=Treatment), alpha=0.75, data=to.plot.means) +
    geom_point(aes(fill=Treatment), size=4, shape=21, stroke=0.2, color="black") +
    geom_errorbar(aes(ymin=expr-sd, ymax=expr+sd), color="black", width=0.4, alpha=0.9, linewidth=0.5, data=to.plot.means) +
    # scale_fill_manual(values=c(brewer.pal(n=9,"OrRd"),"gray70")) +
    # facet_wrap(~Treatment, nrow=1, scales="fixed") +
    # coord_cartesian(ylim=c(0,ylim.max)) +
    labs(x="Passage number", y=sprintf("RNA expression (%s)",ifelse(opts$plot_logcounts_rna, "log cpm", "cpm")), title=i) +
    theme_classic() +
    theme(
      strip.text = element_text(size=rel(1.15)),
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(colour="black",size=rel(1.10)),
      axis.text.y = element_text(colour="black",size=rel(1.10)),
      legend.position = "top",
      legend.title = element_blank(),
    )
  
  # Proteomics
  
  to.plot <- data.table(
    alias = colnames(proteomics.sce),
    expr = assay(proteomics.sce,"zscores")[i,]
  ) %>% merge(sample_metadata_proteomics.dt,by="alias")
  
  # to.plot[,Passage:=factor(Passage,levels=opts$passages)]
  to.plot.means <- to.plot[,.(expr=mean(expr), sd=sd(expr)), by=c("Passage","Treatment")]
  
  ylim.max <- max(5,max(to.plot$expr))
  
  p2 <- ggplot(to.plot, aes(x=Passage, y=expr)) +
    geom_line(aes(color=Treatment), alpha=0.75, data=to.plot.means) +
    geom_point(aes(fill=Treatment), size=4, shape=21, stroke=0.2, color="black") +
    geom_errorbar(aes(ymin=expr-sd, ymax=expr+sd), color="black", width=0.4, alpha=0.9, linewidth=0.5, data=to.plot.means) +
    # coord_cartesian(ylim=c(0,ylim.max)) +
    labs(x="Passage number", y="Protein abundance (z-score)", title=i) +
    theme_classic() +
    theme(
      strip.text = element_text(size=rel(1.15)),
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(colour="black",size=rel(1.10)),
      axis.text.y = element_text(colour="black",size=rel(1.10)),
      legend.position = "top",
      legend.title = element_blank(),
    )
  
  p <- cowplot::plot_grid(plotlist=list(p1,p2), nrow=1)
  
  pdf(file.path(io$outdir,sprintf("%s_expr_vs_proteomics_lineplots_Chan2022_senescence.pdf",i)), width=10, height=5)
  print(p)
  dev.off()
}


