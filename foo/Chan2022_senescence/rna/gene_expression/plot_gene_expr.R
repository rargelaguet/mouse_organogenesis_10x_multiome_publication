source(here::here("settings.R"))
source(here::here("utils.R"))

library(pheatmap)

#####################
## Define settings ##
#####################

# I/O
# io$metadata <- file.path(io$basedir,"processed/rna/sample_metadata.txt")
io$sce <- file.path(io$basedir,"rna/processed/SingleCellExperiment.rds")
io$outdir <- file.path(io$basedir,"rna/results/gene_expression"); dir.create(io$outdir, showWarnings = F, recursive = T)

# Options
opts$plot_logcounts <- TRUE

opts$conditions <- c(
  "hTERT",
  "Replicative_Senescence",
  # "Cell_Density",
  "Radiation_Senescence"
)

###################
## Load metadata ##
###################

sample_metadata.dt <- fread(io$metadata) %>%
  .[condition%in%opts$conditions] %>% .[,condition:=factor(condition,levels=opts$conditions)] %>%
  .[time%in%opts$timepoints] %>% .[,time:=factor(time,levels=opts$timepoints)]

#########################
## Load RNA expression ##
#########################

sce <- readRDS(io$sce)[,sample_metadata.dt$sample]
# colnames(sce) <- sample_metadata.dt$sample
rownames(sce) <- rowData(sce)$gene_name

colData(sce) <- sample_metadata.dt %>% tibble::column_to_rownames("sample") %>% DataFrame

########################
## Load gene metadata ##
########################

# gene_metadata.dt <- fread(io$gene_metadata)

io$limma_results_cell_density <- file.path(io$basedir,"rna/results/limma//limma_results_cell_density.txt.gz")
io$limma_results_hTERT <- file.path(io$basedir,"rna/results/limma//limma_results_hTERT.txt.gz")
io$limma_results_radiation_senescence <- file.path(io$basedir,"rna/results/limma//limma_results_radiation_senescence.txt.gz")
io$limma_results_replicative_senescence <- file.path(io$basedir,"rna/results/limma//limma_results_replicative_senescence.txt.gz")
io$limma_results <- file.path(io$basedir,"results/rna/limma/limma_results.rds")

limma_list <- list()
limma_list[["cell_density"]] <- fread(io$limma_results_cell_density) %>% .[is.na(symbol),symbol:=ens_id] %>% setnames("symbol","feature")
limma_list[["hTERT"]] <- fread(io$limma_results_hTERT) %>% .[is.na(symbol),symbol:=ens_id] %>% setnames("symbol","feature")
limma_list[["radiation_senescence"]] <- fread(io$limma_results_radiation_senescence) %>% .[is.na(symbol),symbol:=ens_id] %>% setnames("symbol","feature")
limma_list[["replicative_senescence"]] <- fread(io$limma_results_replicative_senescence) %>% .[is.na(symbol),symbol:=ens_id] %>% setnames("symbol","feature")

# Load gene signatures
source(here::here("rna/gene_signatures/load_gene_signatures.R"))

###################
## Normalisation ##
###################

if (opts$plot_logcounts) {
  assay(sce,"logcounts") <- log2(1e6*(sweep(assay(sce,"counts"),2,colSums(counts(sce)),"/"))+1)
} else { 
  assay(sce,"logcounts") <- 1e6*(sweep(assay(sce,"counts"),2,colSums(counts(sce)),"/"))
}

############################
## Time-course line plots ##
############################

genes.to.plot <- limma_list[["replicative_senescence"]] %>% .[abs(coef)>0.05 & pval_adj<=0.0001,feature]

genes.to.plot <- unique(unlist(gene.sets[grep("senescence",names(gene.sets), ignore.case = TRUE)]))

genes.to.plot <- grep("^JUN|^FOS",rownames(sce),value=T)
genes.to.plot <- c("HSPA1A", "HSPA1B", "PDLIM1", "IGFBP3", "CLIC1", "IQGAP1", "CXCR2", "TP53", "PLEC", "HMGA1", "CDKN1A", "IL6", "ACTB", "CXCL1", "NFKB1", "RELA", "STAT3", "CCL2")
genes.to.plot <- c("TMX4","TMEM123","ITGAV")

genes.to.plot <- genes.to.plot[genes.to.plot%in%rownames(sce)]

# i <- "IL11"
for (i in genes.to.plot) {
  
  to.plot <- data.table(
    sample = colnames(sce),
    expr = assay(sce,"logcounts")[i,]
  ) %>% merge(sample_metadata.dt,by="sample")
  
  to.plot.means <- to.plot[,.(expr=mean(expr), sd=sd(expr)), by=c("condition","time")]
  
  ylim.max <- max(5,max(to.plot$expr))
  
  p <- ggplot(to.plot, aes(x=time, y=expr)) +
    geom_line(aes(group=condition), color="gray40", linewidth=1, alpha=0.75, data=to.plot.means) +
    geom_point(color="gray40", size=2.50) +
    geom_errorbar(aes(ymin=expr-sd, ymax=expr+sd), width=0.4, alpha=0.9, linewidth=0.5, data=to.plot.means) +
    # geom_line(aes(group=condition, color=condition), linewidth=0.75, alpha=0.75, data=to.plot.means) +
    # geom_point(aes(color=condition), size=1.50) +
    # geom_errorbar(aes(ymin=expr-sd, ymax=expr+sd, color=condition), width=0.4, alpha=0.9, linewidth=0.5, data=to.plot.means) +
    # scale_fill_manual(values=c(brewer.pal(n=9,"OrRd"),"gray70")) +
    facet_wrap(~condition, nrow=1, scales="free_x") +
    # coord_cartesian(ylim=c(0,2)) +
    labs(x="Time", y=sprintf("RNA expression (%s)",ifelse(opts$plot_logcounts, "log cpm", "cpm")), title=i) +
    theme_classic() +
    theme(
      strip.text = element_text(size=rel(1.10)),
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(colour="black",size=rel(1.0)),
      axis.title = element_text(colour="black",size=rel(1.10)),
      axis.text.y = element_text(colour="black",size=rel(1.10)),
      legend.position = "none",
      legend.title = element_blank(),
    )
  
  pdf(file.path(io$outdir,sprintf("%s_expr_lineplots_Chan2022_senescence.pdf",i)), width=7, height=5.75)
  print(p)
  dev.off()
}


###############
## Box plots ##
###############

# genes.to.plot <- genes.to.plot[genes.to.plot%in%rownames(sce) ]
# 
# samples.to.plot <- c("p8_Control","p8_TG","p11_Control","p11_TG","p24_Control","p24_TG","p29_Control","p29_TG")
# 
# # i <- "ST6GAL1"
# for (i in genes.to.plot) {
#   
#   to.plot <- data.table(
#     sample = colnames(sce),
#     expr = assay(sce,"logcounts")[i,]
#   ) %>% merge(sample_metadata.dt,by="sample") %>% .[Group%in%samples.to.plot]
#   
#   to.plot.means <- to.plot[,.(expr=mean(expr), sd=sd(expr)), by=c("Passage","Treatment")]
#   
#   p <- ggplot(to.plot, aes(x=Treatment, y=expr, fill=Treatment)) +
#     # geom_bar(stat="identity", color="black", data=to.plot.means) +
#     # geom_errorbar(aes(ymin=signature-sd, ymax=signature+sd), color="black", width=0.4, alpha=0.9, linewidth=0.5, data=to.plot.means) +
#     geom_boxplot() +
#     geom_point(size=3, shape=21, stroke=0.2, color="black") +
#     facet_wrap(~passage_string, nrow=1) +
#     labs(x="", y=sprintf("RNA expression (%s)",ifelse(opts$plot_logcounts, "log cpm", "cpm")), title=i) +
#     theme_classic() +
#     theme(
#       strip.text = element_text(size=rel(1)),
#       plot.title = element_text(hjust = 0.5),
#       axis.text.x = element_text(colour="black",size=rel(1.20)),
#       axis.text.y = element_text(colour="black",size=rel(1.0)),
#       axis.title.y = element_text(colour="black",size=rel(1.15)),
#       legend.position = "none",
#       legend.title = element_blank(),
#     )
#   
#   pdf(file.path(io$outdir,sprintf("%s_expr_boxplots_Chan2022_senescence.pdf",i)), width=7, height=5)
#   print(p)
#   dev.off()
# }
