source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$sce <- file.path(io$basedir,"results/pseudobulk/SingleCellExperiment_pseudobulk.rds")
io$outdir <- file.path(io$basedir,"results/gene_expression/pseudobulk"); dir.create(io$outdir, showWarnings = F, recursive = T)

# Options
opts$groups <- c(
  "hTert_t2", "hTert_t4", "hTert_t5", "hTert_t6", "hTert_t7", 
  "PDL_25", "PDL_29", "PDL_33", "PDL_37", "PDL_46", "PDL_50"
)

###################
## Load metadata ##
###################

sample_metadata.dt <- fread(io$sample_metadata) %>%
  .[,sample:=NULL] %>% setnames("alias","sample") %>%
  .[group%in%opts$groups] %>% .[,group:=factor(group, levels=opts$groups)] 

#########################
## Load RNA expression ##
#########################

sce <- readRDS(io$sce)[,sample_metadata.dt$sample]

# update colData
colData(sce) <- sample_metadata.dt %>% tibble::column_to_rownames("sample") %>% DataFrame

# update rowData
# rownames(sce) <- rowData(sce)$external_gene_name
# sce <- sce[!duplicated(rownames(sce)) & !is.na(rownames(sce)),]

########################
## Load gene metadata ##
########################

gene_metadata.dt <- fread(io$gene_metadata)

##########################
## Load gene signatures ##
##########################

# source(here::here("gene_signatures/load_gene_signatures.R"))
# gene_signatures_list_generic <- gene_signatures_list

# source(here::here("gene_signatures/load_manual_gene_signatures.R"))
# gene_signatures_list_manual <- gene_signatures_list

###################
## Normalisation ##
###################

# assay(sce,"cpm") <- 1e6*(sweep(assay(sce,"counts"),2,colSums(assay(sce,"counts")),"/"))
# assay(sce,"log_cpm") <- log2(1e6*(sweep(assay(sce,"counts"),2,colSums(assay(sce,"counts")),"/"))+1)

##########################
## Plot gene expression ##
##########################

assay.to.plot <- "log_cpm"

genes.to.plot <- c("CD70","PI16","CD13","KLF4")
genes.to.plot <- grep("_gene",rownames(sce), value=T)

genes.to.plot <- genes.to.plot[genes.to.plot%in%rownames(sce)]

# i <- "ACTA2"
for (i in genes.to.plot) {
  
  to.plot <- data.table(
    sample = colnames(sce),
    expr = assay(sce,assay.to.plot)[i,]
  ) %>% merge(sample_metadata.dt,by="sample")

  if (assay.to.plot == "log_cpm") {
    ylim.max.default <- 5
  } else if (assay.to.plot == "cpm") {
    ylim.max.default <- 10
  }

  p <- ggplot(to.plot, aes(x=day, y=expr, fill=genetic_background)) +
    geom_point(color="black", alpha=0.90, size=5.5, shape=21, stroke=0.2) +
    # geom_line(color="black", alpha=0.90, size=5.5, shape=21, stroke=0.2) +
    # facet_wrap(~treatment, nrow=1) +
    # facet_grid(~treatment, scales = "free_x", space = "free_x") +
    labs(x="Population doublings", y=sprintf("Gene expression (%s)",assay.to.plot), title=sprintf("%s", i)) +
    coord_cartesian(ylim=c(NA,max(ylim.max.default,max(to.plot$expr)))) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(colour="black",size=rel(1)),
      # axis.text.x = element_blank(),
      axis.text.y = element_text(colour="black",size=rel(1)),
      legend.position = "right",
      legend.title = element_blank()
    )

  pdf(file.path(io$outdir,sprintf("%s_expr_pseudobulk_Chan2022_senescence.pdf",i)), width=7, height=5.5)
  print(p)
  dev.off()
}
