source(here::here("settings.R"))
source(here::here("utils.R"))

library(scater)

#####################
## Define settings ##
#####################

# I/O
io$sce <- file.path(io$basedir,"results/rna/metacells/SingleCellExperiment_metacells.rds")
io$metadata <- file.path(io$basedir,"results/rna/metacells/metacells_metadata.txt.gz")
io$outdir <- paste0(io$basedir,"/results/rna/gene_expression/metacells")

# Options
opts$aliases <- c( 
  "D0",
  "D2",
  "D8",
  "D12",
  "D16_TRA-1-60_negative",
  "D16_TRA-1-60_positive"
)

########################
## Load cell metadata ##
########################

cell_metadata.dt <- fread(io$metadata) %>%
  .[,sample:=NULL] %>% setnames("alias","sample") %>%
  .[sample%in%opts$aliases] %>% .[,sample:=factor(sample,levels=opts$aliases)] %>%
  .[,stage:=factor(stage, levels=opts$stages, ordered = T)] %>%
  setnames("metacell","cell")

table(cell_metadata.dt$sample)

###############################
## Load SingleCellExperiment ##
###############################

sce <- load_SingleCellExperiment(io$sce, cells=cell_metadata.dt$cell, normalise = TRUE)
colData(sce) <- cell_metadata.dt %>% tibble::column_to_rownames("cell") %>% DataFrame

##########################
## Plot gene expression ##
##########################

# genes.to.plot <- fread("/Users/rargelaguet/data/Xing2020_bj_oskm/results/differential/diff_snClust_mix_vs_fm5.txt.gz") %>% 
#   .[abs(logFC)>=2 & transcript_type=="gene",gene]

give.n <- function(x, ymax=NULL){
  if (!is.null(ymax)) {
    return( data.frame(y = ymax, label = paste0("n=",length(x)) ) ) 
  } else {
    return( data.frame(y = max(x), label = paste0("n=",length(x)) ) ) 
  }
}
  
# i <- "ELN" 
for (i in genes.to.plot) {
  
  to.plot <- data.table(
    cell = colnames(sce),
    expr = logcounts(sce)[i,]
  ) %>% merge(cell_metadata.dt,by="cell")
  
  p <- ggplot(to.plot, aes(y=sample, x=expr)) +
    geom_boxplot(aes(fill=sample), outlier.shape=NA) +
    # stat_summary(fun.data = give.n, fun.args = list(ymax = max(to.plot$expr)+1), size=2.5, geom = "text", position = position_dodge(width=0.75)) +
    # ggrastr::geom_point_rast(aes(fill=alias), shape=21, size=0.75, stroke=0.05, alpha=0.80, raster.dpi=150) +
    # scale_fill_gradient(low = "gray80", high = "purple") +
    # facet_wrap(~snClust, scales = "fixed") +
    # guides(x = guide_axis(angle = 90)) +
    labs(y="", x="RNA expression", title=i) +
    theme_classic() +
    # ggplot_theme_NoAxes() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.y = element_text(color="black", size=rel(1.10)),
      axis.text = element_text(color="black", size=rel(1.0)),
      legend.title = element_blank(),
      legend.position = "none"
    )
  
  pdf(file.path(io$outdir,sprintf("%s_expr_boxplots_Xing2020_cells.pdf",i)), width=6, height=5)
  print(p)
  dev.off()
}

