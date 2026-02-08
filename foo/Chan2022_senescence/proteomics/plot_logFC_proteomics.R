source(here::here("settings.R"))
source(here::here("utils.R"))

library(pheatmap)

#####################
## Define settings ##
#####################

# I/O
io$outdir <- file.path(io$basedir,"proteomics/results/protein_levels"); dir.create(io$outdir, showWarnings = F, recursive = T)

# Options

###################
## Load metadata ##
###################

# sample_metadata.dt <- fread(io$metadata)

###############
## Load data ##
###############

# proteomics_logfc.dt <- fread("/Users/rargelaguet/data/Chan2022_senescence/proteomics/proteomics_logFC.txt.csv")
proteomics_logfc.mtx <- fread("/Users/rargelaguet/data/Chan2022_senescence/proteomics/proteomics_logFC.txt.csv") %>% matrix.please
colnames(proteomics_logfc.mtx) <- gsub("_fc","",colnames(proteomics_logfc.mtx))

################
## Line plots ##
###############

all.genes <- rownames(proteomics_logfc.mtx)

genes.to.plot <- grep("^MKI",all.genes,value=T)
genes.to.plot <- grep("CD",all.genes,value=T)

# i <- "MKI67"
for (i in genes.to.plot) {
  
  to.plot <- data.table(
    sample = colnames(proteomics_logfc.mtx),
    feature = i,
    value = proteomics_logfc.mtx[i,]
  )
  to.plot[,class:=as.character(NA)]
  to.plot[grepl("PDL",sample),class:="Replicative senescence"]
  to.plot[grepl("hTERT",sample),class:="hTERT"]
  
  # ylim.max <- max(5,max(to.plot$expr))
  
  p <- ggplot(to.plot, aes(x=sample, y=value)) +
    geom_line(aes(group=class), alpha=0.75) +
    geom_point(size=4, shape=21, stroke=0.2, color="black", fill="gray60") +
    labs(x="", y="Protein abundance (relative to control, logFC)", title=i) +
    geom_hline(yintercept=0, linetype="dashed") +
    facet_wrap(~class, nrow=1, scales="free_x") +
    theme_classic() +
    theme(
      strip.text = element_text(size=rel(1.15)),
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(colour="black",size=rel(1.05)),
      axis.text.y = element_text(colour="black",size=rel(1.05)),
      legend.position = "top",
      legend.title = element_blank(),
    )
  
  pdf(file.path(io$outdir,sprintf("%s_logFC_proteomics_lineplots_Chan2022_senescence.pdf",i)), width=8.5, height=5.5)
  print(p)
  dev.off()
}

