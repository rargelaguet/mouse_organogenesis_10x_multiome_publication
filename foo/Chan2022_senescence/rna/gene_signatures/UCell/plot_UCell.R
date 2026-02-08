source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(UCell))

#####################
## Define settings ##
#####################

# I/O
# io$metadata <- file.path(io$basedir,"processed/rna/sample_metadata.txt")
# io$sce <- paste0(io$basedir,"/processed/rna/SingleCellExperiment.rds")
io$outdir <- file.path(io$basedir,"rna/results/gene_signatures/UCell"); dir.create(io$outdir, showWarnings = F, recursive = T)

# Options
opts$conditions <- c(
  "hTERT",
  "Replicative_Senescence",
  "Cell_Density",
  "Radiation_Senescence"
)

##########################
## Load sample metadata ##
##########################

sample_metadata.dt <- fread(io$metadata) %>%
  .[condition%in%opts$conditions] %>% .[,condition:=factor(condition,levels=opts$conditions)] %>%
  .[time%in%opts$timepoints] %>% .[,time:=factor(time,levels=opts$timepoints)] %>% droplevels

########################
## Load gene metadata ##
########################

gene_metadata.dt <- fread(io$gene_metadata)

###############
## Load data ##
###############

sce <- readRDS(io$sce)[,sample_metadata.dt$sample]
colData(sce) <- sample_metadata.dt %>% tibble::column_to_rownames("sample") %>% DataFrame

# Load precomputed UCell scores
# altExp(sce,'UCell') <- readRDS(file.path(io$outdir,"UCell_SingleCellExperiment.rds"))[,sample_metadata.dt$sample]
altExp(sce,'UCell') <- readRDS(file.path(io$outdir,"UCell_SingleCellExperiment_senescence.rds"))[,sample_metadata.dt$sample]
# altExp(sce,'UCell') <- readRDS(file.path(io$outdir,"UCell_SingleCellExperiment_Yamanaka_regulons.rds"))[,sample_metadata.dt$sample]

##########################
## Load gene signatures ##
##########################

# source(here::here("rna/gene_signatures/load_gene_signatures.R"))

#############################
## Lineplots of signatures ##
#############################

gene.sets.to.plot <- rownames(altExp(sce, 'UCell'))
# gene.sets.to.plot <- grep("CREBBP",rownames(assay(altExp(sce, 'UCell'))), value=T)

# i <- gene.sets.to.plot[1]
for (i in gene.sets.to.plot) {
  
  to.plot <- data.table(
    sample = colnames(sce),
    signature = assay(altExp(sce,'UCell'))[i,]
  ) %>% merge(sample_metadata.dt,by="sample")
  
  # to.plot[,Passage:=factor(Passage,levels=opts$passages)]
  to.plot.means <- to.plot[,.(signature=mean(signature), sd=sd(signature)), by=c("condition","time")]
  
  # ylim.max <- max(5,max(to.plot$signature))
  
  p <- ggplot(to.plot, aes(x=time, y=signature)) +
    geom_line(aes(group=condition), color="gray40", linewidth=1, alpha=0.75, data=to.plot.means) +
    geom_point(color="gray40", size=2.50) +
    geom_errorbar(aes(ymin=signature-sd, ymax=signature+sd), width=0.4, alpha=0.9, linewidth=0.5, data=to.plot.means) +
    # scale_fill_manual(values=c(brewer.pal(n=9,"OrRd"),"gray70")) +
    facet_wrap(~condition, nrow=2, scales="free_x") +
    labs(x="Timepoint", y="UCell signature", title=i) +
    theme_classic() +
    theme(
      strip.text = element_text(size=rel(1.15)),
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(colour="black",size=rel(1.10)),
      axis.text.y = element_text(colour="black",size=rel(1.10)),
      legend.position = "top",
      legend.title = element_blank(),
    )
  
  pdf(file.path(io$outdir,sprintf("%s_signature_lineplots_Chan2022.pdf",gsub("[ /#]","_",i))), width=9, height=7)
  print(p)
  dev.off()
}
 
############################
## Boxplots of signatures ##
############################

gene.sets.to.plot <- rownames(altExp(sce, 'UCell'))
gene.sets.to.plot <- grep("hr3",rownames(assay(altExp(sce, 'UCell'))), value=T)
gene.sets.to.plot <- grep("reactome#G1",rownames(assay(altExp(sce, 'UCell'))), value=T)

samples.to.plot <- c("p8_Control","p8_TG","p11_Control","p11_TG","p24_Control","p24_TG","p29_Control","p29_TG")

# i <- gene.sets.to.plot[2]
for (i in gene.sets.to.plot) {
  
  to.plot <- data.table(
    sample = colnames(sce),
    signature = assay(altExp(sce,'UCell'))[i,]
  ) %>% merge(sample_metadata.dt,by="sample") %>%
    .[Group%in%samples.to.plot]
  
  # to.plot[,Passage:=factor(Passage,levels=opts$passages)]
  to.plot.means <- to.plot[,.(signature=mean(signature), sd=sd(signature)), by=c("passage_string","Treatment")]
  
  # ylim.max <- max(5,max(to.plot$signature))
  
  p <- ggplot(to.plot, aes(x=Treatment, y=signature, fill=Treatment)) +
    # geom_bar(stat="identity", color="black", data=to.plot.means) +
    # geom_errorbar(aes(ymin=signature-sd, ymax=signature+sd), color="black", width=0.4, alpha=0.9, linewidth=0.5, data=to.plot.means) +
    geom_boxplot() +
    geom_point(size=3, shape=21, stroke=0.2, color="black") +
    facet_wrap(~passage_string, nrow=1) +
    labs(x="", y="Pathway transcriptional score", title=strsplit(i,split="#") %>% map_chr(2) %>% gsub("_UCell","",.)) +
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
  
  pdf(file.path(io$outdir,sprintf("%s_signature_boxplots_Chan2022_senescence.pdf",gsub("/","-",i))), width=7, height=5)
  print(p)
  dev.off()
}

