#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

# I/O
# io$rna.pseudobulk.sce <- paste0(io$basedir,"/results/rna/pseudobulk/celltype.mapped/SingleCellExperiment_pseudobulk.rds")
io$outdir <- paste0(io$basedir,"/results/rna/plot_individual_genes/celltypes/pseudobulk"); dir.create(io$outdir, showWarnings = F)

# Options
opts$celltypes = c(
	"Epiblast",
	"Primitive_Streak",
	"Caudal_epiblast",
	"PGC",
	"Anterior_Primitive_Streak",
	"Notochord",
	"Def._endoderm",
	"Gut",
	"Nascent_mesoderm",
	"Mixed_mesoderm",
	"Intermediate_mesoderm",
	"Caudal_Mesoderm",
	"Paraxial_mesoderm",
	"Somitic_mesoderm",
	"Pharyngeal_mesoderm",
	"Cardiomyocytes",
	"Allantois",
	"ExE_mesoderm",
	"Mesenchyme",
	"Haematoendothelial_progenitors",
	"Endothelium",
	"Blood_progenitors_1",
	"Blood_progenitors_2",
	"Erythroid1",
	"Erythroid2",
	"Erythroid3",
	"NMP",
	"Rostral_neurectoderm",
	"Caudal_neurectoderm",
	"Neural_crest",
	"Forebrain_Midbrain_Hindbrain",
	"Spinal_cord",
	"Surface_ectoderm",
	"Visceral_endoderm",
	"ExE_endoderm",
	"Parietal_endoderm",
	"ExE_ectoderm"
)

###############
## Load data ##
###############

# Load SingleCellExperiment object
rna.sce <- readRDS(io$rna.pseudobulk.sce)[,opts$celltypes]

# Load gene metadata
gene_metadata <- fread(io$gene_metadata) %>%
  .[symbol%in%rownames(rna.sce)]

####################################
## Plot (style 1: lines and dots) ##
####################################

# genes.to.plot <- c("Lefty1","Cd34","Tmsb4x","Fgf3","Spata7","Cer1","Spink1","Dppa4","Dppa5a","Prc1","Lefty2","Ube2c","Hba-x","Hbb-y","Hba-a1","Hbb-bh1")
# genes.to.plot <- c("Vegfa","Vegfb","Vegfc","Vegfd","Kdr","Flt1","Tal1","Runx1","Etv2)
# genes.to.plot <- c("Tet1","Tet2","Tet3","Dnmt1","Dnmt3a","Dnmt3b","Dnmt3l")
# genes.to.plot <- rownames(rna.sce)[grep("Zscan",rownames(rna.sce))]
genes.to.plot <- c("Gata6")

celltype.order <- opts$celltypes

for (i in 1:length(genes.to.plot)) {
  gene <- genes.to.plot[i]
  print(sprintf("%s/%s: %s",i,length(genes.to.plot),gene))
    
  to.plot <- data.table(
    celltype = colnames(rna.sce),
    expr = logcounts(rna.sce[gene,])[1,]
  ) %>% .[,celltype:=factor(celltype,levels=rev(celltype.order))] %>%
    .[,minmax_expr:=minmax.normalisation(expr)]
  
  p <- ggplot(to.plot, aes(x=celltype, y=expr, fill=celltype)) +
    # geom_bar(stat="identity", width=0.45, alpha=0.9) +
    # geom_point(shape=21, size=3.5) +
    geom_point(aes(size=expr), shape=21) +
    scale_size_continuous(range = c(0.25,5)) + 
    scale_fill_manual(values=opts$celltype.colors) +
    theme_classic() +
    coord_flip(ylim=c(0,10)) +
    labs(x="",y="RNA expression") +
    theme(
      axis.line = element_blank(),
      # axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_text(colour="black",size=rel(1.0)),
      axis.text.y = element_text(colour="black",size=rel(1.0)),
      # axis.text.x = element_blank(),
      axis.title.x = element_text(colour="black",size=rel(1.2)),
      legend.position="none"
    )
  
  pdf(sprintf("%s/%s.pdf",io$outdir,gene), width=4, height=10)
  print(p)
  dev.off()
}

##############################
## Plot (style 2: barplots) ##
##############################

genes.to.plot <- c("Mesp2")

celltype.order <- opts$celltypes

for (i in genes.to.plot) {

  print(sprintf("%s/%s: %s",match(i,genes.to.plot),length(genes.to.plot),i))
  
  to.plot <- data.table(
    celltype = colnames(rna.sce),
    expr = logcounts(rna.sce[i,])[1,]
  ) %>% .[,celltype:=factor(celltype,levels=celltype.order)]
  
  p <- ggplot(to.plot, aes(x=celltype, y=expr, fill=celltype)) +
    geom_bar(stat="identity", color="black") +
    scale_fill_manual(values=opts$celltype.colors) +
    labs(x="",y=sprintf("%s expression",i)) +
    guides(x = guide_axis(angle = 90)) +
    theme_classic() +
    theme(
      axis.text.x = element_text(colour="black",size=rel(1)),
      axis.text.y = element_text(colour="black",size=rel(0.9)),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    )
  
  # pdf(sprintf("%s/%s.pdf",io$outdir,i), width=5, height=3.5, useDingbats = F)
  print(p)
  # dev.off()
}


##########
## Test ##
##########

gene <- "Foxa2"

to.plot <- data.table(
  celltype = colnames(rna.sce),
  expr = logcounts(rna.sce[gene,])[1,]
) %>% .[,celltype:=factor(celltype,levels=rev(celltype.order))] %>%
  .[,expr:=minmax.normalisation(expr)]

p <- ggplot(to.plot, aes(x=celltype, y=expr, fill=celltype)) +
  geom_bar(stat="identity", width=0.45, alpha=0.9) +
  geom_point(shape=21, size=3.5) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme_classic() +
  coord_flip() +
  labs(x="",y="RNA expression") +
  theme(
    axis.line = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    # axis.text.x = element_text(colour="black",size=rel(1.0)),
    axis.text.x = element_blank(),
    axis.title.x = element_text(colour="black",size=rel(1.2)),
    legend.position="none"
  )

pdf(sprintf("%s/test.pdf",io$outdir), width=2.5, height=8)
print(p)
dev.off()