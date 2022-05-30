#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/utils.R")
}

# I/O
io$outdir <- paste0(io$basedir,"/results/rna/individual_genes/pseudobulk"); dir.create(io$outdir, showWarnings = F)

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

###############
## Load PAGA ##
###############

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/load_paga_graph.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/load_paga_graph.R")
} else {
  stop("Computer not recognised")
}

# Plot graph structure
p <- ggnet2(
  net = net.paga,
  mode = c("x", "y"),
  node.size = 0,
  edge.size = 0.15,
  edge.color = "grey",
  label = FALSE,
  label.size = 2.3
)


##########
## Plot ##
##########

# Define color scale
rna.col.seq <- chromvar.col.seq <- round(seq(0,1,0.1), 2)
rna.colors <- colorRampPalette(c("gray92", "darkgreen"))(length(rna.col.seq))

# Define genes to plot
# genes.to.plot <- rownames(rna.sce)[grep("Gata",rownames(rna.sce))]
genes.to.plot <- c("Foxa2","Tfap2a","Mesp1")

for (i in 1:length(genes.to.plot)) {
  gene <- genes.to.plot[i]
  print(sprintf("%s/%s: %s",i,length(genes.to.plot),gene))
    
  expr.values <- logcounts(rna.sce[gene,])[1,] %>% minmax.normalisation()
  expr.colors <- round(expr.values,1) %>% map(~ rna.colors[which(rna.col.seq == .)]) %>% unlist
  
  p.rna <- p + geom_text(label = "\u25D0", aes(x=x, y=y), color=expr.colors, size=20, family = "Arial Unicode MS",
                      data = p$data[,c("x","y")] %>% dplyr::mutate(expr=expr.colors)) +
    scale_colour_manual(values=expr.colors) + 
    labs(title=gene) +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
  
  
  png(sprintf("%s/%s_rna_expression_paga.png",io$outdir,gene), width = 350, height = 400)
  # pdf(sprintf("%s/%s_rna_expression_paga.pdf",io$outdir,i), width=5, height=3.5)
  print(p.rna)
  dev.off()
}


