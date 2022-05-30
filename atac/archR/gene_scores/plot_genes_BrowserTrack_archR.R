
#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/utils.R")
} else {
  stop("Computer not recognised")
}

# I/O
# io$archR.diff.dir <- paste0(io$basedir,"/results/atac/archR/differential/PeakMatrix")
io$outdir <- paste0(io$basedir,"/results/atac/archR/BrowserTrack/genes"); dir.create(io$outdir, showWarnings = F)

# Options
opts$celltypes <- c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  # "PGC",
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
  "ExE_ectoderm",
  "Parietal_endoderm"
)

opts$celltypes.to.rename <- c(
  "Anterior_Primitive_Streak" = "Primitive_Streak",
  "Erythroid1" = "Erythroid",
  "Erythroid2" = "Erythroid",
  "Erythroid3" = "Erythroid",
  "Blood_progenitors_1" = "Blood_progenitors",
  "Blood_progenitors_2" = "Blood_progenitors",
  "Caudal_neurectoderm" = "Neurectoderm",
  "Rostral_neurectoderm" = "Neurectoderm",
  # "Parietal_endoderm" = "ExE_endoderm",
  # "Visceral_endoderm" = "ExE_endoderm",
  # "Allantois" = "ExE_mesoderm",
  "Mixed_mesoderm" = "Nascent_mesoderm"
  # "Paraxial_mesoderm" = "Somitic_mesoderm"
)


########################
## Load ArchR project ##
########################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
} else {
  stop("Computer not recognised")
}

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & doublet_call==FALSE] %>%
  .[sample%in%opts$samples & celltype.predicted%in%opts$celltypes] %>%
  .[,celltype.predicted:=stringr::str_replace_all(celltype.predicted,opts$celltypes.to.rename)]
  # .[,celltype.predicted:=factor(celltype.predicted,levels=opts$celltypes)]

#########################
## Subset ArchR object ##
#########################

ArchRProject.filt <- ArchRProject[sample_metadata$cell]

# add celltype.predicted to ArchR's CellColData
tmp <- sample_metadata %>% 
  .[cell%in%rownames(ArchRProject.filt)] %>% setkey(cell) %>% .[rownames(ArchRProject.filt)] %>%
  as.data.frame() %>% tibble::column_to_rownames("cell")
stopifnot(all(tmp$TSSEnrichment_atac == getCellColData(ArchRProject.filt, "TSSEnrichment")[[1]]))
ArchRProject.filt <- addCellColData(
  ArchRProject.filt,
  data = tmp$celltype.predicted, 
  name = "celltype.predicted",
  cells = rownames(tmp),
  force = TRUE
)

#######################
## Load marker genes ##
#######################

marker_genes.dt <- fread(io$rna.atlas.marker_genes)

# GenomicRanges
genes.gr <- getGeneAnnotation(ArchRProject.filt)[["genes"]]
genes.gr <- genes.gr[genes.gr$symbol%in%unique(marker_genes.dt$gene)]

####################
## Browser tracks ##
####################

genes.to.plot <- genes.gr$symbol %>% head(n=5) %>% unname
# i <- "chr2:39483639-39484239"

# opts$extend.upstream <- 1e4
# opts$extend.downstream <- 1e4
# opts$tileSize <- 50

opts$tileSize <- 30

# Ugly hack
celltype.order = c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
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
  "Blood_progenitors",
  "Erythroid",
  "NMP",
  "Neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm",
  "Visceral_endoderm",
  "ExE_endoderm",
  "Parietal_endoderm",
  "ExE_ectoderm"
)
stopifnot(sort(celltype.order)==sort(unique(ArchRProject.filt$celltype.predicted)))
rename <- paste(1:length(celltype.order),celltype.order,sep="_")
names(rename) <- celltype.order
ArchRProject.filt$celltype.predicted2 <- stringr::str_replace_all(ArchRProject.filt$celltype.predicted,rename)
opts$celltype.colors2 <- opts$celltype.colors[names(opts$celltype.colors)%in%unique(ArchRProject.filt$celltype.predicted)]
names(opts$celltype.colors2) <- stringr::str_replace_all(names(opts$celltype.colors2),rename)

for (i in genes.to.plot) {
  
  to.plot.gr <- genes.gr[genes.gr$symbol==i]
  gene.length <- abs(end(to.plot.gr) - start(to.plot.gr))
  start(to.plot.gr) <- start(to.plot.gr) - round(gene.length/1.5)
  end(to.plot.gr) <- end(to.plot.gr) + round(gene.length/1.5)
  
  # Plot
  p <- plotBrowserTrack(
    ArchRProj = ArchRProject.filt, 
    region = to.plot.gr,
    geneSymbol = i,
    # useMatrix = "GeneScoreMatrix",
    groupBy = "celltype.predicted2", 
    tileSize = opts$tileSize,
    # upstream = opts$extend.upstream,
    # downstream = opts$extend.downstream,
    pal = opts$celltype.colors2,
    plotSummary = c("bulkTrack", "featureTrack", "geneTrack"),
    sizes = c(13, 1, 1),
  )
  
  pdf(sprintf("%s/test_%s_BrowserTrack.pdf",io$outdir,i), width = 9, height = 5)
  grid::grid.draw(p[[1]])
  dev.off()
}
