# https://www.ArchRProject.com/bookdown/identifying-marker-peaks-with-archr.html

########################
## Load ArchR Project ##
########################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
} else {
  stop("Computer not recognised")
}

#####################
## Define settings ##
#####################

# I/O
# io$metadata <- paste0(io$basedir,"/processed/atac/archR/sample_metadata_after_archR.txt.gz")
io$outdir <- paste0(io$basedir,"/results/atac/archR/marker_peaks/all_peaks")

# Options
opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E8.5_rep1",
  "E8.5_rep2"
)

opts$celltypes <- c(
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
  "ExE_ectoderm",
  "Parietal_endoderm"
)

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE] %>%
  .[sample%in%opts$samples & celltype.predicted%in%opts$celltypes]
stopifnot(sample_metadata$archR_cell %in% rownames(ArchRProject))

# subset celltypes with sufficient number of cells
opts$min.cells <- 100
sample_metadata <- sample_metadata %>%
  .[,N:=.N,by=c("celltype.predicted")] %>% .[N>opts$min.cells] %>% .[,N:=NULL]

##################
## Subset ArchR ##
##################

ArchRProject.filt <- ArchRProject[sample_metadata$archR_cell,]
table(getCellColData(ArchRProject.filt,"Sample")[[1]])
table(getCellColData(ArchRProject.filt,"celltype.predicted")[[1]])

###################
## Sanity checks ##
###################

getAvailableMatrices(ArchRProject.filt)
length(getPeakSet(ArchRProject.filt))

######################################
## Identify celltype-specific peaks ##
######################################

markersPeaks <- getMarkerFeatures(
  ArchRProj = ArchRProject.filt, 
  useMatrix = "PeakMatrix", 
  groupBy = "celltype.predicted",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
saveRDS(markersPeaks, paste0(io$outdir,"/markersPeaks_summarized_experiment.rds"))
markersPeaks <- readRDS(paste0(io$outdir,"/markersPeaks_summarized_experiment.rds"))

####################
## Prepare output ##
####################

# returns a SummarizedExperiment object
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= Inf & Log2FC >= -Inf")
# markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
names(markerList)
head(markerList[[1]])

dt <- names(markerList) %>%
  map(function(i) as.data.table(markerList[[i]]) %>% .[,celltype:=i]) %>% 
  rbindlist %>%
  setnames("seqnames","chr") %>%
  .[,c("strand","width","idx"):=NULL]
head(dt)

# save 
fwrite(dt, paste0(io$outdir,"/markersPeaks.tsv.gz"), sep="\t")


################
## START TEST ##
################

R.utils::sourceDirectory("/Users/ricard/git/ArchR/R/", verbose=T, modifiedOnly=FALSE)

# peaks.bed <- fread(io$archR.peakSet.bed) %>% .[,foo:=sprintf("%s_%s_%s",V1,V2,V3)]
# foo <- sprintf("%s_%s_%s",seqnames(ArchRProject.filt@peakSet),start(ArchRProject.filt@peakSet),end(ArchRProject.filt@peakSet))

ArchRProj = ArchRProject.filt
groupBy = "celltype.predicted"
useGroups = NULL
bgdGroups = NULL
useMatrix = "PeakMatrix"
bias = c("TSSEnrichment","log10(nFrags)")
normBy = NULL
testMethod = "wilcoxon"
maxCells = 500
scaleTo = 10^4
threads = 1
k = 100
bufferRatio = 0.8
binarize = FALSE
useSeqnames = NULL
verbose = TRUE
logFile = createLogFile("getMarkerFeatures")

##############
## END TEST ##
##############