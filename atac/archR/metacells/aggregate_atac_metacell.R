here::i_am("atac/archR/metacells/aggregate_atac_metacell.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(ArchR))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--archr_directory',     type="character", help='ArchR directory')
p$add_argument('--metadata',    type="character",  help='Cell metadata file')
p$add_argument('--matrices',     type="character",       nargs="+",   help='Matrices to pseudobulk')
p$add_argument('--cell2metacell',    type="character",  nargs="+", help='Metacell results')
p$add_argument('--metacell_min_frags',     type="integer",    default=1e5,     help='Minimum number of fragments per metacell')
p$add_argument('--outdir',          type="character",               help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

###################
## Load settings ##
###################

# START TEST ##
# args$archr_directory <- file.path(io$basedir,"processed/atac/archR")
# args$metadata <- file.path(io$basedir,"results/atac/archR/celltype_assignment/sample_metadata_after_celltype_assignment.txt.gz")
# # args$cell2metacell <- file.path(io$basedir,sprintf("results/rna/metacells/%s/cell2metacell_assignment.txt.gz",opts$samples))
# args$cell2metacell <- file.path(io$basedir,"results/rna/metacells/trajectories/nmp/cell2metacell_assignment.txt.gz")
# args$matrices <- c("PeakMatrix") # c("PeakMatrix", "GeneScoreMatrix_distal", "GeneScoreMatrix_TSS")
# args$metacell_min_frags <- 2e5
# args$outdir <- file.path(io$basedir,"results/atac/archR/metacells/trajectories/nmp/test")
## END TEST ##

# START TEST ##
# args$archr_directory <- file.path(io$basedir,"processed/atac/archR")
# args$metadata <- file.path(io$basedir,"results/atac/archR/celltype_assignment/sample_metadata_after_celltype_assignment.txt.gz")
# args$cell2metacell <- file.path(io$basedir,sprintf("results/rna/metacells/all_cells/%s/cell2metacell_assignment.txt.gz",opts$samples))
# args$matrices <- "GeneScoreMatrix_TSS" # c("PeakMatrix", "GeneScoreMatrix_distal", "GeneScoreMatrix_TSS")
# args$metacell_min_frags <- 2e5
# args$outdir <- file.path(io$basedir,sprintf("results/atac/archR/metacells/all_cells/%s/test",args$matrices))
## END TEST ##

stopifnot(file.exists(args$cell2metacell))
dir.create(args$outdir, showWarnings = F, recursive = T)

###########################
## Load metacell results ##
###########################

cell2metacell.dt <- args$cell2metacell %>% map(~ fread(.)) %>% rbindlist

###################
## Load metadata ##
###################

sample_metadata <- fread(args$metadata) %>%
  .[pass_atacQC==TRUE] %>% # temporary 
  # .[cell%in%cell2metacell.dt$cell] %>%
  # .[,c("cell","sample","stage","genotype","celltype.mapped","closest.cell")] %>%
  merge(cell2metacell.dt,"cell")

print(sprintf("Number of metacells before filtering by minimum number of fragments: %s", length(unique(sample_metadata$metacell))))

######################
## Filter metacells ##
######################

nfrags_per_metacell.dt <- sample_metadata[,.(nfrags=sum(nFrags_atac)),by="metacell"]

# Plot
p <- gghistogram(nfrags_per_metacell.dt[nfrags<=5e6], x="nfrags", y="..density..", bins=100, fill="gray70", color="gray50") +
  geom_vline(xintercept=args$metacell_min_frags, linetype="dashed") +
  labs(x="Number of ATAC fragments per metacell", y="Density") +
  theme(
    axis.text =  element_text(size=rel(0.8)),
    axis.title.x = element_blank()
  )

pdf(file.path(args$outdir,"nfrags_metacell_threshold.pdf"), width=8, height=4)
print(p)
dev.off()

metacells.to.use <- intersect(sample_metadata$cell, nfrags_per_metacell.dt[nfrags>=args$metacell_min_frags,metacell])
sample_metadata <- sample_metadata[metacell%in%metacells.to.use]
cell2metacell.dt <- cell2metacell.dt[metacell%in%metacells.to.use]

print(sprintf("Number of metacells after filtering by total number of fragments: %s", length(unique(sample_metadata$metacell))))

########################
## Load ArchR Project ##
########################

setwd(args$archr_directory)

addArchRGenome("mm10")
addArchRThreads(threads = 1) 

ArchRProject <- loadArchRProject(args$archr_directory)

# Load ArchR projectMetadata
# ArchRProject@projectMetadata <- readRDS(io$archR.projectMetadata)

# Load peaks
# ArchRProject <- addPeakSet(ArchRProject, peakSet = readRDS(args$peak_calls), force = TRUE)
# names(ArchRProject@peakSet) <- sprintf("%s:%s-%s",seqnames(ArchRProject@peakSet), start(ArchRProject@peakSet), end(ArchRProject@peakSet))

# Subset metacells
ArchRProject <- ArchRProject[sample_metadata$cell]

# Sanity checks
stopifnot(args$matrices%in%getAvailableMatrices(ArchRProject))

#############################################
## Add metacell assignment to ArchR object ##
#############################################

tmp <- sample_metadata %>% 
  .[cell%in%rownames(ArchRProject)] %>% setkey(cell) %>% .[rownames(ArchRProject)] %>%
  as.data.frame() %>% tibble::column_to_rownames("cell")
stopifnot(all(rownames(tmp) == rownames(getCellColData(ArchRProject))))
ArchRProject <- addCellColData(
  ArchRProject,
  data = tmp$metacell,
  name = "metacell",
  cells = rownames(tmp),
  force = TRUE
)

# table(ArchRProject$metacell)

###################################################
## Pseudobulk into a SummarizedExperiment object ##
###################################################

# i <- "PeakMatrix"
for (i in args$matrices) {
  print(sprintf("Calculating pseudobulk matrix for %s",i))
  
  # summarise
  atac_metacells.se <- getGroupSE(ArchRProject, groupBy = "metacell", useMatrix = i, divideN = FALSE)
  
  # rename features
  if (grepl("peak",tolower(i),ignore.case=T)) {
    rownames(atac_metacells.se) <- rowData(atac_metacells.se) %>% as.data.table %>% .[,idx:=sprintf("%s:%s-%s",seqnames,start,end)] %>% .$id
  }
  if (grepl("gene",tolower(i),ignore.case=T)) {
    rownames(atac_metacells.se) <- rowData(atac_metacells.se)$name
  }
  
  # TO-DO: ADD CELL METADATA LIKE CELLTYPE
  # colData(atac_metacells.se) <- metacell_metadata.dt %>% tibble::column_to_rownames("metacell") %>% DataFrame
  
  # save
  saveRDS(atac_metacells.se, file.path(args$outdir,sprintf("%s_summarized_experiment_metacells.rds",i)))
}

##############################
## Define metacell metadata ##
##############################

metacell_metadata.dt <- sample_metadata %>%
  .[cell%in%colnames(atac_metacells.se)] %>% setkey(cell) %>% .[colnames(atac_metacells.se)] %>%
  .[,c("metacell", "cell", "sample", "stage", "genotype",  "celltype", "closest.cell")] %>%
  merge(cell2metacell.dt[,.(ncells=.N),by="metacell"],by="metacell")

# Add QC stats
tmp <- data.table(
  metacell = colnames(atac_metacells.se), 
  nFrags_atac = colSums(assay(atac_metacells.se))
)
metacell_metadata.dt <- metacell_metadata.dt %>% merge(tmp,by="metacell")

####################################
## Plot coverage before and after ##
####################################

to.plot <- rbind(
  metacell_metadata.dt[,c("metacell","nFrags_atac")] %>% .[,class:="metacell"],
  sample_metadata[,c("metacell","nFrags_atac")] %>% .[,class:="cell"]
) %>% .[,log_nFrags_atac:=log10(nFrags_atac)] 

p <- gghistogram(to.plot, x="log_nFrags_atac", fill="class", bins=100) +
  labs(x="Number of ATAC fragments (log10)") +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size=rel(0.75), color="black")
  )

pdf(file.path(args$outdir,"nfrags_before_vs_after.pdf"), width=7, height=4)
print(p)
dev.off()

##########
## Save ##
##########

fwrite(metacell_metadata.dt, file.path(args$outdir,"metacells_metadata.txt.gz"), sep="\t", quote=F, na="NA")

##########
## TEST ##
##########

# atac.se <- readRDS("/Users/argelagr/data/gastrulation_multiome_10x/results/atac/archR/metacells/trajectories/nmp/PeakMatrix_summarized_experiment_metacells.rds")
# hist(assay(atac.se))
# 
# ncells <- atac.se@colData$nCells
# 
# # Plot number of cells vs accessibility before normalisation
# to.plot <-data.table(
#   foo = colSums(assay(atac.se)),
#   ncells = ncells
# )
# ggscatter(to.plot, x="foo", y="ncells")
# 
# # depth normalisation
# atac.mtx <- t(t(assay(atac.se)) / colSums(assay(atac.se))) * 1e6
# 
# # divide by number of cells
# atac.mtx <- t(t(assay(atac.se)) / as.vector(ncells))
# 
# # test normalisation
# # atac.mtx <- log2(t(t(assay(atac.se)) / colSums(assay(atac.se))) * 1e6)
# atac.mtx <- log2(1e6*(sweep(assay(atac.se),2,colSums(assay(atac.se)),"/"))+0.5)
# 
# # plot histogram
# hist(atac.mtx)
# 
# # Plot number of cells vs accessibility after normalisation
# to.plot <- data.table(
#   nreads = colSums(atac.mtx),
#   ncells = ncells
# )
# ggscatter(to.plot, x="nreads", y="ncells")
# 
# 
# # plot mean vs variance
# to.plot <-data.table(
#   mean = rowMeans(atac.mtx),
#   var = rowVars(atac.mtx)
# ) %>% .[sample.int(nrow(.), size=nrow(.)/3)]
# ggscatter(to.plot, x="mean", y="var")
