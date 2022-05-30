# https://www.ArchRProject.com/bookdown/calling-peaks-with-archr.html
# Note: this requires the creation of pseudobulk replicates with 'addGroupCoverages'
# see /.../atac/archR/pseudobulk/archR_pseudobulk_celltypes.R

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

#####################
## Define settings ##
#####################

# io$metadata <- paste0(io$basedir,"/sample_metadata.txt.gz")
io$outdir <- paste0(io$basedir,"/results/atac/archR/peak_calling")

opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E8.0_rep1",
  "E8.0_rep2",
  "E8.5_rep1",
  "E8.5_rep2"
)

opts$min.score <- 20

########################
## Load cell metadata ##
########################

# sample_metadata <- fread(io$metadata) %>%
#   .[pass_atacQC==TRUE] %>%
#   .[sample%in%opts$samples]
# stopifnot(sample_metadata$cell %in% rownames(ArchRProject))

##################
## Subset ArchR ##
##################

# ArchRProject <- ArchRProject[sample_metadata$cell]
# table(getCellColData(ArchRProject,"Sample")[[1]])

##########################
## Load peak annotation ##
##########################

peakSet.gr <- readRDS(io$archR.peakSet.granges)

# Filter
peakSet.gr.filt <- peakSet.gr[peakSet.gr$score>opts$min.score,]

# Merge (I haven't tested the code)
# hits <- findOverlaps(x, y)
# xhits <- x[queryHits(hits)]
# yhits <- y[subjectHits(hits)]
# frac <- width(pintersect(xhits, yhits)) / pmin(width(xhits), width(yhits))
# merge <- frac >= minfrac
# c(reduce(c(xhits[merge], yhits[merge])),
#   xhits[!merge], yhits[!merge],
#   x[-queryHits(hits)], y[-subjectHits(hits)])

##################
## Add to ArchR ##
##################

ArchRProject <- addPeakSet(ArchRProject, peakSet.gr.filt, force = T)
length(ArchRProject@peakSet)

##########
## Save ##
##########

# Save PeakSet
saveRDS(peakSet.gr, paste0(io$archR.directory,"/peakSet_original.rds"))
saveRDS(peakSet.gr.filt, io$archR.peakSet.granges)

# fetch peaks in data.table format
dt <- peakSet.gr.filt %>% as.data.table() %>% setnames(c("seqnames"),c("chr"))

# Save peak metadata
outfile <- paste0(io$archR.directory,"/PeakCalls/peak_metadata.tsv.gz")
fwrite(dt, outfile, sep="\t")

# save peaks in bed format
to.save <- dt[,c("chr","start","end")]
fwrite(to.save, paste0(io$archR.directory,"/PeakCalls/peaks_archR_macs2.bed.gz"), sep="\t", col.names = F)

#####################
## Add peak matrix ##
#####################

ArchRProject <- addPeakMatrix(ArchRProject, binarize = FALSE, force = TRUE)

# Save matrix.mtx.gz
# mtx <- getMatrixFromProject(ArchRProject, useMatrix = "PeakMatrix")@assays@data[[1]]
# outfile <- paste0(io$archR.directory,"/PeakMatrix/matrix.mtx.gz")
# Matrix::writeMM(mtx, file=outfile)

# Save barcodes.tsv.gz
# outfile <- paste0(io$archR.directory,"/PeakMatrix/barcodes.tsv.gz")
# fwrite(as.data.table(gsub("#","_", colnames(mtx))), outfile, col.names=F)

# save features.tsv.gz
# foo <- dt[,c("chr","start","end")] %>% .[,foo:=sprintf("%s_%s_%s",chr,start,end)] %>% .[,"foo"]
# outfile <- paste0(io$archR.directory,"/PeakMatrix/features.tsv.gz")
# fwrite(foo, outfile, sep="\t", col.names = F)

########################
## Save ARchR project ##
########################

saveArchRProject(ArchRProject)