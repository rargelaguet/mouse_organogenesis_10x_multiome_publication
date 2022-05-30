#####################
## Define settings ##
#####################

here::i_am("atac/archR/processing/update_archR_metadata.R")

source(here::here("settings.R"))
source(here::here("utils.R"))


########################
## Load ArchR project ##
########################

source(here::here("atac/archR/load_archR_project.R"))

############################################
## Merge archR metadata with RNA metadata ##
############################################

# Fetch pre-computed archR's metadata
io$archr.metadata <- paste0(io$basedir,"/processed/atac/archR/sample_metadata_after_archR.txt.gz")
archr_metadata <- fread(io$archr.metadata)
stopifnot(all(rownames(ArchRProject) %in% archr_metadata$cell))
# cols.to.rename <- c("TSSEnrichment","ReadsInTSS","PromoterRatio","NucleosomeRatio","nFrags","BlacklistRatio")
# idx.cols.to.rename <- which(colnames(archr_metadata)%in%cols.to.rename)
# colnames(archr_metadata)[idx.cols.to.rename] <- paste0(colnames(archr_metadata)[idx.cols.to.rename], "_atac")

# Fetch the metadata file of interest
io$updated.metadata <- paste0(io$basedir,"/sample_metadata.txt.gz")
updated_metadata <- fread(io$updated.metadata)
colnames(updated_metadata)

# remove overlapping columns in the archR metadata
overlaping.columns <- intersect(colnames(updated_metadata),colnames(archr_metadata))
overlaping.columns <- overlaping.columns[!overlaping.columns%in%c("sample","cell","barcode")]
archr_metadata <- archr_metadata[,which(!colnames(archr_metadata)%in%overlaping.columns),with=F]

###########
## Merge ##
###########

foo <- updated_metadata %>% 
  merge(archr_metadata, by=c("cell","sample","barcode"), all=TRUE)

#############################
## Update ArchR's metadata ##
#############################

bar <- foo %>% 
  .[cell%in%rownames(ArchRProject)] %>% setkey(cell) %>% .[rownames(ArchRProject)] %>%
  as.data.frame() %>% tibble::column_to_rownames("cell")

stopifnot(bar$cell == rownames(getCellColData(ArchRProject)))

for (i in colnames(bar)) {
  ArchRProject <- addCellColData(
    ArchRProject,
    data = bar[[i]], 
    name = i,
    cells = rownames(bar),
    force = TRUE
  )
}

colnames(getCellColData(ArchRProject))

##########
## Save ##
##########

io$metadata.out <- paste0(io$basedir,"/processed/atac/archR/sample_metadata_after_archR.txt.gz")
fwrite(foo, io$metadata.out, sep="\t", na="NA", quote=F)

saveArchRProject(ArchRProject)
