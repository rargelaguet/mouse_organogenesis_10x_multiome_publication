here::i_am("atac/archR/processing/2_create_archR_metadata.R")

source(here::here("settings.R"))

suppressPackageStartupMessages(library(ArchR))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--archr_directory',    type="character",    help='ArchR directory')
p$add_argument('--metadata',    type="character",    help='metadata file')
p$add_argument('--outfile',     type="character",    help='Output file')

args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
# args$metadata <- file.path(io$basedir,"results/rna/mapping/sample_metadata_after_mapping.txt.gz")
# args$outfile <- file.path(io$basedir,"processed/atac/archR/sample_metadata_after_archR.txt.gz")
## END TEST ##

###################
## Load metadata ##
###################

sample_metadata <- fread(args$metadata)

########################
## Load ArchR project ##
########################

# source(here::here("atac/archR/load_archR_project.R"))

setwd(args$archr_directory)

addArchRGenome("mm10")
addArchRThreads(threads = args$threads)

ArchRProject <- loadArchRProject(args$archr_directory)[sample_metadata$cell]

######################
## Load ArchR stats ##
######################
  
# fetch archR's metadata
archR_metadata <- getCellColData(ArchRProject) %>%
  as.data.table(keep.rownames = T) %>% setnames("rn","cell") %>%
  .[,c("cell", "TSSEnrichment", "ReadsInTSS", "PromoterRatio", "NucleosomeRatio", "nFrags",  "BlacklistRatio")]

cols.to.rename <- c("TSSEnrichment","ReadsInTSS","PromoterRatio","NucleosomeRatio","nFrags","BlacklistRatio")
idx.cols.to.rename <- which(colnames(archR_metadata)%in%cols.to.rename)
colnames(archR_metadata)[idx.cols.to.rename] <- paste0(colnames(archR_metadata)[idx.cols.to.rename], "_atac")

###########
## Merge ##
###########

# print stats
print(sprintf("%s cells have both RNA expression and chromatin accessibility measurements",length(intersect(sample_metadata$cell,archR_metadata$cell))))
print(sprintf("%s cells have RNA expression, but do not have chromatin accessibility measurements",sum(!sample_metadata$cell%in%archR_metadata$cell)))
print(sprintf("%s cells have chromatin accessibility, but do not have RNA expression measurements",sum(!archR_metadata$cell%in%sample_metadata$cell)))

# merge
sample_metadata_tosave <- sample_metadata %>% 
  merge(archR_metadata,by="cell", all=TRUE) 

# Fill missing entries for cells that did not pass QC for the RNA
sample_metadata_tosave %>%
  .[is.na(sample),sample:=strsplit(cell,"#") %>% map_chr(1)] %>%
  .[is.na(stage),stage:=strsplit(sample,"_") %>% map_chr(1)] %>%
  .[is.na(barcode),barcode:=strsplit(cell,"#") %>% map_chr(2)]

# round
sample_metadata_tosave[,c("TSSEnrichment_atac","NucleosomeRatio_atac","PromoterRatio_atac","BlacklistRatio_atac"):=list(round(TSSEnrichment_atac,2),round(NucleosomeRatio_atac,2),round(PromoterRatio_atac,2),round(BlacklistRatio_atac,2))]
sample_metadata_tosave[,c("ribosomal_percent_RNA","mitochondrial_percent_RNA"):=list(round(ribosomal_percent_RNA,2),round(mitochondrial_percent_RNA,2))]

# sanity checks
stopifnot(all(!is.na(sample_metadata_tosave$sample)))
stopifnot(all(!is.na(sample_metadata_tosave$stage)))
stopifnot(all(!is.na(sample_metadata_tosave$barcode)))

# print stats
table(sample_metadata$sample)
table(sample_metadata$stage)

#############################
## Update ArchR's metadata ##
#############################

metadata.to.archR <- sample_metadata_tosave %>% 
  .[cell%in%rownames(ArchRProject)] %>% setkey(cell) %>% .[rownames(ArchRProject)] %>%
  as.data.frame() %>% tibble::column_to_rownames("cell")

# stopifnot(all(metadata.to.archR$TSSEnrichment_atac == getCellColData(ArchRProject,"TSSEnrichment")[[1]]))

for (i in colnames(metadata.to.archR)) {
  ArchRProject <- addCellColData(
    ArchRProject,
    data = metadata.to.archR[[i]], 
    name = i,
    cells = rownames(metadata.to.archR),
    force = TRUE
  )
}

head(getCellColData(ArchRProject))

##########
## Save ##
##########

fwrite(sample_metadata_tosave, args$outfile, sep="\t", na="NA", quote=F)
saveArchRProject(ArchRProject, drop=TRUE)
