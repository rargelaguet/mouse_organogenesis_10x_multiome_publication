source(here::here("settings.R"))
source(here::here("utils.R"))

library(GenomicRanges)

#####################
## Define settings ##
#####################

io$feature_metadata <- "/Users/rargelaguet/data/BJ_fibroblasts_reprogramming_10x_multiome_CHAMP/processed/atac/archR/PeakCalls/peak_metadata.tsv.gz"
io$deeptols.output <- file.path(io$basedir,"results/atac/deeptools/multiBigwigSummary_BJ_reprogramming.tab")
io$outfile <- file.path(io$basedir,"results/atac/deeptools/SingleCellExperiment_BJ_reprogramming.rds")

# opts$groups <- c(
#   "p8_Control", 
#   "p11_Control", 
#   "p14_Control", 
#   "p19_Control", 
#   "p22_Control", 
#   "p24_Control", 
#   "p29_Control", 
#   "p8_TG", 
#   "p11_TG", 
#   "p14_TG", 
#   "p19_TG", 
#   "p22_TG", 
#   "p24_TG", 
#   "p29_TG"
# )

##########################
## Load sample metadata ##
##########################

# sample_metadata.dt <- fread(io$sample_metadata)

###############
## Load data ##
###############

atac.dt <- fread(io$deeptols.output, skip = 1)

colnames(atac.dt) <- c("chr", "start", "end", "Control_p11", "Control_p14", "Control_p19", "Control_p22", "Control_p24", "Control_p29", "Control_p8", "TG_p11", "TG_p14", "TG_p19", "TG_p22", "TG_p24", "TG_p29", "TG_p8")
atac.mtx <- copy(atac.dt) %>%
  .[,feature:=sprintf("%s:%s-%s",chr,start,end)] %>% 
  setcolorder(c('feature',setdiff(names(dt),'feature'))) %>%
  .[,c("chr","start","end"):=NULL] %>%
  matrix.please# %>% .[,opts$samples]

duplicated.sites <- rownames(atac.mtx)[duplicated(rownames(atac.mtx))]

atac.mtx <- atac.mtx[!rownames(atac.mtx)%in%duplicated.sites,]

############################
## Create sample metadata ##
############################

#############################
## Create feature metadata ##
#############################

feature_metadata.dt <- fread(io$feature_metadata, select=c(1,2,3,5,13)) %>%
  setnames(c("chr","start","end","score","peak_type")) %>%
  .[,feature:=sprintf("%s:%s-%s",gsub("chr","",chr),start,end)] %>%
  .[,c("feature","chr","start","end","score","peak_type")] %>%
  .[!feature%in%duplicated.sites] %>% setkey(feature) %>% .[rownames(atac.mtx)]

#################################
## Create SingleCellExperiment ##
#################################

# coldata_to_sce <- DataFrame(sample_metadata.dt[sample%in%samples.to.use] %>% setkey(sample) %>% .[colnames(counts.mtx)] %>% tibble::column_to_rownames("sample"))
# stopifnot(colnames(counts.mtx)==rownames(coldata_to_sce))

# gr <- GRanges(seqnames=feature_metadata.dt$chr, ranges=IRanges(start=feature_metadata.dt$start, end=feature_metadata.dt$end))
# gr <- GRanges(seqnames=feature_metadata.dt$chr, ranges=IRanges(start=feature_metadata.dt$start, end=feature_metadata.dt$end))
gr <- makeGRangesFromDataFrame(feature_metadata.dt, keep.extra.columns = T)

# rowdata_to_sce <- DataFrame(feature_metadata.dt[,c("feature","annotation")] %>% tibble::column_to_rownames("feature"))
# stopifnot(rownames(atac.mtx)==rownames(rowdata_to_sce))

sce <- SingleCellExperiment(
  assays = list("atac"=atac.mtx), 
  # rowData = rowdata_to_sce,
  rowRanges = gr
  # colData = coldata_to_sce
)

##########
## Save ##
##########

saveRDS(sce, io$outfile)
