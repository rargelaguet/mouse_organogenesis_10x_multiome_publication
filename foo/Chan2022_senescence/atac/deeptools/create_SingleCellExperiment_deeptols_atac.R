source(here::here("settings.R"))
source(here::here("utils.R"))

library(GenomicRanges)

#####################
## Define settings ##
#####################

io$feature_metadata <- "/Users/rargelaguet/data/hg38_regulation/chromHMM/HDF_6_dense.merged.bed"
io$deeptols.output <- file.path(io$basedir,"results/atac/deeptools/multiBigwigSummary_chromHMM.tab")
io$outfile <- file.path(io$basedir,"results/atac/deeptools/SingleCellExperiment_chromHMM.rds")

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
  .[,idx:=sprintf("%s:%s-%s",chr,start,end)] %>% 
  setcolorder(c('idx',setdiff(names(dt),'idx'))) %>%
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

feature_metadata.dt <- fread(io$feature_metadata, select=c(1,2,3,4)) %>%
  setnames(c("chr","start","end","annotation")) %>%
  .[,idx:=sprintf("%s:%s-%s",chr,start,end)] %>%
  .[,c("idx","annotation","chr","start","end")] %>%
  .[!idx%in%duplicated.sites] %>% setkey(idx) %>% .[rownames(atac.mtx)]

#################################
## Create SingleCellExperiment ##
#################################

# coldata_to_sce <- DataFrame(sample_metadata.dt[sample%in%samples.to.use] %>% setkey(sample) %>% .[colnames(counts.mtx)] %>% tibble::column_to_rownames("sample"))
# stopifnot(colnames(counts.mtx)==rownames(coldata_to_sce))

# gr <- GRanges(seqnames=feature_metadata.dt$chr, ranges=IRanges(start=feature_metadata.dt$start, end=feature_metadata.dt$end))
# gr <- GRanges(seqnames=feature_metadata.dt$chr, ranges=IRanges(start=feature_metadata.dt$start, end=feature_metadata.dt$end))
gr <- makeGRangesFromDataFrame(feature_metadata.dt, keep.extra.columns = T)

# rowdata_to_sce <- DataFrame(feature_metadata.dt[,c("idx","annotation")] %>% tibble::column_to_rownames("idx"))
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

