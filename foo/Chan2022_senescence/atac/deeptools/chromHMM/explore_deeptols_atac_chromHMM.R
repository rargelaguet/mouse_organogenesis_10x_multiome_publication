# here::i_am("dimensionality_reduction/dimensionality_reduction.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(scran))

#####################
## Define settings ##
#####################

# I/O
io$metadata <- file.path(io$basedir,"processed/atac/sample_metadata.txt")
io$feature_metadata <- "/Users/rargelaguet/data/hg38_regulation/chromHMM/HDF_6_dense.merged.bed"
io$sce <- file.path(io$basedir,"processed/atac/SingleCellExperiment.rds")
io$outdir <- file.path(io$basedir,"results/atac/dimensionality_reduction"); dir.create(io$outdir, showWarnings = F)



# Options
opts$groups <- c(
  "p8_Control", 
  "p11_Control", 
  "p14_Control", 
  "p19_Control", 
  "p22_Control", 
  "p24_Control", 
  "p29_Control", 
  "p8_TG", 
  "p11_TG", 
  "p14_TG", 
  "p19_TG", 
  "p22_TG", 
  "p24_TG", 
  "p29_TG"
)

##########################
## Load sample metadata ##
##########################

sample_metadata.dt <- fread(io$metadata) %>%
  .[Group%in%opts$groups]

table(sample_metadata.dt$Group)

###############
## Load data ##
###############

sce <- readRDS(io$sce)[,sample_metadata.dt$sample]

# Add sample metadata as colData
colData(sce) <- sample_metadata.dt %>% tibble::column_to_rownames("sample") %>% DataFrame



ggplot(to.plot, aes(x=population, y=value, fill=population)) +
  geom_boxplot(outlier.shape=NA) +
  facet_wrap(~cell_line, scales = "fixed", nrow=1) +
  # scale_fill_manual(values=opts$celltype.colors) +
  coord_cartesian(ylim=c(0,0.005)) +
  labs(x="", y="ATAC signal") +
  theme_classic() +
  theme(
    axis.text.y = element_text(size=rel(0.90), color="black"),
    axis.text.x = element_text(colour="black",size=rel(1.0), angle=90, hjust=1, vjust=0.5),
    # axis.title = element_text(colour="black",size=rel(0.75)),
    # axis.ticks.x = element_blank(),
    legend.position = "right"
  )


atac.dt <- fread("/Users/rargelaguet/data/adult_fibroblast_reprogramming_10x_multiome_reik/results/atac/archR/deeptools_pseudobulk/multiBigwigSummary_chromHMM.tab", skip = 1)
colnames(atac.dt) <- c("chr", "start", "end", "25F",  "25F_day0_GAPDH_high",  "25F_day0_GAPDH_low",  "28F",  "28F_day0_GAPDH_high",  "28F_day0_GAPDH_low",  "52F",  "52F_day0_GAPDH_high",  "52F_day0_GAPDH_low",  "54F",  "54F_day0_GAPDH_high",  "54F_day0_GAPDH_low",  "63F",  "63F_day0_GAPDH_high",  "63F_day0_GAPDH_low",  "63F_day13_Maturation",  "63F_day13_Refractory",  "63F_day27_Unsorted_Reversion",  "63F_day2_Initiation",  "63F_day2_Refractory",  "63F_day7_Initiation",  "63F_day7_Maturation",  "63F_day7_Refractory",  "63F_day7_re6_unknown",  "70M",  "70M_day0_GAPDH_high",  "70M_day0_GAPDH_low",  "BJ_day0_NoReprogramming",  "BJ_day13_NoReprogramming",  "BJ_day27_Unsorted_Reversion",  "BJ_day2_NoReprogramming",  "BJ_day2_Reprogramming",  "BJ_day7_NoReprogramming",  "FSFA",  "FSFA_day0_GAPDH_high",  "FSFB",  "FSFB_day0_GAPDH_high",  "FSFB_day0_GAPDH_low",  "FSF_day13_Refractory",  "FSF_day27_FTR_Reversion",  "FSF_day27_TR_Reversion",  "FSF_day27_Unsorted_Reversion",  "FSF_day2_Initiation",  "FSF_day2_Refractory",  "FSF_day7_Initiation",  "FSF_day7_Maturation",  "FSF_day7_Refractory",  "FSF_day7_re6_unknown")
atac.dt[,idx:=sprintf("%s:%s-%s",chr,start,end)]

samples.to.plot <- c("25F", "28F", "52F", "54F", "63F", "70M", "FSFA", "FSFB")

chromHMM.dt <- fread("/Users/rargelaguet/data/hg38_regulation/chromHMM/HDF_6_dense.merged.bed", select=c(1,2,3,4)) %>%
  setnames(c("chr","start","end","annotation"))
chromHMM.dt[,idx:=sprintf("chr%s:%s-%s",chr,start,end)]
chromHMM.dt[,c("chr","start","end"):=NULL]
atac.dt <- atac.dt %>% merge(chromHMM.dt, by="idx")

to.plot <- atac.dt %>%
  .[,c("chr","start","end","annotation",samples.to.plot), with = FALSE] %>% 
  melt(id.vars=c("chr","start","end","annotation"), variable.name="sample")

tmp <- atac.dt[,c("chr","start","end","25F","70M")] %>%
  .[,diff:=`70M`-`25F`]

