#####################
## Define settings ##
#####################

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

# I/O
io$archR.diff.dir <- file.path(io$basedir,"results_new/atac/archR/differential/GeneScoreMatrix_TSS")
io$outdir <- file.path(io$basedir,"results_new/atac/archR/differential/GeneScoreMatrix_TSS/markers"); dir.create(io$outdir, showWarnings = F)

# Options
# opts$groups <- strsplit(list.files(io$diff.dir, pattern="*.gz"),"_vs_") %>% map(~ .[[1]]) %>% unlist %>% unique
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
)# %>% head(n=4)

opts$min.MeanDiff <- 0.10
opts$fdr <- 0.01

# Minimum fraction of significant differential pairwise comparisons
opts$score <- 0.75

##################
## Load results ##
##################

dt <- opts$celltypes %>% map(function(i) { opts$celltypes %>% map(function(j) {
  file <- sprintf("%s/GeneScoreMatrix_TSS_%s_vs_%s.txt.gz", io$archR.diff.dir,i,j)
  if (file.exists(file)) {
    fread(file) %>% .[,c("celltypeA","celltypeB"):=list(as.factor(i),as.factor(j))] %>%
      return
  }
}) %>% rbindlist }) %>% rbindlist %>% 
  .[,sig:=FALSE] %>% .[abs(MeanDiff)>=opts$min.MeanDiff & FDR<=opts$fdr,sig:=TRUE] %>%
  .[,direction:=c("up","down")[as.numeric(MeanDiff<0)+1]]  # up = higher accessibility in celltype A

ncelltypes <- length(intersect(unique(dt$celltypeA),unique(dt$celltypeB)))

#########################
## Define marker genes ##
#########################

foo <- dt[,.(score=sum(sig==T & direction=="up")), by=c("celltypeA","name")] %>% setnames("celltypeA","celltype")
bar <- dt[,.(score=sum(sig==T & direction=="down")), by=c("celltypeB","name")] %>% setnames("celltypeB","celltype")
  
markers_peaks.dt <- merge(foo,bar,by=c("celltype","name"), all=TRUE) %>% .[,score:=score.x+score.y] %>%
  .[,c("score.x","score.y"):=NULL] %>%
  .[,score:=round(score/(ncelltypes+1),2)] %>%
  .[score>=opts$score] %>%
  setorder(celltype,-score)
rm(foo,bar)

length(unique(markers_peaks.dt$name))

##########
## Save ##
##########

fwrite(markers_peaks.dt, file.path(io$outdir,"marker_genes.txt.gz"))

