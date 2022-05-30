here::i_am("atac/archR/differential/analysis/PeakMatrix/define_marker_peaks.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#############
## Options ##
#############

opts$matrix <- "PeakMatrix"
opts$group_variable <- "celltype.mapped"
opts$min.MeanDiff <- 0.10
opts$fdr <- 0.01
opts$score <- 0.75 # Minimum fraction of significant differential pairwise comparisons

#########
## I/O ##
#########

io$diff.dir <- file.path(io$basedir,sprintf("results/atac/archR/differential/%s/%s",opts$group_variable,opts$matrix))
io$outdir <- file.path(io$basedir,sprintf("results/atac/archR/differential/%s/%s/markers",opts$group_variable,opts$matrix)); dir.create(io$outdir, showWarnings = F)

##################
## Load results ##
##################

source(here::here("atac/archR/differential/analysis/load_data.R"))

###################
## Sanity checks ##
###################

# Load stats
diff_stats.dt <- fread(file.path(io$diff.dir,"diff_stats.txt")) %>% setnames(c("celltypeA","celltypeB","N_groupA","N_groupB"))

# check if some DA comparison is missing
tmp <- diff_stats.dt %>% 
  merge(diff.dt[,c("celltypeA","celltypeB")] %>% unique %>% .[,done:=TRUE], all.x=TRUE, by=c("celltypeA","celltypeB")) %>%
  .[is.na(done),done:=FALSE]
stopifnot(tmp$done==TRUE)

#########################
## Define marker genes ##
#########################

ncelltypes <- unique(c(as.character(unique(diff.dt$celltypeA)),as.character(unique(diff.dt$celltypeB)))) %>% length

foo <- diff.dt[,.(score=sum(sig==T & direction=="up")), by=c("celltypeA","idx")] %>% setnames("celltypeA","celltype")
bar <- diff.dt[,.(score=sum(sig==T & direction=="down")), by=c("celltypeB","idx")] %>% setnames("celltypeB","celltype")
  
markers_peaks.dt <- merge(foo,bar,by=c("celltype","idx"), all=TRUE) %>% .[,score:=score.x+score.y] %>%
  .[,c("score.x","score.y"):=NULL] %>%
  # .[,score:=round(score/(ncelltypes+1),2)] %>%
  .[,score:=round(score/(ncelltypes-1),2)] %>%
  setorder(celltype,-score)
# rm(foo,bar)

stopifnot(max(markers_peaks.dt$score,na.rm=T)==1)


##############################################
## Add MeanDiff values from pseudobulk data ##
##############################################

diff_pseudobulk.dt <- file.path(io$basedir,"results/atac/archR/differential/pseudobulk/celltype.mapped/PeakMatrix/differential_atac_PeakMatrix_pseudobulk_summary.txt.gz") %>% fread
markers_peaks.dt <- markers_peaks.dt %>% merge(diff_pseudobulk.dt, by=c("celltype","idx"))

##########
## Save ##
##########

# Save marker score for all combination of genes and cell types
length(unique(markers_peaks.dt$idx))
length(unique(markers_peaks.dt$celltype))
fwrite(markers_peaks.dt, file.path(io$outdir,"marker_peaks_upregulated_all.txt.gz"), sep="\t")

# Save marker score for strong markers
markers_peaks_filt.dt <- markers_peaks.dt %>% .[score>=opts$score & diff>=opts$min.MeanDiff]
length(unique(markers_peaks_filt.dt$idx))
length(unique(markers_peaks_filt.dt$celltype))
fwrite(markers_peaks_filt.dt, file.path(io$outdir,"marker_peaks_upregulated_filtered.txt.gz"), sep="\t")

