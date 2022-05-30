here::i_am("atac/archR/differential/metacells/celltype/analysis/PeakMatrix/define_marker_peaks.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# Options
opts$matrix <- "PeakMatrix"
opts$group_variable <- "celltype"
opts$min.diff <- 0.50
opts$fdr <- 0.10
opts$score <- 0.75 # Minimum fraction of significant differential pairwise comparisons

# I/O
io$basedir <- file.path(io$basedir,"test")
io$diff.dir <- file.path(io$basedir,sprintf("results/atac/archR/differential/metacells/%s/%s",opts$group_variable,opts$matrix))
io$diff.stats <- file.path(io$diff.dir,"diff_stats.txt")
io$outdir <- file.path(io$basedir,sprintf("results/atac/archR/differential/metacells/%s/%s/markers",opts$group_variable,opts$matrix)); dir.create(io$outdir, showWarnings = F)

##################
## Load results ##
##################

# source(here::here("atac/archR/differential/analysis/load_data.R"))

# i <- opts$celltypes[1]; j <- opts$celltypes[2]
diff.dt <- opts$celltypes %>% map(function(i) { opts$celltypes %>% map(function(j) {
  file <- sprintf("%s/%s_%s_vs_%s.txt.gz", io$diff.dir,opts$matrix,i,j)
  if (file.exists(file)) {
    fread(file, select=c(1,2,3)) %>% 
      .[,c("celltypeA","celltypeB"):=list(as.factor(i),as.factor(j))] %>% 
      setnames("logFC","diff") %>%
      return
  }
}) %>% rbindlist }) %>% rbindlist %>% 
  .[is.na(diff),c("diff","padj_fdr"):=list(0,1)] %>%
  # .[,diff:=-diff] %>% # change sign to keep the groupB - groupA consistency
  .[,sig:=ifelse(abs(diff)>=opts$min.diff & padj_fdr<=opts$fdr,TRUE,FALSE)] %>%
  .[,direction:=c("up","down")[as.numeric(diff>0)+1]]  # up = higher accessibility in celltype A

# diff.dt[is.na(sig),sig:=FALSE]

###################
## Sanity checks ##
###################

sum(is.na(diff.dt$sig))

# Load stats
diff_stats.dt <- fread(io$diff.stats) %>% 
  setnames(c("celltypeA","celltypeB","N_groupA","N_groupB"))

# check if some DA comparison is missing
tmp <- diff_stats.dt %>% 
  merge(diff.dt[,c("celltypeA","celltypeB")] %>% unique %>% .[,done:=TRUE], all.x=TRUE, by=c("celltypeA","celltypeB")) %>%
  .[is.na(done),done:=FALSE]
stopifnot(tmp$done==TRUE)

#########################
## Define marker genes ##
#########################

ncelltypes <- unique(c(as.character(unique(diff.dt$celltypeA)),as.character(unique(diff.dt$celltypeB)))) %>% length

foo <- diff.dt[,.(N=.N, diff=mean(diff), score=sum(sig==T & direction=="up")), by=c("celltypeA","feature")] %>% setnames("celltypeA","celltype")
bar <- diff.dt[,.(N=.N, diff=mean(diff), score=sum(sig==T & direction=="down")), by=c("celltypeB","feature")] %>% setnames("celltypeB","celltype")
  
markers_peaks.dt <- merge(foo[,c("celltype","feature","score")],bar[,c("celltype","feature","score")],by=c("celltype","feature"), all=TRUE) %>% 
  .[is.na(score.x),score.x:=0] %>% .[is.na(score.y),score.y:=0] %>%
  .[,score:=score.x+score.y] %>%
  .[,c("score.x","score.y"):=NULL] %>%
  # .[,score:=round(score/(ncelltypes+1),2)] %>%
  .[,score:=round(score/(ncelltypes-1),2)] %>%
  setorder(celltype,-score)
# rm(foo,bar)

stopifnot(max(markers_peaks.dt$score,na.rm=T)==1)

# diff.dt[feature=="chr11:11901425-11902025" & celltypeA=="Cardiomyocytes"]
# diff.dt[feature=="chr11:11901425-11902025" & celltypeB=="Cardiomyocytes"] %>% View
# diff_stats.dt[celltypeA=="Blood_progenitors_1"]

##############################################
## Add MeanDiff values from pseudobulk data ##
##############################################

tmp <- merge(foo[,c("celltype","feature","N","diff")], bar[,c("celltype","feature","N","diff")], by=c("celltype","feature"), all=TRUE) %>% 
  .[is.na(N.x),N.x:=0] %>% .[is.na(N.y),N.y:=0] %>%
  .[,diff:=weighted.mean(c(-diff.x,diff.y), c(N.x,N.y)),by=c("celltype","feature")]

markers_peaks.dt <- markers_peaks.dt %>% 
  merge(tmp[,c("celltype","feature","diff")], by=c("celltype","feature"))

##########
## Save ##
##########

# Save marker score for all combination of genes and cell types
length(unique(markers_peaks.dt$feature))
length(unique(markers_peaks.dt$celltype))
fwrite(markers_peaks.dt, file.path(io$outdir,sprintf("%s_markers_upregulated_all.txt.gz",opts$matrix)), sep="\t")

# Save marker score for strong markers
markers_peaks_filt.dt <- markers_peaks.dt %>% .[score>=0.90 & diff>=2]
length(unique(markers_peaks_filt.dt$feature))
length(unique(markers_peaks_filt.dt$celltype))
table(markers_peaks_filt.dt$celltype)
fwrite(markers_peaks_filt.dt, file.path(io$outdir,sprintf("%s_markers_upregulated_filt.txt.gz",opts$matrix)), sep="\t")

##########
## TEST ##
##########

# foo[celltype=="Mesenchyme"]
# bar[celltype=="Mesenchyme"]
# markers_peaks.dt[celltype=="Mesenchyme"]
# markers_peaks_filt.dt[celltype=="Epiblast"] %>% View
