
#############################################
## Load results from differential analysis ##
#############################################

diff.dt <- opts$celltypes %>% map(function(j) {
  file <- sprintf("%s/%s_WT_vs_KO.txt.gz", io$indir,j)
  if (file.exists(file)) {
    fread(file, select=c(1,2,3)) %>% .[,c("celltype"):=list(j)]
  }
}) %>% rbindlist %>%
  # .[,MeanDiff:=-MeanDiff] # change sign to keep the groupB - groupA consistency
  .[, sign := ifelse(MeanDiff>0,"Upregulated in KO","Downregulated in KO")] %>%
  .[, sig := (FDR<=opts$threshold_fdr & abs(MeanDiff)>=opts$min.MeanDiff)]


# Print stats
print(sprintf("Number of celltypes: %s",length(unique(diff.dt$celltype))))
print(sprintf("Number of features: %s",length(unique(diff.dt$idx))))

