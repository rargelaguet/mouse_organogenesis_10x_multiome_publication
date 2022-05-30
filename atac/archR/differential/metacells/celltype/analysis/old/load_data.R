# i <- opts$celltypes[2]; j <- opts$celltypes[1]

diff.dt <- opts$celltypes %>% map(function(i) { opts$celltypes %>% map(function(j) {
  file <- sprintf("%s/%s_%s_vs_%s.txt.gz", io$diff.dir,opts$matrix,i,j)
  if (file.exists(file)) {
    fread(file) %>% .[,c("celltypeA","celltypeB"):=list(as.factor(i),as.factor(j))] %>% return
  }
}) %>% rbindlist }) %>% rbindlist %>% 
  .[,MeanDiff:=-MeanDiff] %>% # change sign to keep the groupB - groupA consistency
  .[,sig:=FALSE] %>% .[abs(MeanDiff)>=opts$min.MeanDiff & FDR<=opts$fdr,sig:=TRUE] %>%
  .[,direction:=c("up","down")[as.numeric(MeanDiff>0)+1]]  # up = higher accessibility in celltype A

# ad hoc
# if ("name"%in%colnames(atac_diff_cells.dt)) {
#   atac_diff_cells.dt[,idx:=NULL] %>% setnames("name","idx")
# }
