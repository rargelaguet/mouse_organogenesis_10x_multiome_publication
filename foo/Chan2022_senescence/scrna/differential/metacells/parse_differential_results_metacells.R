# here::i_am("rna/differential/cells/parse_differential_results.R")

# Load default settings
source(here::here("settings.R"))

#####################
## Define settings ##
#####################

# I/O
io$diff_results_dir <- file.path(io$basedir,"results/rna/differential/metacells/alias")
io$outdir <- file.path(io$basedir,"results/rna/differential/metacells/alias/parsed"); dir.create(io$outdir, showWarnings = F, recursive=T)

# Options
opts$groups <- opts$aliases
opts$group_variable <- "alias"

#########################
## Load and parse data ##
#########################

# stats.dt <- data.table(groupA=as.character(NA), groupB=as.character(NA), groupA_N=as.integer(NA), groupB_N=as.integer(NA), included=as.logical(NA))
stats.dt <- data.table(groupA=as.character(NA), groupB=as.character(NA), included=as.logical(NA))
diff_results_list <- list()

# i <- "SRR24317836_GSM7230012_hCiPSCs"; j <- "SRR24317837_GSM7230011_StageIII_Day12"
for (i in opts$groups) {
  for (j in opts$groups) {
    
    if (i!=j) {
      file <- file.path(io$diff_results_dir,sprintf("diff_%s_%s_vs_%s.txt.gz",opts$group_variable,i,j))
      if (file.exists(file)) {
        tmp <- fread(file) %>% .[,c("groupA","groupB"):=list(i,j)]
        
        # Empty file (not enough cells to do DE)
        if (nrow(tmp)>1) {
          # stats.dt <- rbind(stats.dt, data.table(groupA=i, groupB=j, groupA_N=tmp$groupA_N[1], groupB_N=tmp$groupB_N[1], included=TRUE))
          stats.dt <- rbind(stats.dt, data.table(groupA=i, groupB=j, included=TRUE))
          diff_results_list[[sprintf("%s_vs_%s",i,j)]] <- tmp
        } else {
          stats.dt <- rbind(stats.dt, data.table(groupA=i, groupB=j, included=FALSE)) 
        }
      } else {
        print(sprintf("%s not found...",file))
      }
    }
  }
}

##########
## Save ##
##########

fwrite(stats.dt[-1], file.path(io$outdir,"diff_expr_stats.txt.gz"), sep="\t", quote=F, na="NA")
fwrite(rbindlist(diff_results_list), file.path(io$outdir,"diff_expr_results.txt.gz"), sep="\t", quote=F, na="NA")

