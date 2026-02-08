source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
if (grepl("^(BAI-|SDI-|CAM-)", Sys.info()['nodename'])) {
  io$script <- "/Users/rargelaguet/altos_github/reprogramming/Xing2020_bj_oskm/differential/metacells/differential_gene_expr_metacells.R"
} else {
  io$script <- "/share/rargelaguet/altos_github/reprogramming/Xing2020_bj_oskm/differential/metacells/differential_gene_expr_metacells.R"
} 
io$metadata <- file.path(io$basedir,"results/rna/metacells/scTE/metacells_metadata.txt.gz")
io$sce <- file.path(io$basedir,"results/rna/metacells/scTE/SingleCellExperiment_metacells.rds")
io$outdir <- file.path(io$basedir,"results/rna/differential/metacells/alias"); dir.create(io$outdir, showWarnings = F, recursive = T)

# Options
opts$test_mode <- FALSE
opts$group_variable <- "alias"

########################
## Load cell metadata ##
########################

cell_metadata.dt <- fread(io$metadata)

# print stats
group.stats <- table(cell_metadata.dt[[opts$group_variable]])
print(group.stats)

groups.to.use <- unique(cell_metadata.dt[[opts$group_variable]])

#########
## Run ##
#########

if (opts$test_mode) {
  print("Test mode activated, running only a few comparisons...")
  groups.to.use <- groups.to.use %>% head(n=3)
}

stats.dt <- data.table(groupA=as.character(NA), groupB=as.character(NA), N_A=as.integer(NA), N_B=as.integer(NA))

# i <- "c1"; j <- "c2"
for (i in groups.to.use) {
  for (j in groups.to.use) {
    if (i!=j) {
      
      stats.dt <- rbind(stats.dt, data.table(groupA=i, groupB=j, N_A=group.stats[[i]], N_B=group.stats[[j]]))

      outfile <- sprintf("%s/diff_%s_%s_vs_%s.txt.gz", io$outdir,opts$group_variable,i,j)
      if (!file.exists(outfile)) {

        # Define command
        cmd <- sprintf("nohup Rscript %s --metadata %s --sce %s --groupA %s --groupB %s --group_variable %s --outfile %s &> nohup_%s_vs_%s.out &", 
          io$script, io$metadata, io$sce, i, j, opts$group_variable, outfile, i, j)
        # cmd <- sprintf("sbatch -n 4 --mem 30G -o slurm_%s_vs_%s.out --wrap 'Rscript %s --metadata %s --sce %s --groupA %s --groupB %s --group_variable %s --outfile %s'", 
        #   i, j, io$script, io$metadata, io$sce, i, j, opts$group_variable, outfile)

        # Run
        print(cmd)
        system(cmd)

        
      } else {
        print(sprintf("%s already exists...",outfile))
      }
    }
  }
}

# Save stats
fwrite(stats.dt[-1], file.path(io$outdir,"diff_stats.txt"), sep="\t", quote=F)
