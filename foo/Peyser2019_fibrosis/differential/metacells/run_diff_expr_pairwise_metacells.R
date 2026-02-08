source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
# if (grepl("^(BAI-|SDI-|CI-)", Sys.info()['nodename'])) {
#   io$script <- "/Users/rargelaguet/altos_github/reprogramming/HDF_oskm_scrna_CHAMP/rna/differential/metacells/differential_gene_expr_metacells.R"
# } else {
#   io$script <- "/share/rargelaguet/altos_github/reprogramming/HDF_oskm_scrna_CHAMP/rna/differential/metacells/differential_gene_expr_metacells.R"
# } 
io$script <- "/share/rargelaguet/altos_github/fibrosis/Peyser2019_fibrosis/differential/metacells/diff_gene_expr_metacells.R"
io$metadata <- file.path(io$basedir,"results/metacells/per_condition/metacells_metadata.txt.gz")
io$sce <- file.path(io$basedir,"results/metacells/per_condition/SingleCellExperiment_metacells.rds")


# Options
opts$test_mode <- FALSE
opts$group_variable <- "condition"

io$outdir <- file.path(io$basedir,sprintf("results/differential/metacells/%s",opts$group_variable))
dir.create(io$outdir, showWarnings = F, recursive = T)

########################
## Load cell metadata ##
########################

cell_metadata.dt <- fread(io$metadata)

cell_metadata.dt <- cell_metadata.dt[!is.na(cell_metadata.dt[[opts$group_variable]])]

# print stats
group.stats <- table(cell_metadata.dt[[opts$group_variable]])
print(group.stats)

groups.to.use <- unique(cell_metadata.dt[[opts$group_variable]])

stats.dt <- data.table(groupA=as.character(NA), groupB=as.character(NA), N_A=as.integer(NA), N_B=as.integer(NA))

#########
## Run ##
#########

groupA <- "Bleomycin_day11"; groupB <- "Bleomycin_nintedanib_day11"
# groupA <- "Saline_Control"; groupB <- "Bleomycin_day11"
# groupA <- "Saline_Control"; groupB <- "Bleomycin_nintedanib_day11"

celltypes.to.use <- c("Type_II_alveolar", "dendritic_cells", "monocytes", "Type_I_alveolar", "neutrophils", "mesothelium", "endothelial", "macrophages", "B_cells", "fibroblasts")
# celltypes.to.use <- c("fibroblasts")

# i <- celltypes.to.use[1]
for (i in celltypes.to.use) {
    
    outfile <- file.path(io$outdir,sprintf("%s_%s_vs_%s.txt.gz",i,groupA,groupB))
    if (!file.exists(outfile)) {

      # Define command
      cmd <- sprintf("nohup Rscript %s --metadata %s --sce %s --celltypes %s --groupA %s --groupB %s --group_variable %s --outfile %s &> nohup_%s_vs_%s.out &",
        io$script, io$metadata, io$sce, i, groupA, groupB, opts$group_variable, outfile, groupA, groupB)

      # Run
      print(cmd)
      system(cmd)
      
      # save stats
      # tmp <- cell_metadata.dt[cell_type==i]
      # group.stats <- table(cell_metadata.dt[[opts$group_variable]])
      # stats.dt <- rbind(stats.dt, data.table(groupA=groupA, groupB=groupB, N_A=group.stats[[groupA]], N_B=group.stats[[groupB]]))
      
    } else {
      print(sprintf("%s already exists...",outfile))
    }
}

# Save stats
# fwrite(stats.dt[-1], file.path(io$outdir,"diff_stats.txt"), sep="\t", quote=F)
