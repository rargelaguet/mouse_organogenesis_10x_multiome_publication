source(here::here("settings.R"))
source(here::here("utils.R"))

library(enrichR)

######################
## Define arguments ##
######################

p <- ArgumentParser(description = "GSEA (enrichR) on pseudobulk DE (perturbation vs SCRAMBLE)")
p$add_argument("--cell_type", type = "character", default = "CCL-210", help = "Cell type to run (single)")
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

# I/O
io$diff_results_pseudobulk <- file.path(io$basedir, "results/differential/pseudobulk/perturbation")
io$outdir <- file.path(io$basedir, "results/gene_signatures/gsea_differential")
dir.create(io$outdir, showWarnings = F, recursive = T)

# Options
opts$cell_type <- args$cell_type
opts$min_logFC <- 1

# Foreground = all perturbations except SCRAMBLE; background = SCRAMBLE
opts$background_perturbation <- "SCRAMBLE"
opts$foreground_perturbations <- setdiff(opts$perturbations, opts$background_perturbation)
# opts$foreground_perturbations <- opts$foreground_perturbations[1:3]

opts$databases <- c("MSigDB_Hallmark_2020")
# opts$databases <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023", "Reactome_2022", "MSigDB_Hallmark_2020")

##########################################
## Load differential expression results ##
##########################################
# One file per (cell_type, perturbation): diff_expr_{cell_type}_SCRAMBLE_vs_{perturbation}.txt.gz

diff_list <- list()
for (pert in opts$foreground_perturbations) {
  f <- file.path(io$diff_results_pseudobulk, sprintf("diff_expr_%s_SCRAMBLE_vs_%s.txt.gz", opts$cell_type, pert))
  if (file.exists(f)) {
    diff_list[[pert]] <- fread(f)[, perturbation := pert]
  }
}
if (length(diff_list) == 0) stop(sprintf("No DE files found for cell_type=%s in %s", opts$cell_type, io$diff_results_pseudobulk))

diff.dt <- rbindlist(diff_list) %>% setnames("gene", "feature")

diff.dt[, sig := abs(logFC) >= opts$min_logFC]

########################
## GSEA using EnrichR ##
########################
# One contrast per foreground perturbation (vs SCRAMBLE)

gsea_enrichr_list <- list()

for (db in opts$databases) {
  print(db)
  gsea_enrichr_list[[db]] <- list()

  for (pert in opts$foreground_perturbations) {

    tmp <- diff.dt[perturbation == pert]
    if (nrow(tmp) == 0) next

    genes.up <- tmp[logFC > 0 & sig == TRUE, feature]
    genes.down <- tmp[logFC < 0 & sig == TRUE, feature]

    contrast_name <- sprintf("%s_vs_SCRAMBLE", pert)
    gsea_enrichr_list[[db]][[contrast_name]] <- list()

    if (length(genes.up) > 5) {
      gsea_enrichr_list[[db]][[contrast_name]][["up"]] <- enrichr(
        genes = genes.up,
        database = db
      )[[1]] %>%
        as.data.table() %>%
        .[, .(Term, Overlap, Adjusted.P.value, Odds.Ratio, Genes)] %>%
        setnames(c("gene_signature", "overlap", "pvalue", "odds_ratio", "genes"))
    } else {
      gsea_enrichr_list[[db]][[contrast_name]][["up"]] <- NULL
    }

    if (length(genes.down) > 5) {
      gsea_enrichr_list[[db]][[contrast_name]][["down"]] <- enrichr(
        genes = genes.down,
        database = db
      )[[1]] %>%
        as.data.table() %>%
        .[, .(Term, Overlap, Adjusted.P.value, Odds.Ratio, Genes)] %>%
        setnames(c("gene_signature", "overlap", "pvalue", "odds_ratio", "genes"))
    } else {
      gsea_enrichr_list[[db]][[contrast_name]][["down"]] <- NULL
    }
  }
}

saveRDS(gsea_enrichr_list, file.path(io$outdir, sprintf("gsea_enrichr_list_%s_vs_SCRAMBLE.rds", opts$cell_type)))


