source(here::here("settings.R"))
source(here::here("utils.R"))

library(clusterProfiler)
library(org.Hs.eg.db)
library(furrr)

######################
## Define arguments ##
######################

p <- ArgumentParser(description = "GSEA on pseudobulk DE (perturbation vs SCRAMBLE)")
p$add_argument("--cell_type", type = "character", default = "CCL-210", help = "Cell type to run (single)")
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

# I/O
io$diff_results_pseudobulk <- file.path(io$basedir, "results/differential/pseudobulk/perturbation")
io$outdir <- file.path(io$basedir, "results/gene_signatures/gsea_differential/perturbation")
dir.create(io$outdir, showWarnings = F, recursive = T)

# Options (cell_type from argument)
opts$cell_type <- args$cell_type
# opts$pval_threshold <- 0.05   # not used when DE has no p.value; sig from logFC only
opts$min_logFC <- 1

# Foreground = all perturbations except SCRAMBLE; background = SCRAMBLE
opts$background_perturbation <- "SCRAMBLE"
opts$foreground_perturbations <- setdiff(opts$perturbations, opts$background_perturbation)
# opts$foreground_perturbations <- opts$foreground_perturbations[1:3]

# Parallel processing
plan(multisession, workers = 6)
options(future.globals.maxSize = 2 * 1024^3)

##########################################
## Load differential expression results ##
##########################################
# One file per (cell_type, perturbation): diff_expr_{cell_type}_SCRAMBLE_vs_{perturbation}.txt.gz

ct <- opts$cell_type
diff_list <- list()
for (pert in opts$foreground_perturbations) {
  f <- file.path(io$diff_results_pseudobulk, sprintf("diff_expr_%s_SCRAMBLE_vs_%s.txt.gz", ct, pert))
  if (file.exists(f)) {
    diff_list[[pert]] <- fread(f)[, perturbation := pert]
  }
}
if (length(diff_list) == 0) stop(sprintf("No DE files found for cell_type=%s in %s", ct, io$diff_results_pseudobulk))

diff.dt <- rbindlist(diff_list) %>% setnames("gene", "feature")

diff.dt[, sig := abs(logFC) >= opts$min_logFC]

##########################
## Load gene signatures ##
##########################

source(here::here("gene_signatures/load_gene_signatures.R"))
source(here::here("gene_signatures/load_manual_gene_signatures.R"))

tmp <- list()
for (i in names(gene_signatures_list)) {
  tmp[[i]] <- data.table(pathway = i, gene = gene_signatures_list[[i]])
}
gs_manual_t2g.dt <- rbindlist(tmp)

#############################
## Run GSEA using GO terms  ##
#############################
# One contrast per foreground perturbation (vs SCRAMBLE)

# gsea_go_out_list <- setNames(vector("list", length(opts$foreground_perturbations)), opts$foreground_perturbations)

# # pert <- opts$foreground_perturbations[1]
# for (pert in opts$foreground_perturbations) {

#   diff_results_filt.dt <- diff.dt[perturbation == pert]
#   if (nrow(diff_results_filt.dt) == 0) next

#   genes_list <- list()
#   genes_list[["up"]] <- diff_results_filt.dt[logFC > 0 & sig == T, feature]
#   genes_list[["down"]] <- diff_results_filt.dt[logFC < 0 & sig == T, feature]
#   genes_list[["all"]] <- c(genes_list[["down"]], genes_list[["up"]])
#   all.genes <- diff_results_filt.dt$feature

#   gsea_out <- list()
#   for (j in c("down", "up", "all")) {
#     if (length(genes_list[[j]]) > 10) {
#       gsea_out[[j]] <- enrichGO(
#         gene = genes_list[[j]],
#         universe = all.genes,
#         OrgDb = org.Hs.eg.db,
#         keyType = "SYMBOL",
#         ont = "BP",
#         pAdjustMethod = "BH",
#         qvalueCutoff = 0.10,
#         minGSSize = 10,
#         maxGSSize = 500
#       )
#     } else {
#       gsea_out[[j]] <- NULL
#     }
#   }
#   gsea_go_out_list[[pert]] <- gsea_out
# }

# saveRDS(gsea_go_out_list, file.path(io$outdir, sprintf("gsea_go_list_%s_vs_SCRAMBLE.rds", ct)))

##############
## Run GSEA ##
##############

# opts$databases <- c("MSigDB", "MitoCarta", "Reactome", "Manual")
opts$databases <- c("MSigDB", "Reactome", "Manual")

for (db in opts$databases) {

  if (db == "MSigDB") {
    term2gene.dt <- msigdbr_hallmark_t2g.dt
  } else if (db == "MitoCarta") {
    term2gene.dt <- mito_t2g.dt
  } else if (db == "Reactome") {
    term2gene.dt <- reactome_t2g.dt
  } else if (db == "Manual") {
    term2gene.dt <- gs_manual_t2g.dt
  } else {
    stop("Database not recognised")
  }

  gsea_out_list <- setNames(vector("list", length(opts$foreground_perturbations)), opts$foreground_perturbations)

  for (pert in opts$foreground_perturbations) {

    diff_results_filt.dt <- diff.dt[perturbation == pert]
    if (nrow(diff_results_filt.dt) == 0) next

    genes_list <- list()
    genes_list[["up"]] <- diff_results_filt.dt[logFC > 0 & sig == T, feature]
    genes_list[["down"]] <- diff_results_filt.dt[logFC < 0 & sig == T, feature]
    genes_list[["all"]] <- c(genes_list[["down"]], genes_list[["up"]])
    all.genes <- diff_results_filt.dt$feature

    gsea_out <- list()
    for (j in c("down", "up", "all")) {
      if (length(genes_list[[j]]) > 10) {
        gsea_out[[j]] <- enricher(
          gene = genes_list[[j]],
          universe = all.genes,
          TERM2GENE = term2gene.dt,
          pvalueCutoff = 0.01,
          pAdjustMethod = "BH",
          minGSSize = 10,
          maxGSSize = 500,
          qvalueCutoff = 0.1
        )
      } else {
        gsea_out[[j]] <- NULL
      }
    }
    gsea_out_list[[pert]] <- gsea_out
  }

  saveRDS(gsea_out_list, file.path(io$outdir, sprintf("gsea_%s_list_%s_vs_SCRAMBLE.rds", db, ct)))
}
