source(here::here("settings.R"))
source(here::here("utils.R"))

library(pheatmap)

#####################
## Define settings ##
#####################

# Options
opts$cell_types <- c("CCL-210", "BEAS-2B", "A549", "BEAS-2B_co-culture", "CCL-210_co-culture")
opts$perturbations <- setdiff(opts$perturbations, "SCRAMBLE")
# opts$perturbations <- c("nintedanib", "TGFBR1", "TGFBR2", "TGFBR1_TRADD","CD47", "CRLF1", "CD248", "TNS1", "RUNX1", "STAT3", "APOE", "STAT6")
opts$log_qval_threshold <- 3
# opts$databases <- c("MSigDB", "Manual")
opts$databases <- c("MSigDB")

# I/O
io$gsea_results_dir <- file.path(io$basedir, "results/gene_signatures/gsea_differential/perturbation")
io$outdir <- file.path(io$gsea_results_dir, "pdf")
dir.create(io$outdir, showWarnings = F, recursive = T)
dir.create(file.path(io$outdir, "all_celltypes"), showWarnings = F, recursive = T)

#############################################
## Load GSEA results across all cell types ##
#############################################
# Perturbation GSEA: one file per cell type, gsea_{database}_list_{cell_type}_vs_SCRAMBLE.rds
# Each file: list keyed by perturbation, each with up/down/all

extract_pert_from_res <- function(res, pert) {
  if (!pert %in% names(res)) return(NULL)
  out <- list()
  for (direction in c("up", "down")) {
    tmp <- res[[pert]][[direction]]
    if (is.null(tmp)) next
    dt <- as.data.table(tmp@result)
    if (nrow(dt) == 0) next
    dt[, log_qval := -log10(qvalue)]
    dt[, Description := gsub("HALLMARK_", "", Description)]
    # Up = positive, down = negative (blue to red scale)
    dt[, signed_log_qval := (if (direction == "up") 1 else -1) * log_qval]
    out[[direction]] <- dt[, .(gene_set = Description, log_qval, direction, signed_log_qval)]
  }
  if (length(out) == 0) return(NULL)
  rbindlist(out)
}

for (db in opts$databases) {

  all_dt <- list()
  for (ct in opts$cell_types) {
    f <- file.path(io$gsea_results_dir, sprintf("gsea_%s_list_%s_vs_SCRAMBLE.rds", db, ct))
    if (!file.exists(f)) next
    res <- readRDS(f)
    for (pert in opts$perturbations) {
      dt <- extract_pert_from_res(res, pert)
      if (!is.null(dt)) {
        dt[, c("cell_type", "perturbation") := .(ct, pert)]
        all_dt[[paste(ct, pert, sep = "_")]] <- dt
      }
    }
  }
  if (length(all_dt) == 0) {
    message("No GSEA results found for database ", db, " in ", io$gsea_results_dir)
    next
  }
  all_dt <- rbindlist(all_dt)

  #############################################
  ## One heatmap per perturbation: x = cell type, y = gene set, color = signed_log_qval ##
  #############################################
  # Union of gene sets with |signed_log_qval| > 5 in at least one cell type (per perturbation)
  # Up = positive (red), down = negative (blue)

  max.val <- 25

  for (pert in opts$perturbations) {
    dt_pert <- all_dt[perturbation == pert]
    if (nrow(dt_pert) == 0) next

    gene_sets_union <- unique(dt_pert[abs(log_qval) > opts$log_qval_threshold, gene_set])
    
    if (length(gene_sets_union) < 3) next

    dt_pert <- dt_pert[gene_set %in% gene_sets_union]
    # If gene set appears in both up and down, keep the value with larger magnitude
    heatmap_dt <- dt_pert[, .(signed_log_qval = signed_log_qval[which.max(abs(signed_log_qval))]), by = .(gene_set, cell_type)]
    heatmap_wide <- dcast(heatmap_dt, gene_set ~ cell_type, value.var = "signed_log_qval")
    mat <- as.matrix(heatmap_wide[, -1, with = FALSE])
    rownames(mat) <- heatmap_wide$gene_set
    mat[is.na(mat)] <- 0
    mat[mat > max.val] <- max.val
    mat[mat < -max.val] <- -max.val

    if (nrow(mat) == 0 || ncol(mat) == 0) next

    breaks <- seq(-max.val, max.val, length.out = 101)

    pheatmap(
      mat,
      breaks = breaks,
      color = colorRampPalette(c("blue", "gray90", "red"))(100),
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_colnames = TRUE,
      show_rownames = TRUE
      scale = "none",
      main = sprintf("GSEA signed -log10(qvalue) (%s vs SCRAMBLE, %s)", pert, db),
      border_color = "gray30",
      fontsize = 8,
      fontsize_row = 7,
      fontsize_col = 9,
      angle_col = 45,
      filename = file.path(io$outdir, "all_celltypes", sprintf("%s_%s_gsea_heatmap.pdf", pert, db)),
      width = max(8, ncol(mat) * 0.4),
      height = max(5, nrow(mat) * 0.4)
    )
  }
}
