source(here::here("settings.R"))
source(here::here("utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description = "Plot GSEA results (perturbation vs SCRAMBLE)")
p$add_argument("--cell_type", type = "character", default = "CCL-210", help = "Cell type to plot")
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

# Options
opts$cell_type <- args$cell_type
opts$perturbations <- setdiff(opts$perturbations, "SCRAMBLE")

# I/O
io$gsea_results_dir <- file.path(io$basedir, "results/gene_signatures/gsea_differential/perturbation")
# io$gsea_go_results <- file.path(io$gsea_results_dir, sprintf("gsea_go_list_%s_vs_SCRAMBLE.rds", opts$cell_type))
io$gsea_msigdb_results <- file.path(io$gsea_results_dir, sprintf("gsea_MSigDB_list_%s_vs_SCRAMBLE.rds", opts$cell_type))
io$gsea_mitocarta_results <- file.path(io$gsea_results_dir, sprintf("gsea_MitoCarta_list_%s_vs_SCRAMBLE.rds", opts$cell_type))
io$gsea_enrichr_results <- file.path(io$gsea_results_dir, sprintf("gsea_enrichr_list_%s_vs_SCRAMBLE.rds", opts$cell_type))
io$gsea_reactome_results <- file.path(io$gsea_results_dir, sprintf("gsea_Reactome_list_%s_vs_SCRAMBLE.rds", opts$cell_type))
io$gsea_manual_results <- file.path(io$gsea_results_dir, sprintf("gsea_Manual_list_%s_vs_SCRAMBLE.rds", opts$cell_type))

io$outdir <- file.path(io$gsea_results_dir, "pdf")
dir.create(io$outdir, showWarnings = F, recursive = T)

#######################
## Load GSEA results ##
#######################
# Lists keyed by perturbation (clusterProfiler) or by contrast "pert_vs_SCRAMBLE" (enrichR)

# go_results_list <- if (file.exists(io$gsea_go_results)) readRDS(io$gsea_go_results) else list()
msigdb_results_list <- if (file.exists(io$gsea_msigdb_results)) readRDS(io$gsea_msigdb_results) else list()
mitocarta_results_list <- if (file.exists(io$gsea_mitocarta_results)) readRDS(io$gsea_mitocarta_results) else list()
enrichr_results_list <- if (file.exists(io$gsea_enrichr_results)) readRDS(io$gsea_enrichr_results) else list()
reactome_results_list <- if (file.exists(io$gsea_reactome_results)) readRDS(io$gsea_reactome_results) else list()
manual_results_list <- if (file.exists(io$gsea_manual_results)) readRDS(io$gsea_manual_results) else list()

dir.create(file.path(io$outdir, opts$cell_type), showWarnings = F, recursive = T)

#############
## Plot GO ##
#############

# if (length(go_results_list) > 0) {
#   # pert <- opts$perturbations[1]; k <- "up"
#   for (pert in opts$perturbations) {
#     p_list <- list()
#     for (k in c("up", "down")) {
#       tmp <- go_results_list[[pert]][[k]]
#       if (is.null(tmp)) {
#         to.plot <- data.table(Description = "NA", Count = 1, qvalue = 1)
#       } else {
#         to.plot <- tmp@result %>% as.data.table() %>% head(n = 25) %>%
#           setorder(-qvalue) %>%
#           .[, Description := factor(Description, levels = Description)]
#       }
#       to.plot[, log_qval := -log10(qvalue)]
#       p_list[[k]] <- ggplot(to.plot, aes(x = log_qval, y = Description)) +
#         geom_point(aes(fill = log_qval), size = 4, shape = 21) +
#         scale_fill_gradient(low = "gray80", high = "red") +
#         labs(x = "Enrichment score (log10 q-value)", y = "", title = sprintf("%s-regulated in %s (background SCRAMBLE)", k, pert)) +
#         theme_bw() +
#         theme(axis.text.x = element_text(color = "black"), axis.text.y = element_text(color = "black", size = rel(0.85)))
#     }
#     p <- cowplot::plot_grid(plotlist = p_list, nrow = 1)
#     pdf(file.path(io$outdir, opts$cell_type, sprintf("GO_gsea_differential_%s_vs_SCRAMBLE_%s.pdf", pert, opts$cell_type)), width = 15, height = 6)
#     print(p)
#     dev.off()
#   }
# }


#################
## Plot msigdb ##
#################

if (length(msigdb_results_list) > 0) {
  
  for (pert in opts$perturbations) {
    p_list <- list()
    for (k in c("up", "down")) {
      tmp <- msigdb_results_list[[pert]][[k]]
      if (is.null(tmp)) {
        to.plot <- data.table(Description = "NA", Count = 1, qvalue = 1)
      } else {
        to.plot <- tmp@result %>% as.data.table() %>% head(n = 25) %>%
          setorder(-qvalue) %>%
          .[, Description := gsub("HALLMARK_", "", Description)] %>%
          .[, Description := factor(Description, levels = Description)]
      }
      to.plot[, log_qval := -log10(qvalue)]
      p_list[[k]] <- ggplot(to.plot, aes(x = log_qval, y = Description)) +
        geom_point(aes(fill = log_qval), size = 4, shape = 21) +
        scale_fill_gradient(low = "gray80", high = "red") +
        labs(x = "Enrichment score (log10 q-value)", y = "", title = sprintf("%s-regulated in %s (background SCRAMBLE)", k, pert)) +
        theme_bw() +
        coord_cartesian(xlim = c(0, max(5, max(to.plot$log_qval, na.rm = TRUE)))) +
        theme(axis.text.x = element_text(color = "black"), axis.text.y = element_text(color = "black", size = rel(0.85)))
    }
    p <- cowplot::plot_grid(plotlist = p_list, nrow = 1)
    pdf(file.path(io$outdir, opts$cell_type, sprintf("%s_%s_msigdb_gsea_differential.pdf", pert, opts$cell_type)), width = 15, height = 6)
    print(p)
    dev.off()
  }
}

####################
## Plot mitoCarta ##
####################

# if (length(mitocarta_results_list) > 0) {
# for (pert in opts$perturbations) {
#   p_list <- list()
#   for (k in c("up", "down")) {
#     tmp <- mitocarta_results_list[[pert]][[k]]
#     if (is.null(tmp)) {
#       to.plot <- data.table(Description = "NA", Count = 1, qvalue = 1)
#     } else {
#       to.plot <- tmp@result %>% as.data.table() %>% head(n = 25) %>%
#         setorder(-qvalue) %>%
#         .[, Description := factor(Description, levels = Description)]
#     }
#     to.plot[, log_qval := -log10(qvalue)]
#     p_list[[k]] <- ggplot(to.plot, aes(x = log_qval, y = Description)) +
#       geom_point(aes(fill = log_qval), size = 4, shape = 21) +
#       scale_fill_gradient(low = "gray80", high = "red") +
#       labs(x = "Enrichment score (log10 q-value)", y = "", title = sprintf("%s-regulated in %s (background SCRAMBLE)", k, pert)) +
#       theme_bw() +
#       coord_cartesian(xlim = c(0, max(5, max(to.plot$log_qval, na.rm = TRUE)))) +
#       theme(axis.text.x = element_text(color = "black"), axis.text.y = element_text(color = "black", size = rel(0.85)))
#   }
#   p <- cowplot::plot_grid(plotlist = p_list, nrow = 1)
#   pdf(file.path(io$outdir, opts$cell_type, sprintf("%s_%s_mitocarta_gsea_differential.pdf", pert, opts$cell_type)), width = 15, height = 6)
#   print(p)
#   dev.off()
# }
# }

##################
## Plot EnrichR ##
##################
# enrichr list: [[db]][[pert_vs_SCRAMBLE]][[up/down]]

opts$databases_enrichr <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023", "Reactome_2022", "MSigDB_Hallmark_2020")

for (db in opts$databases_enrichr) {
  if (length(enrichr_results_list) == 0 || !db %in% names(enrichr_results_list)) next
  for (pert in opts$perturbations) {
    contrast_name <- sprintf("%s_vs_SCRAMBLE", pert)
    if (!contrast_name %in% names(enrichr_results_list[[db]])) next
    p_list <- list()
    for (k in c("up", "down")) {
      tmp <- enrichr_results_list[[db]][[contrast_name]][[k]]
      if (is.null(tmp)) {
        to.plot <- data.table(Description = "NA", Count = 1, pvalue = 1)
      } else {
        to.plot <- tmp %>% head(n = 25) %>%
          setnames("gene_signature", "Description") %>%
          setorder(-pvalue) %>%
          .[, Description := factor(Description, levels = Description)]
      }
      to.plot[, log_qval := -log10(pvalue)]
      p_list[[k]] <- ggplot(to.plot, aes(x = log_qval, y = Description)) +
        geom_point(aes(fill = log_qval), size = 4, shape = 21) +
        scale_fill_gradient(low = "gray80", high = "red") +
        labs(x = "Enrichment score (log10 q-value)", y = "", title = sprintf("%s-regulated in %s (background SCRAMBLE)", k, pert)) +
        theme_bw() +
        coord_cartesian(xlim = c(0, max(5, max(to.plot$log_qval, na.rm = TRUE)))) +
        theme(axis.text.x = element_text(color = "black"), axis.text.y = element_text(color = "black", size = rel(0.85)))
    }
    p <- cowplot::plot_grid(plotlist = p_list, nrow = 1)
    pdf(file.path(io$outdir, opts$cell_type, sprintf("%s_%s_enrichr_gsea_differential.pdf", pert, opts$cell_type)), width = 15, height = 6)
    print(p)
    dev.off()
  }
}

############################
## Plot Reactome pathways ##
############################

for (pert in opts$perturbations) {
  if (length(reactome_results_list) == 0) next
  p_list <- list()
  for (k in c("up", "down")) {
    tmp <- reactome_results_list[[pert]][[k]]
    if (is.null(tmp)) {
      to.plot <- data.table(Description = "NA", Count = 1, qvalue = 1)
    } else {
      to.plot <- tmp@result %>% as.data.table() %>% head(n = 25) %>%
        setorder(-qvalue) %>%
        .[, Description := factor(Description, levels = Description)]
    }
    to.plot[, log_qval := -log10(qvalue)]
    p_list[[k]] <- ggplot(to.plot, aes(x = log_qval, y = Description)) +
      geom_point(aes(fill = log_qval), size = 4, shape = 21) +
      scale_fill_gradient(low = "gray80", high = "red") +
      labs(x = "Enrichment score (log10 q-value)", y = "", title = sprintf("%s-regulated in %s (background SCRAMBLE)", k, pert)) +
      theme_bw() +
      coord_cartesian(xlim = c(0, max(5, max(to.plot$log_qval, na.rm = TRUE)))) +
      theme(axis.text.x = element_text(color = "black"), axis.text.y = element_text(color = "black", size = rel(0.85)))
  }
  p <- cowplot::plot_grid(plotlist = p_list, nrow = 1)
  pdf(file.path(io$outdir, opts$cell_type, sprintf("%s_%s_reactome_gsea_differential.pdf", pert, opts$cell_type)), width = 15, height = 6)
  print(p)
  dev.off()
}


############################
## Plot manual signatures ##
############################

for (pert in opts$perturbations) {
  if (length(manual_results_list) == 0) next
  p_list <- list()
  for (k in c("up", "down")) {
    tmp <- manual_results_list[[pert]][[k]]
    if (is.null(tmp)) {
      to.plot <- data.table(Description = "NA", Count = 1, qvalue = 1)
    } else {
      to.plot <- tmp@result %>% as.data.table() %>% head(n = 25) %>%
        setorder(-qvalue) %>%
        .[, Description := factor(Description, levels = Description)]
    }
    to.plot[, log_qval := -log10(qvalue)]
    p_list[[k]] <- ggplot(to.plot, aes(x = log_qval, y = Description)) +
      geom_point(aes(fill = log_qval), size = 4, shape = 21) +
      scale_fill_gradient(low = "gray80", high = "red") +
      labs(x = "Enrichment score (log10 q-value)", y = "", title = sprintf("%s-regulated in %s\n(background SCRAMBLE)", k, pert)) +
      theme_bw() +
      coord_cartesian(xlim = c(0, max(5, max(to.plot$log_qval, na.rm = TRUE)))) +
      theme(axis.text.x = element_text(color = "black"), axis.text.y = element_text(color = "black", size = rel(0.85)), plot.title = element_text(color = "black", size = rel(0.95)))
  }
  p <- cowplot::plot_grid(plotlist = p_list, nrow = 1)
  pdf(file.path(io$outdir, opts$cell_type, sprintf("%s_%s_manual_gsea_differential.pdf", pert, opts$cell_type)), width = 15, height = 6)
  print(p)
  dev.off()
}


################################################################
## Plot p-value across all conditions for individual pathways ##
################################################################

# outdir <- file.path(io$outdir,"all_conditions"); dir.create(outdir, showWarnings = F)
# 
# # i <- opts$contrasts[1]; j <- "down"
# tmp.dt <- list()
# for (i in opts$contrasts) {
#   for (j in c("up","down")) {
#     # tmp <- hr3_results_list[[i]][[j]]
#     tmp <- msigdb_results_list[[i]][[j]]
#     # tmp <- go_results_list[[i]][[j]]
#     # tmp <- mitocarta_results_list[[i]][[j]]
#     if (is.null(tmp)) { 
#       tmp.dt[[sprintf("%s_%s",i,j)]] <- data.table(Description="NA", Count=1, qvalue=1, log_qval=0) %>%
#         .[,class:=sprintf("Genes %s-regulated",j)] %>% .[,contrast:=i]
#     } else {
#       tmp.dt[[sprintf("%s_%s",i,j)]] <- as.data.table(tmp@result) %>% .[,log_qval:=-log10(qvalue)] %>% 
#         .[,Description:=gsub("HALLMARK_","",Description)] %>%
#         .[,class:=sprintf("Genes %s-regulated",j)] %>% .[,contrast:=i]
#     }
#   }
# }
# 
# tmp.dt <- rbindlist(tmp.dt,fill=TRUE) %>% setnames("Description","pathway")
# tmp.dt[,contrast:=factor(contrast,levels=opts$contrasts)]
# tmp.dt[,class:=factor(class, levels=c("Genes up-regulated","Genes down-regulated"))]
# 
# pathways.to.plot <- unique(tmp.dt$pathway)
# # pathways.to.plot <- grep("RIBOSOME",pathways.to.plot,value=T)
# 
# # i <- pathways.to.plot[1]
# for (i in pathways.to.plot) {
#   
#   to.plot <- data.table(expand.grid(opts$constructs,c("Genes up-regulated","Genes down-regulated"))) %>%
#     setnames(c("contrast","class")) %>% .[,pathway:=i] %>%
#     merge(tmp.dt[pathway==i,c("pathway","contrast","class","log_qval")], all.x=TRUE, by=c("pathway","contrast","class")) %>%
#     .[is.na(log_qval) | log_qval==0,log_qval:=0.01]
#   
#   xlim <- max(5,max(to.plot[pathway==i]$log_qval))
#   
#   p <- ggplot(to.plot[pathway==i], aes(x=log_qval, y=contrast)) +
#     # geom_bar(aes(fill=log_qval), stat="identity", color="black") +
#     geom_bar(fill="gray70", stat="identity", color="black") +
#     # scale_fill_gradient(low = "gray80", high = "purple") +
#     # labs(x="Number of diff. genes in the pathway", y="", title=sprintf("%s-regulated in %s",j,i)) +
#     geom_vline(xintercept=2, linetype="dashed") +
#     facet_wrap(~class, nrow=1, scales="fixed") +
#     scale_y_discrete(drop=FALSE) +
#     coord_cartesian(xlim=c(0,xlim)) +
#     labs(x="Enrichment score (log10 q-value)", y="", title=i) +
#     theme_classic() +
#     theme(
#       plot.title = element_text(hjust=0.5, size=rel(1.10), color="black"),
#       legend.position = "none",
#       axis.text.x = element_text(color="black", size=rel(1)),
#       strip.text = element_text(color="black", size=rel(1.05)),
#       axis.title.x = element_text(color="black", size=rel(1.0)),
#       axis.text.y = element_text(color="black", size=rel(1.0)),
#     )
#   
#   pdf(file.path(outdir,sprintf("%s_gsea_differential_all_conditions_HDF_oskm_mptr_rna_reik.pdf",i)), width=6, height=5.25)
#   print(p)
#   dev.off()
# }


####################
## Selected plots ##
####################

# if (length(msigdb_results_list) > 0 && length(opts$perturbations) > 0) {
#   i <- opts$perturbations[1]
#   j <- "down"
#   tmp <- msigdb_results_list[[i]][[j]]
#   if (!is.null(tmp)) {
#     to.plot <- tmp@result %>% as.data.table() %>% head(n = 15) %>%
#       setorder(-qvalue) %>%
#       .[, log_qval := -log10(qvalue)] %>%
#       .[, Description := gsub("HALLMARK_", "", Description)] %>%
#       .[, Description := factor(Description, levels = Description)]

#     p <- ggplot(to.plot, aes(x = log_qval, y = Description)) +
#       geom_point(aes(fill = log_qval), size = 4, shape = 21) +
#       scale_fill_gradient(low = "gray80", high = "red") +
#       coord_cartesian(xlim = c(0, max(5, max(to.plot$log_qval, na.rm = TRUE)))) +
#       labs(x = "Enrichment score (log10 q-value)", y = "", title = sprintf("%s, %s-regulated genes (%s)", i, j, opts$cell_type)) +
#       theme_bw() +
#       theme(
#         legend.position = "none",
#         plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(color = "black"),
#         axis.text.y = element_text(color = "black", size = rel(1.0)),
#       )

#     pdf(file.path(io$outdir, opts$cell_type, sprintf("msigdb_selected_%s_%s_%s.pdf", i, j, opts$cell_type)), width = 15, height = 6)
#     print(p)
#     dev.off()
#   }
# }
