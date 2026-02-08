source(here::here("settings.R"))
source(here::here("utils.R"))

library(enrichplot)

#####################
## Define settings ##
#####################

# I/O
io$gsea_results <- file.path(io$basedir,"results/rna/gene_signatures/GSEA_differential/enrichR/gsea_enrichr_list.rds")
io$outdir <- file.path(io$basedir,"results/rna/gene_signatures/GSEA_differential/enrichR/pdf"); dir.create(io$outdir, showWarnings = F)

# Options
opts$group_variable <- "sample"

opts$groups <- c( 
  # "OS_8OFF",
  # "OS_6OFF-2ON",
  # "OS_4OFF-4ON",
  # "OS_8ON",
  # "OS_2OFF-2ON-4OFF",
  # "OS_4ON-4OFF",
  # "OSK_8OFF",
  # "OSK_4OFF-4ON",
  # "OSK_8ON",
  # "OSK_4ON-4OFF",
  "OSKM_8OFF",
  # "OSKM_6OFF-2ON",
  "OSKM_4OFF-4ON"
  # "OSKM_8ON",
  # "OSKM_2OFF-2ON-4OFF",
  # "OSKM_4ON-4OFF"
)

opts$databases <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023", "Reactome_2022","MSigDB_Hallmark_2020", "BioCarta_2016", "KEGG_2021_Human")

#######################
## Load GSEA results ##
#######################

gsea_results_list <- readRDS(io$gsea_results)
names(gsea_results_list)
names(gsea_results_list[[1]])
names(gsea_results_list[[1]][[1]])

#######################
## Plot GSEA results ##
#######################

# db <- opts$databases[1]; i <- "OSKM_8OFF"; j <- "OSKM_4OFF-4ON"; k <- "up"
for (db in opts$databases) {
  
  outdir <- file.path(io$outdir,db); dir.create(outdir, showWarnings = F)
  
  for (i in opts$groups) {
    for (j in opts$groups) {
      if (i!=j) {
        p_list <- list()
        for (k in c("up","down")) {
          tmp <- gsea_results_list[[db]][[sprintf("%s_vs_%s",i,j)]][[k]]
          if (is.null(tmp)) { 
            to.plot <- data.table(Description="NA", Count=1, pval_adj=1)
          } else {
            to.plot <- tmp %>% as.data.table %>% head(n=25) %>%
              setorder(-pval_adj) %>% setnames("gene_signature","Description") %>%
              .[,Description:=factor(Description, levels=Description)] %>%
              .[,log_pval:=-log10(pval_adj)]
          }
          p_list[[k]] <- ggplot(to.plot, aes(x=log_pval, y=Description)) +
            geom_point(aes(fill=log_pval), size=4, shape=21) +
            scale_fill_gradient(low = "gray80", high = "red") + 
            labs(x="FDR-adjusted p-value (log10)", y="", title=sprintf("%s-regulated in %s",k,j)) +
            theme_bw() +
            theme(
              legend.position = "none",
              axis.text.x = element_text(color="black"),
              axis.text.y = element_text(color="black", size=rel(0.80)),
            )
          
        }
        p <- cowplot::plot_grid(plotlist=p_list, nrow=1)
        
        pdf(file.path(outdir,sprintf("GSEA_enrichR_%s_vs_%s.pdf",i,j)), width=16, height=7)
        print(p)
        dev.off()
      }
    }
  }
}

##########################################
## Scatterplots (comparing OSK vs OSKM) ##
##########################################

conditions.to.plot <- c(
  "4OFF-4ON",
  "8ON",
  "4ON-4OFF"
)

background_condition <- "8OFF"

genes.to.highlight <- c("LY6D","ST6GAL1","PI16","CD70","CD177","CD9","PODXL","ANPEP","A4GALT","CD24")

# i <- conditions.to.plot[1]
for (i in conditions.to.plot) {
  
  gsea_diff_1.dt <- gsea_diff.dt[groupA==sprintf("OSK_%s",background_condition) & groupB==sprintf("OSK_%s",i),c("feature","transcript_type","logFC")]
  gsea_diff_2.dt <- gsea_diff.dt[groupA==sprintf("OSKM_%s",background_condition) & groupB==sprintf("OSKM_%s",i),c("feature","transcript_type","logFC")]
  to.plot <- merge(gsea_diff_1.dt, gsea_diff_2.dt, by=c("feature","transcript_type"), suffixes=c("_A","_B")) %>%
    .[abs(logFC_A)>0.01 & abs(logFC_B)>0.01] %>%
    .[,tmp:=abs(logFC_A)] %>% setorder(-tmp) %>% .[,tmp:=NULL]
  
  to.plot <- to.plot[!feature%in%c("POU5F1","SOX2","KLF4","MYC","exogenous_UTR")]
  
  to.plot <- to.plot[transcript_type%in%c("Surface markers")]
  
  ylim.max <- max(abs(to.plot$logFC_A))
  xlim.max <- max(abs(to.plot$logFC_B))
  lim.max <- max(xlim.max,ylim.max)
  
  p <- ggplot(to.plot, aes(x=logFC_A, y=logFC_B)) +
    geom_jitter(size=0.25, width=0.01, height=0.01, color="gray30", data=to.plot[!(abs(logFC_B)>=1.5 | abs(logFC_A)>=1.5)][!feature%in%genes.to.highlight]) +
    geom_jitter(size=0.55, width=0.01, height=0.01, color="gray30", data=to.plot[abs(logFC_B)>=1.5 | abs(logFC_A)>=1.5][!feature%in%genes.to.highlight]) +
    # geom_jitter(size=0.55, width=0.01, height=0.01, color="purple", data=to.plot[feature%in%genes.to.highlight]) +
    stat_cor(method = "pearson") +
    coord_cartesian(xlim=c(-lim.max,lim.max), ylim=c(-lim.max,lim.max)) +
    geom_segment(x=-xlim.max, xend=xlim.max, y=0, yend=0, color="orange", linewidth=0.1) +
    geom_segment(x=0, xend=0, y=-ylim.max, yend=ylim.max, color="orange", linewidth=0.1) +
    # ggrepel::geom_text_repel(data=head(to.plot[transcript_type=="Genes" & logFC_B>=1.5 & logFC_A>=1.5 & !feature%in%genes.to.highlight],n=15), aes(x=logFC_A, y=logFC_B, label=feature), color="black", size=3, segment.size=0.1, max.overlaps = Inf) +
    # ggrepel::geom_text_repel(data=head(to.plot[transcript_type=="Genes" & logFC_B<(-1.5) & logFC_A<=(-1.5) & !feature%in%genes.to.highlight],n=15), aes(x=logFC_A, y=logFC_B, label=feature), color="black", size=3, segment.size=0.1, max.overlaps = Inf) +
    # ggrepel::geom_text_repel(data=to.plot[feature%in%genes.to.highlight], aes(x=logFC_A, y=logFC_B, label=feature), color="purple", size=3, segment.size=0.1, max.overlaps = Inf) +
    labs(x=sprintf("logFC (OSK, %s vs %s)",i,background_condition), y=sprintf("logFC (OSKM, %s vs %s)",i,background_condition), title="") +
    # facet_wrap(~transcript_type, nrow=1) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(size=rel(1.0), color="black"),
      axis.text.y = element_text(size=rel(1.0), color="black"),
      legend.position = "none",
      legend.title = element_blank()
    )
  
  # pdf(file.path(io$outdir.scatter,sprintf("scatterplot_OSK_vs_OSKM_%s_HDF_multiome_CHAMP.pdf",i)), width=16, height=7)
  pdf(file.path(io$outdir.scatter,sprintf("scatterplot_surface_markers_OSK_vs_OSKM_%s_HDF_multiome_CHAMP.pdf",i)), width=8, height=7)
  print(p)
  dev.off()
}


