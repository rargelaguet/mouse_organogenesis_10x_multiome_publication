source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# Options
opts$condition <- c(
  # "hTERT", 
  "replicative_senescence"
  # "radiation_senescence"
)

# I/O
io$limma_results <- file.path(io$basedir,sprintf("rna/results/limma/limma_results_%s.txt.gz",opts$condition))
io$gsea_go_results <- file.path(io$basedir,sprintf("rna/results/gene_signatures/gsea_limma/%s/gsea_go_list.rds",opts$condition))
io$gsea_msigdb_results <- file.path(io$basedir,sprintf("rna/results/gene_signatures/gsea_limma/%s/gsea_msigdb_list.rds",opts$condition))
io$gsea_mitocarta_results <- file.path(io$basedir,sprintf("rna/results/gene_signatures/gsea_limma/%s/gsea_mitocarta_list.rds",opts$condition))
io$gsea_enrichr_results <- file.path(io$basedir,sprintf("rna/results/gene_signatures/gsea_limma/%s/gsea_enrichr_list.rds",opts$condition))
io$gsea_hr3_results <- file.path(io$basedir,sprintf("rna/results/gene_signatures/gsea_limma/%s/gsea_hr3_list.rds",opts$condition))
io$gsea_reactome_results <- file.path(io$basedir,sprintf("rna/results/gene_signatures/gsea_limma/%s/gsea_reactome_list.rds",opts$condition))
io$gsea_decoupler_results <- file.path(io$basedir,sprintf("rna/results/gene_signatures/gsea_limma/%s/gsea_decoupler_list.rds",opts$condition))
io$outdir <- file.path(io$basedir,sprintf("rna/results/gene_signatures/gsea_limma/%s/pdf",opts$condition)); dir.create(io$outdir, showWarnings = F)

########################
## Load limma results ##
########################

limma.list <- list("time"=fread(io$limma_results) %>% setnames("symbol","feature"))

opts$contrasts <- names(limma.list)
# opts$contrasts <- names(limma.list)[2]

#######################
## Load GSEA results ##
#######################

go_results_list <- readRDS(io$gsea_go_results)
msigdb_results_list <- readRDS(io$gsea_msigdb_results)
mitocarta_results_list <- readRDS(io$gsea_mitocarta_results)
enrichr_results_list <- readRDS(io$gsea_enrichr_results)
hr3_results_list <- readRDS(io$gsea_hr3_results)
reactome_results_list <- readRDS(io$gsea_reactome_results)
decoupler_results_list <- readRDS(io$gsea_decoupler_results)

#############
## Plot GO ##
#############

# i <- opts$contrasts[1]; j <- "down"
for (i in opts$contrasts) {
  p_list <- list()
  for (j in c("up","down")) {
    tmp <- go_results_list[[i]][[j]]
    if (is.null(tmp)) { 
      to.plot <- data.table(Description="NA", Count=1, qvalue=1)
    } else {
      to.plot <- tmp@result %>% as.data.table %>% head(n=25) %>%
        setorder(-qvalue) %>%
        .[,Description:=factor(Description, levels=Description)] %>%
        .[,log_qval:=-log10(qvalue)]
    }
    p_list[[j]] <- ggplot(to.plot, aes(x=log_qval, y=Description)) +
      geom_point(aes(fill=log_qval), size=4, shape=21) +
      scale_fill_gradient(low = "gray80", high = "red") + 
      # labs(x="Number of diff. genes in the pathway", y="", title=sprintf("%s-regulated in %s",j,i)) +
      labs(x="Enrichment score (log10 q-value)", y="", title=sprintf("%s-regulated in %s",j,i)) +
      theme_bw() +
      theme(
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black", size=rel(0.85)),
      )
    
  }
  p <- cowplot::plot_grid(plotlist=p_list, nrow=1)
  
  pdf(file.path(io$outdir,sprintf("GO_gsea_limma_%s_Chan2022_senescence.pdf",i)), width=15, height=6)
  print(p)
  dev.off()
}


#################
## Plot msigdb ##
#################

# i <- opts$contrasts[1]; j <- "down"
for (i in opts$contrasts) {
  p_list <- list()
  for (j in c("up","down")) {
    tmp <- msigdb_results_list[[i]][[j]]
    if (is.null(tmp)) { 
      to.plot <- data.table(Description="NA", Count=1, qvalue=1)
    } else {
      to.plot <- tmp@result %>% as.data.table %>% head(n=25) %>%
        setorder(-qvalue) %>%
        .[,log_qval:=-log10(qvalue)] %>%
        .[,Description:=gsub("HALLMARK_","",Description)] %>%
        .[,Description:=factor(Description, levels=Description)] 
    }
    p_list[[j]] <- ggplot(to.plot, aes(x=log_qval, y=Description)) +
      geom_point(aes(fill=log_qval), size=4, shape=21) +
      scale_fill_gradient(low = "gray80", high = "red") + 
      # labs(x="Number of diff. genes in the pathway", y="", title=sprintf("%s-regulated in %s",j,i)) +
      labs(x="Enrichment score (log10 q-value)", y="", title=sprintf("%s-regulated in %s",j,i)) +
      theme_bw() +
      theme(
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black", size=rel(0.85)),
      )
    
  }
  p <- cowplot::plot_grid(plotlist=p_list, nrow=1)
  
  pdf(file.path(io$outdir,sprintf("msigdb_gsea_limma_%s_Chan2022_senescence.pdf",i)), width=15, height=6)
  print(p)
  dev.off()
}


####################
## Plot mitoCarta ##
####################

# i <- opts$contrasts[1]; j <- "down"
for (i in opts$contrasts) {
  p_list <- list()
  for (j in c("up","down")) {
    tmp <- mitocarta_results_list[[i]][[j]]
    if (is.null(tmp)) { 
      to.plot <- data.table(Description="NA", Count=1, qvalue=1)
    } else {
      to.plot <- tmp@result %>% as.data.table %>% head(n=30) %>%
        setorder(-qvalue) %>%
        .[,log_qval:=-log10(qvalue)] %>%
        .[,Description:=factor(Description, levels=Description)]
    }
    p_list[[j]] <- ggplot(to.plot, aes(x=log_qval, y=Description)) +
      geom_point(aes(fill=log_qval), size=4, shape=21) +
      scale_fill_gradient(low = "gray80", high = "red") + 
      # labs(x="Number of diff. genes in the pathway", y="", title=sprintf("%s-regulated in %s",j,i)) +
      labs(x="Enrichment score (log10 q-value)", y="", title=sprintf("%s-regulated in %s",j,i)) +
      theme_bw() +
      theme(
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black", size=rel(0.85)),
      )
    
  }
  p <- cowplot::plot_grid(plotlist=p_list, nrow=1)
  
  pdf(file.path(io$outdir,sprintf("mitocarta_gsea_limma_%s_Chan2022_senescence.pdf",i)), width=17, height=7)
  print(p)
  dev.off()
}


##################
## Plot EnrichR ##
##################

opts$databases <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023", "Reactome_2022", "MSigDB_Hallmark_2020", "BioCarta_2016", "KEGG_2021_Human")

# i <- opts$contrasts[1]; k <- opts$databases[1]; j <- "down"
for (k in opts$databases) { 
  for (i in opts$contrasts) {
    p_list <- list()
    for (j in c("up","down")) {
      tmp <- enrichr_results_list[[k]][[i]][[j]]
      if (is.null(tmp)) { 
        to.plot <- data.table(Description="NA", Count=1, qvalue=1)
      } else {
        to.plot <- tmp %>% head(n=30) %>%
          setorder(-pval_adj) %>% setnames("gene_signature","Description") %>%
          .[,Description:=factor(Description, levels=Description)] %>%
          .[,log_pval:=-log10(pval_adj)]
      }
      p_list[[j]] <- ggplot(to.plot, aes(x=log_pval, y=Description)) +
        geom_point(aes(fill=log_pval), size=4, shape=21) +
        scale_fill_gradient(low = "gray80", high = "red") + 
        labs(x="Enrichment score (log10 q-value)", y="", title=sprintf("%s-regulated in %s",j,i)) +
        theme_bw() +
        theme(
          axis.text.x = element_text(color="black"),
          axis.text.y = element_text(color="black", size=rel(0.85)),
        )
      
    }
    p <- cowplot::plot_grid(plotlist=p_list, nrow=1)
    
    pdf(file.path(io$outdir,sprintf("enrichr_%s_gsea_limma_%s_Chan2022_senescence.pdf",i,k)), width=17, height=7)
    print(p)
    dev.off()
  }
} 

#######################
## Plot HR3 pathways ##
#######################

# i <- opts$contrasts[1]; j <- "down"
for (i in opts$contrasts) {
  p_list <- list()
  for (j in c("up","down")) {
    tmp <- hr3_results_list[[i]][[j]]
    if (is.null(tmp)) { 
      to.plot <- data.table(Description="NA", Count=1, qvalue=1)
    } else {
      to.plot <- tmp@result %>% as.data.table %>% head(n=30) %>%
        setorder(-qvalue) %>%
        .[,log_qval:=-log10(qvalue)] %>%
        .[,Description:=factor(Description, levels=Description)]
    }
    p_list[[j]] <- ggplot(to.plot, aes(x=log_qval, y=Description)) +
      geom_point(aes(fill=log_qval), size=4, shape=21) +
      scale_fill_gradient(low = "gray80", high = "red") + 
      # labs(x="Number of diff. genes in the pathway", y="", title=sprintf("%s-regulated in %s",j,i)) +
      labs(x="Enrichment score (log10 q-value)", y="", title=sprintf("%s-regulated in %s",j,i)) +
      theme_bw() +
      theme(
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black", size=rel(0.85)),
      )
    
  }
  p <- cowplot::plot_grid(plotlist=p_list, nrow=1)
  
  pdf(file.path(io$outdir,sprintf("hr3_gsea_limma_%s_Chan2022_senescence.pdf",i)), width=17, height=7)
  print(p)
  dev.off()
}


############################
## Plot Reactome pathways ##
############################

# i <- opts$contrasts[1]; j <- "down"
for (i in opts$contrasts) {
  p_list <- list()
  for (j in c("up","down")) {
    tmp <- reactome_results_list[[i]][[j]]
    if (is.null(tmp)) { 
      to.plot <- data.table(Description="NA", Count=1, qvalue=1)
    } else {
      to.plot <- tmp@result %>% as.data.table %>% head(n=30) %>%
        setorder(-qvalue) %>%
        .[,log_qval:=-log10(qvalue)] %>%
        .[,Description:=factor(Description, levels=Description)]
    }
    p_list[[j]] <- ggplot(to.plot, aes(x=log_qval, y=Description)) +
      geom_point(aes(fill=log_qval), size=4, shape=21) +
      scale_fill_gradient(low = "gray80", high = "red") + 
      # labs(x="Number of diff. genes in the pathway", y="", title=sprintf("%s-regulated in %s",j,i)) +
      labs(x="Enrichment score (log10 q-value)", y="", title=sprintf("%s-regulated in %s",j,i)) +
      theme_bw() +
      theme(
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black", size=rel(0.85)),
      )
    
  }
  p <- cowplot::plot_grid(plotlist=p_list, nrow=1)
  
  pdf(file.path(io$outdir,sprintf("reactome_gsea_limma_%s_Chan2022_senescence.pdf",i)), width=17, height=7)
  print(p)
  dev.off()
}


#############################
## Plot DecoupleR pathways ##
#############################

# i <- opts$contrasts[1]; j <- "down"
for (i in opts$contrasts) {
  p_list <- list()
  for (j in c("up","down")) {
    tmp <- decoupler_results_list[[i]][[j]]
    if (is.null(tmp)) { 
      to.plot <- data.table(Description="NA", Count=1, qvalue=1)
    } else {
      to.plot <- tmp@result %>% as.data.table %>% head(n=30) %>%
        setorder(-qvalue) %>%
        .[,log_qval:=-log10(qvalue)] %>%
        .[,Description:=factor(Description, levels=Description)]
    }
    p_list[[j]] <- ggplot(to.plot, aes(x=log_qval, y=Description)) +
      geom_point(aes(fill=log_qval), size=4, shape=21) +
      scale_fill_gradient(low = "gray80", high = "red") + 
      # labs(x="Number of diff. genes in the pathway", y="", title=sprintf("%s-regulated in %s",j,i)) +
      labs(x="Enrichment score (log10 q-value)", y="", title=sprintf("%s-regulated in %s",j,i)) +
      theme_bw() +
      theme(
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black", size=rel(0.85)),
      )
    
  }
  p <- cowplot::plot_grid(plotlist=p_list, nrow=1)
  
  pdf(file.path(io$outdir,sprintf("decoupler_gsea_limma_%s_Chan2022_senescence.pdf",i)), width=17, height=7)
  print(p)
  dev.off()
}


####################
## Selected plots ##
####################

# # i <- opts$contrasts[2]; j <- "up"
# # i <- opts$contrasts[1]; j <- "down"
# i <- opts$contrasts[3]; j <- "up"
# # tmp <- go_results_list[[i]][[j]]
# tmp <- hr3_results_list[[i]][[j]]
# # tmp <- enrichr_results_list[["MSigDB_Hallmark_2020"]][[i]][[j]]
# 
# # to.plot <- tmp %>% head(n=15) %>%
# #   setorder(-pval_adj) %>% setnames("gene_signature","Description") %>%
# #   .[,Description:=factor(Description, levels=Description)] %>%
# #   .[,log_pval:=-log10(pval_adj)]
# 
# to.plot <- tmp@result %>% as.data.table %>% head(n=13) %>%
#   setorder(-qvalue) %>%
#   .[,Description:=factor(Description, levels=Description)] %>%
#   .[,log_pval:=-log10(qvalue)]
# 
# p <- ggplot(to.plot, aes(x=log_pval, y=Description)) +
#   geom_point(aes(fill=log_pval), size=4, shape=21) +
#   scale_fill_gradient(low = "gray80", high = "red") + 
#   labs(x="Enrichment score (log10 q-value)", y="", title="") +
#   theme_bw() +
#   theme(
#     legend.position = "none",
#     axis.text.x = element_text(color="black"),
#     axis.text.y = element_text(color="black", size=rel(1.05)),
#     axis.title = element_text(color="black", size=rel(1.05))
#   )
# 
# pdf(file.path(io$outdir,"test.pdf"), width=6, height=6)
# print(p)
# dev.off()
