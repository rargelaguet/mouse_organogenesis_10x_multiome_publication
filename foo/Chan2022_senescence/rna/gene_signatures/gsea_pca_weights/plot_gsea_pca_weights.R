source(here::here("settings.R"))
source(here::here("utils.R"))

library(enrichplot)

#####################
## Define settings ##
#####################

# I/O
io$pca_weights <- file.path(io$basedir,"results/rna/dimensionality_reduction/PCA_weights.txt.gz")
io$gsea_go_results <- file.path(io$basedir,"results/rna/gene_signatures/pca_weights/gsea_go_list.rds")
io$gsea_msigdb_results <- file.path(io$basedir,"results/rna/gene_signatures/pca_weights/gsea_mgsidb_list.rds")
io$gsea_mitocarta_results <- file.path(io$basedir,"results/rna/gene_signatures/pca_weights/gsea_mitocarta_list.rds")
io$gsea_enrichr_results <- file.path(io$basedir,"results/rna/gene_signatures/pca_weights/gsea_enrichr_list_top_150_genes.rds")
io$outdir <- file.path(io$basedir,"results/rna/gene_signatures/pca_weights/pdf"); dir.create(io$outdir, showWarnings = F)

# Options
opts$PCs.to.use <- c("PC1","PC2")

######################
## Load PCA weights ##
######################

pca_weights.dt <- fread(io$pca_weights)

#######################
## Load GSEA results ##
#######################

msigdb_results_list <- readRDS(io$gsea_msigdb_results)
go_results_list <- readRDS(io$gsea_go_results)
mitocarta_results_list <- readRDS(io$gsea_mitocarta_results)
enrichr_results_list <- readRDS(io$gsea_enrichr_results)

#############
## Plot GO ##
#############

# i <- PCs.to.use[1]; j <- "up"
for (i in opts$PCs.to.use) {
  p_list <- list()
  for (j in c("up","down")) {
    tmp <- go_results_list[[i]][[j]]
    if (is.null(tmp)) { 
      to.plot <- data.table(Description="NA", Count=1, qvalue=1)
    } else {
      to.plot <- tmp@result %>% as.data.table %>% head(n=25) %>%
        setorder(-qvalue) %>%
        .[,Description:=gsub("HALLMARK_","",Description)] %>%
        .[,Description:=factor(Description, levels=Description)]
    }
    p_list[[j]] <- ggplot(to.plot, aes(x=Count, y=Description)) +
      geom_point(aes(fill=-log10(qvalue)), size=4, shape=21) +
      scale_fill_gradient(low = "gray80", high = "red") + 
      labs(x="Number of diff. genes in the pathway", y="", title=sprintf("%s-regulated in %s",j,i)) +
      theme_bw() +
      theme(
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black", size=rel(0.85)),
      )
    
  }
  p <- cowplot::plot_grid(plotlist=p_list, nrow=1)
  
  pdf(file.path(io$outdir,sprintf("GO_%s_BJ_passage_CHAMP.pdf",i)), width=15, height=6)
  print(p)
  dev.off()
}


#################
## Plot msigdb ##
#################

# i <- opts$PCs.to.use[1]; j <- "down"
for (i in opts$PCs.to.use) {
  p_list <- list()
  for (j in c("up","down")) {
    tmp <- msigdb_results_list[[i]][[j]]
    if (is.null(tmp)) { 
      to.plot <- data.table(Description="NA", Count=1, qvalue=1)
    } else {
      to.plot <- tmp@result %>% as.data.table %>% head(n=25) %>%
        setorder(-qvalue) %>%
        .[,Description:=gsub("HALLMARK_","",Description)] %>%
        .[,Description:=factor(Description, levels=Description)]
    }
    p_list[[j]] <- ggplot(to.plot, aes(x=Count, y=Description)) +
      geom_point(aes(fill=-log10(qvalue)), size=4, shape=21) +
      scale_fill_gradient(low = "gray80", high = "red") + 
      labs(x="Number of diff. genes in the pathway", y="", title=sprintf("%s-regulated in %s",j,i)) +
      theme_bw() +
      theme(
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black", size=rel(0.85)),
      )
    
  }
  p <- cowplot::plot_grid(plotlist=p_list, nrow=1)
  
  pdf(file.path(io$outdir,sprintf("msigdb_%s_BJ_passage_CHAMP.pdf",i)), width=15, height=6)
  print(p)
  dev.off()
}

####################
## Plot mitoCarta ##
####################

# i <- opts$PCs.to.use[1]; j <- "down"
for (i in opts$PCs.to.use) {
  p_list <- list()
  for (j in c("up","down")) {
    tmp <- mitocarta_results_list[[i]][[j]]
    if (is.null(tmp)) { 
      to.plot <- data.table(Description="NA", Count=1, qvalue=1)
    } else {
      to.plot <- tmp@result %>% as.data.table %>% head(n=30) %>%
        setorder(-qvalue) %>%
        .[,Description:=gsub("HALLMARK_","",Description)] %>%
        .[,Description:=factor(Description, levels=Description)]
    }
    p_list[[j]] <- ggplot(to.plot, aes(x=Count, y=Description)) +
      geom_point(aes(fill=-log10(qvalue)), size=4, shape=21) +
      scale_fill_gradient(low = "gray80", high = "red") + 
      labs(x="Number of diff. genes in the pathway", y="", title=sprintf("%s-regulated in %s",j,i)) +
      theme_bw() +
      theme(
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black", size=rel(0.85)),
      )
    
  }
  p <- cowplot::plot_grid(plotlist=p_list, nrow=1)
  
  pdf(file.path(io$outdir,sprintf("mitocarta_%s_BJ_passage_CHAMP.pdf",i)), width=17, height=7)
  print(p)
  dev.off()
}


##################
## Plot EnrichR ##
##################

opts$databases <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023", "Reactome_2022", "MSigDB_Hallmark_2020", "BioCarta_2016", "KEGG_2021_Human")

# i <- opts$PCs.to.use[1]; k <- opts$databases[1]; j <- "down"
for (k in opts$databases) { 
  for (i in opts$PCs.to.use) {
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
        labs(x="Log p-value", y="", title=sprintf("%s-regulated in %s",j,i)) +
        theme_bw() +
        theme(
          axis.text.x = element_text(color="black"),
          axis.text.y = element_text(color="black", size=rel(0.85)),
        )
      
    }
    p <- cowplot::plot_grid(plotlist=p_list, nrow=1)
    
    pdf(file.path(io$outdir,sprintf("enrichr_%s_%s_BJ_passage_CHAMP.pdf",i,k)), width=17, height=7)
    print(p)
    dev.off()
  }
} 
