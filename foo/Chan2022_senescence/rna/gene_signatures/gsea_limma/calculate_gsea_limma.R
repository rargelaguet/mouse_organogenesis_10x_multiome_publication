source(here::here("settings.R"))
source(here::here("utils.R"))

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichR)

#####################
## Define settings ##
#####################

# I/O
# io$limma_results <- file.path(io$basedir,"rna/results/limma/limma_results_cell_density.txt.gz")
# io$outdir <- file.path(io$basedir,"rna/results/gene_signatures/gsea_limma/cell_density")

# io$limma_results <- file.path(io$basedir,"rna/results/limma/limma_results_hTERT.txt.gz")
# io$outdir <- file.path(io$basedir,"results/rna/gene_signatures/gsea_limma/hTERT")

io$limma_results <- file.path(io$basedir,"rna/results/limma/limma_results_radiation_senescence.txt.gz")
io$outdir <- file.path(io$basedir,"results/rna/gene_signatures/gsea_limma/radiation_senescence")

# io$limma_results <- file.path(io$basedir,"rna/results/limma/limma_results_replicative_senescence.txt.gz")
# io$outdir <- file.path(io$basedir,"rna/results/gene_signatures/gsea_limma/replicative_senescence")


dir.create(io$outdir, showWarnings = F, recursive = T)
# io$outdir.scatterplots <- file.path(io$basedir,"results/rna/differential/metacells/alias/pdf/scatterplots"); dir.create(io$outdir.scatterplots, showWarnings = F)

# Options
opts$pval <- 1e-4
opts$coef <- 0.05

# Parallel processing
library(furrr)
plan(multisession, workers = 3)

########################
## Load gene metadata ##
########################

gene_metadata.dt <- fread(io$gene_metadata)

########################
## Load limma results ##
########################

limma.list <- list("time"=fread(io$limma_results) %>% setnames("symbol","feature"))

.[,sig:=pval_adj<opts$pval & abs(coef)>=opts$coef] %>%
  setorder(pval_adj)

contrasts <- names(limma.list)

##########################
## Load gene signatures ##
##########################

source(here::here("rna/gene_signatures/load_gene_signatures.R"))

#########################
## GSEA using GO terms ##
#########################

# i <- contrasts[1]
gsea_go_list <- contrasts %>% future_map(function(i) {
  
  gsea_list <- list()
  
  tmp <- limma.list[[i]] %>%
    .[,sig:=pval_adj<opts$pval & abs(coef)>=opts$coef] %>%
    setorder(pval_adj)

  genes.up <- tmp[coef>0 & sig==T,feature]
  genes.down <- tmp[coef<0 & sig==T,feature]
  background.genes <- tmp[sig==F]
  all.genes <- tmp$feature
  
  # Upregulated genes
  gsea_list[[i]][["up"]] <- enrichGO(
    gene = genes.up,
    universe = all.genes,
    OrgDb       = org.Hs.eg.db,
    keyType     = 'SYMBOL',
    ont         = "ALL", # MF, CC
    pAdjustMethod = "BH",
    qvalueCutoff = 0.10,
    minGSSize = 5,
    maxGSSize = 300
  )
  # View(gsea_list[[i]][["up"]]@result)
    
  # Downregulated genes
  gsea_list[[i]][["down"]] <- enrichGO(
    gene = genes.down,
    universe = all.genes,
    OrgDb       = org.Hs.eg.db,
    keyType     = 'SYMBOL',
    ont         = "BP", # MF, CC
    pAdjustMethod = "BH",
    qvalueCutoff = 0.10,
    minGSSize = 5,
    maxGSSize = 300
  )
      
  return(gsea_list)
  }) %>% unlist(., recursive = FALSE)

# save
saveRDS(gsea_go_list, file.path(io$outdir,"gsea_go_list.rds"))

#############################
## GSEA using MSigDB terms ##
#############################

# i <- contrasts[1]
gsea_msigdb_list <- contrasts %>% future_map(function(i) {
  
  gsea_list <- list()
  
  tmp <- limma.list[[i]] %>%
    .[,sig:=pval_adj<opts$pval & abs(coef)>=opts$coef] %>%
    setorder(pval_adj)
  
  genes.up <- tmp[coef>0 & sig==T,feature]
  genes.down <- tmp[coef<0 & sig==T,feature]
  background.genes <- tmp[sig==F]
  all.genes <- tmp$feature
  
  # Upregulated genes
  gsea_list[[i]][["up"]] <- enricher(
    gene = genes.up,
    universe = all.genes,
    TERM2GENE = msigdbr_hallmark_t2g.dt,
    pvalueCutoff = 0.01,
    pAdjustMethod = "BH",
    minGSSize = 5,
    maxGSSize = 300,
    qvalueCutoff = 0.1
  )
  
  # Downregulated genes
  gsea_list[[i]][["down"]] <- enricher(
    gene = genes.down,
    universe = all.genes,
    TERM2GENE = msigdbr_hallmark_t2g.dt,
    pvalueCutoff = 0.01,
    pAdjustMethod = "BH",
    minGSSize = 5,
    maxGSSize = 300,
    qvalueCutoff = 0.1
  )
  
  return(gsea_list)
}) %>% unlist(., recursive = FALSE)

# save
saveRDS(gsea_msigdb_list, file.path(io$outdir,"gsea_msigdb_list.rds"))


##########################
## GSEA using Mitocarta ##
##########################

# i <- contrasts[1]
gsea_mitocarta_list <- contrasts %>% future_map(function(i) {
  
  gsea_list <- list()
  
  tmp <- limma.list[[i]] %>%
    # .[,sig:=pval_adj<0.01 & abs(coef)>=0.50] %>%
    .[,sig:=pval_adj<opts$pval_threshold] %>%
    setorder(pval_adj)
  
  genes.up <- tmp[coef>0 & sig==T,feature]
  genes.down <- tmp[coef<0 & sig==T,feature]
  background.genes <- tmp[sig==F]
  all.genes <- tmp$feature
  
  # Upregulated genes
  gsea_list[[i]][["up"]] <- enricher(
    gene = genes.up,
    universe = all.genes,
    TERM2GENE = mito_t2g.dt,
    pvalueCutoff = 0.01,
    pAdjustMethod = "BH",
    minGSSize = 5,
    maxGSSize = 300,
    qvalueCutoff = 0.1
  )
  
  # Downregulated genes
  gsea_list[[i]][["down"]] <- enricher(
    gene = genes.down,
    universe = all.genes,
    TERM2GENE = mito_t2g.dt,
    pvalueCutoff = 0.01,
    pAdjustMethod = "BH",
    minGSSize = 5,
    maxGSSize = 300,
    qvalueCutoff = 0.1
  )
  
  return(gsea_list)
}) %>% unlist(., recursive = FALSE)

saveRDS(gsea_mitocarta_list, file.path(io$outdir,"gsea_mitocarta_list.rds"))
       
########################
## GSEA using EnrichR ##
########################

# dbs = listEnrichrDbs()
opts$databases <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023", "Reactome_2022", "MSigDB_Hallmark_2020", "BioCarta_2016", "KEGG_2021_Human")

gsea_enrichr_list <- list()

# k <- opts$databases[1]; i <- contrasts[2]
for (k in opts$databases) {
  print(k)
  
  gsea_enrichr_list[[k]] <- list()
  
  for (i in contrasts) {
    gsea_enrichr_list[[k]][[i]] <- list()
    
    tmp <- limma.list[[i]] %>%
      .[,sig:=pval_adj<opts$pval_threshold] %>%
      setorder(pval_adj)
    
    genes.up <- tmp[coef>0 & sig==T,feature]
    genes.down <- tmp[coef<0 & sig==T,feature]
    background.genes <- tmp[sig==F]
    all.genes <- tmp$feature
    
    # Upregulated genes
    gsea_enrichr_list[[k]][[i]][["up"]] <- enrichr(
      genes = genes.up,
      database = k
    )[[1]] %>% as.data.table %>% .[,c("Term","Overlap","Adjusted.P.value","Odds.Ratio","Genes")] %>%
      setnames(c("gene_signature","overlap","pval_adj","odds_ratio","genes"))
    
    # Downregulated genes
    gsea_enrichr_list[[k]][[i]][["down"]] <- enrichr(
      gene = genes.down,
      database = k
    )[[1]] %>% as.data.table %>% .[,c("Term","Overlap","Adjusted.P.value","Odds.Ratio","Genes")] %>%
      setnames(c("gene_signature","overlap","pval_adj","odds_ratio","genes"))
  }
}

saveRDS(gsea_enrichr_list, file.path(io$outdir,"gsea_enrichr_list.rds"))

#############################
## GSEA using HR3 pathways ##
#############################

# i <- contrasts[3]
gsea_hr3_list <- contrasts %>% future_map(function(i) {
  
  gsea_list <- list()
  
  tmp <- limma.list[[i]] %>%
    # .[,sig:=pval_adj<0.01 & abs(coef)>=0.50] %>%
    .[,sig:=pval_adj<opts$pval_threshold] %>%
    setorder(pval_adj)
  
  genes.up <- tmp[coef>0 & sig==T,feature]
  genes.down <- tmp[coef<0 & sig==T,feature]
  background.genes <- tmp[sig==F]
  all.genes <- tmp$feature
  
  # Upregulated genes
  gsea_list[[i]][["up"]] <- enricher(
    gene = genes.up,
    universe = all.genes,
    TERM2GENE = hr3_t2g.dt,
    pvalueCutoff = 0.01,
    pAdjustMethod = "BH",
    minGSSize = 5,
    maxGSSize = 300,
    qvalueCutoff = 0.1
  )
  
  # Downregulated genes
  gsea_list[[i]][["down"]] <- enricher(
    gene = genes.down,
    universe = all.genes,
    TERM2GENE = hr3_t2g.dt,
    pvalueCutoff = 0.01,
    pAdjustMethod = "BH",
    minGSSize = 5,
    maxGSSize = 300,
    qvalueCutoff = 0.1
  )
  
  return(gsea_list)
}) %>% unlist(., recursive = FALSE)

saveRDS(gsea_hr3_list, file.path(io$outdir,"gsea_hr3_list.rds"))


###################################
## GSEA using decoupleR pathways ##
###################################

# i <- contrasts[2]
gsea_decoupler_list <- contrasts %>% future_map(function(i) {
  
  gsea_list <- list()
  
  tmp <- limma.list[[i]] %>%
    # .[,sig:=pval_adj<0.01 & abs(coef)>=0.50] %>%
    .[,sig:=pval_adj<opts$pval_threshold] %>%
    setorder(pval_adj)
  
  genes.up <- tmp[coef>0 & sig==T,feature]
  genes.down <- tmp[coef<0 & sig==T,feature]
  background.genes <- tmp[sig==F]
  all.genes <- tmp$feature
  
  # Upregulated genes
  gsea_list[[i]][["up"]] <- enricher(
    gene = genes.up,
    universe = all.genes,
    TERM2GENE = decoupler_t2g.dt,
    pvalueCutoff = 0.01,
    pAdjustMethod = "BH",
    minGSSize = 5,
    maxGSSize = 300,
    qvalueCutoff = 0.1
  )
  
  # Downregulated genes
  gsea_list[[i]][["down"]] <- enricher(
    gene = genes.down,
    universe = all.genes,
    TERM2GENE = decoupler_t2g.dt,
    pvalueCutoff = 0.01,
    pAdjustMethod = "BH",
    minGSSize = 5,
    maxGSSize = 300,
    qvalueCutoff = 0.1
  )
  
  return(gsea_list)
}) %>% unlist(., recursive = FALSE)

saveRDS(gsea_decoupler_list, file.path(io$outdir,"gsea_decoupler_list.rds"))



##################################
## GSEA using Reactome pathways ##
##################################

# i <- constrasts[2]
gsea_reactome_list <- contrasts %>% future_map(function(i) {
  
  gsea_list <- list()
  
  tmp <- limma.list[[i]] %>%
    # .[,sig:=pval_adj<0.01 & abs(coef)>=0.50] %>%
    .[,sig:=pval_adj<opts$pval_threshold] %>%
    setorder(pval_adj)
  
  genes.up <- tmp[coef>0 & sig==T,feature]
  genes.down <- tmp[coef<0 & sig==T,feature]
  background.genes <- tmp[sig==F]
  all.genes <- tmp$feature
  
  # Upregulated genes
  gsea_list[[i]][["up"]] <- enricher(
    gene = genes.up,
    universe = all.genes,
    TERM2GENE = reactome_t2g.dt,
    pvalueCutoff = 0.01,
    pAdjustMethod = "BH",
    minGSSize = 5,
    maxGSSize = 300,
    qvalueCutoff = 0.1
  )
  
  # Downregulated genes
  gsea_list[[i]][["down"]] <- enricher(
    gene = genes.down,
    universe = all.genes,
    TERM2GENE = reactome_t2g.dt,
    pvalueCutoff = 0.01,
    pAdjustMethod = "BH",
    minGSSize = 5,
    maxGSSize = 300,
    qvalueCutoff = 0.1
  )
  
  return(gsea_list)
}) %>% unlist(., recursive = FALSE)

saveRDS(gsea_reactome_list, file.path(io$outdir,"gsea_reactome_list.rds"))
