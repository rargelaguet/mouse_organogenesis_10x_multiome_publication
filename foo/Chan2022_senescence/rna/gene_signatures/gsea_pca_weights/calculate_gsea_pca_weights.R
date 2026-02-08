source(here::here("settings.R"))
source(here::here("utils.R"))

library(clusterProfiler)
library(org.Hs.eg.db)
# library(enrichplot)
library(enrichR)
library(msigdbr)

#####################
## Define settings ##
#####################

# I/O
io$pca_weights <- file.path(io$basedir,"results/rna/dimensionality_reduction/PCA_weights.txt.gz")
io$outdir <- file.path(io$basedir,"results/rna/gene_signatures/pca_weights"); dir.create(io$outdir, showWarnings = F, recursive = T)

# Options
library(furrr)
# plan(multisession, workers = 3)

########################
## Load gene metadata ##
########################

TEs <- fread(io$TEs, header = F)[[1]]

######################
## Load PCA weights ##
######################

pca_weights.dt <- fread(io$pca_weights) %>%
  .[!gene%in%TEs]
  
pca_weights.dt[,PC1:=-PC1]
pca_weights.dt[,PC2:=-PC2]

##########################
## Load gene signatures ##
##########################

source(here::here("rna/gene_signatures/load_gene_signatures.R"))

#########################
## GSEA using GO terms ##
#########################

PCs.to.use <- c("PC1","PC2")

# i <- PCs.to.use[1]
gsea_go_list <- PCs.to.use %>% future_map(function(i) {
  
  gsea_list <- list()
  
  tmp <- pca_weights.dt[,c("gene",i),with=F] %>% 
    setnames(c("feature","weight")) %>% 
    setorder(-weight)

  genes.up <- tmp[weight>0,feature] %>% head(n=250)
  genes.down <- tmp[weight<0,feature] %>% tail(n=250)
  all.genes <- tmp$feature
  background.genes <- all.genes
  # background.genes <- tmp[abs(weight)<0.01]
  
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
    
  # Downregulated genes
  gsea_list[[i]][["down"]] <- enrichGO(
    gene = genes.down,
    universe = all.genes,
    OrgDb       = org.Hs.eg.db,
    keyType     = 'SYMBOL',
    ont         = "ALL", # MF, CC
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

genes_msigdb_hallmark.dt <- msigdbr(species = "Homo sapiens", category = "H") %>% as.data.table
msigdbr_hallmark_t2g.dt <- genes_msigdb_hallmark.dt %>% .[,c("gs_name", "human_gene_symbol")] %>% setnames(c("pathway","gene"))

PCs.to.use <- c("PC1","PC2")

# i <- PCs.to.use[1]
gsea_msigdb_list <- PCs.to.use %>% future_map(function(i) {
  
  gsea_list <- list()
  
  tmp <- pca_weights.dt[,c("gene",i),with=F] %>% 
    setnames(c("feature","weight")) %>% 
    setorder(-weight)
  
  genes.up <- tmp[weight>0,feature] %>% head(n=150)
  genes.down <- tmp[weight<0,feature] %>% tail(n=150)
  background.genes <- tmp[abs(weight)<0.01]
  all.genes <- tmp$feature
  
  # Upregulated genes
  gsea_list[[i]][["up"]] <- enricher(
    gene = genes.up,
    universe = all.genes,
    TERM2GENE = msigdbr_hallmark_t2g.dt,
    pvalueCutoff = 0.01,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff = 0.1
  )
  
  # Downregulated genes
  gsea_list[[i]][["down"]] <- enricher(
    gene = genes.down,
    universe = all.genes,
    TERM2GENE = msigdbr_hallmark_t2g.dt,
    pvalueCutoff = 0.01,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff = 0.1
  )
  
  return(gsea_list)
}) %>% unlist(., recursive = FALSE)

# save
saveRDS(gsea_msigdb_list, file.path(io$outdir,"gsea_mgsidb_list.rds"))
       

################################
## GSEA using MitoCarta terms ##
################################

PCs.to.use <- c("PC1","PC2")

# i <- PCs.to.use[1]
gsea_mitocarta_list <- PCs.to.use %>% future_map(function(i) {
  
  gsea_list <- list()
  
  tmp <- pca_weights.dt[,c("gene",i),with=F] %>% 
    setnames(c("feature","weight")) %>% 
    setorder(-weight)
  
  genes.up <- tmp[weight>0,feature] %>% head(n=150)
  genes.down <- tmp[weight<0,feature] %>% tail(n=150)
  background.genes <- tmp[abs(weight)<0.01]
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
## GSEA using enrichR ##
########################

# dbs = listEnrichrDbs()
opts$databases <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023", "Reactome_2022", "MSigDB_Hallmark_2020", "BioCarta_2016", "KEGG_2021_Human")

PCs.to.use <- c("PC1","PC2")

gsea_list <- list()

# i <- PCs.to.use[1]; k <- opts$databases[1]
for (k in opts$databases) {
  print(k)
  
  gsea_list[[k]] <- list()
  
  for (i in PCs.to.use) {
      gsea_list[[k]][[i]] <- list()
      
      tmp <- pca_weights.dt[,c("gene",i),with=F] %>% 
        setnames(c("feature","weight")) %>% 
        setorder(-weight)
      
      genes.up <- tmp[weight>0,feature] %>% head(n=150)
      genes.down <- tmp[weight<0,feature] %>% tail(n=150)
      background.genes <- tmp[abs(weight)<0.01]
      all.genes <- tmp$feature
      
      # Upregulated genes
      gsea_list[[k]][[i]][["up"]] <- enrichr(
        genes = genes.up,
        database = k
      )[[1]] %>% as.data.table %>% .[,c("Term","Overlap","Adjusted.P.value","Odds.Ratio","Genes")] %>%
        setnames(c("gene_signature","overlap","pval_adj","odds_ratio","genes"))
      
      # Downregulated genes
      gsea_list[[k]][[i]][["down"]] <- enrichr(
        gene = genes.down,
        database = k
      )[[1]] %>% as.data.table %>% .[,c("Term","Overlap","Adjusted.P.value","Odds.Ratio","Genes")] %>%
        setnames(c("gene_signature","overlap","pval_adj","odds_ratio","genes"))
    }
  }

# save
saveRDS(gsea_list, file.path(io$outdir,"gsea_enrichr_list_top_150_genes.rds"))
