source(here::here("settings.R"))
source(here::here("utils.R"))

library(clusterProfiler)
library(org.Hs.eg.db)
# library(enrichplot)

#####################
## Define settings ##
#####################

# I/O
io$diff_results <- file.path(io$basedir,"results/rna/differential/metacells/alias/parsed/diff_expr_results.txt.gz")
io$outdir <- file.path(io$basedir,"results/rna/gene_signatures/differential"); dir.create(io$outdir, showWarnings = F)
# io$outdir.scatterplots <- file.path(io$basedir,"results/rna/differential/metacells/alias/pdf/scatterplots"); dir.create(io$outdir.scatterplots, showWarnings = F)

# Options
opts$groups <- opts$aliases

opts$group_variable <- "alias"

opts$conditions <- c(
  "8OFF",
  "6OFF-2ON",
  "4OFF-4ON",
  "2OFF-2ON-4OFF",
  "4ON-4OFF",
  "8ON"
)

opts$cell_lines <- c(
  "4F",
  "2F"
)

library(furrr)
plan(multisession, workers = 6)

########################
## Load gene metadata ##
########################

TEs <- fread(io$TEs, header = F)[[1]]

TFs <- fread(io$TFs)[[1]]

##########################################
## Load differential expression results ##
##########################################

diff.dt <- fread(io$diff_results) %>%
  setnames("gene","feature") %>%
  .[,sig:=padj_fdr<0.01 & abs(logFC)>=2]

diff.dt %>% 
  .[,transcript_type:=c("gene","te")[as.numeric(feature%in%TEs)+1]] %>%
  .[feature%in%TFs,transcript_type:="tf"]

# diff.dt[,cell_lineA:=strsplit(groupA,"_") %>% map_chr(1)] %>% .[,cell_lineA:=factor(cell_lineA,levels=opts$cell_lines)]
# diff.dt[,cell_lineB:=strsplit(groupB,"_") %>% map_chr(1)] %>% .[,cell_lineB:=factor(cell_lineB,levels=opts$cell_lines)]
# 
# diff.dt[,conditionA:=strsplit(groupA,"_") %>% map_chr(2)] %>% .[,conditionA:=factor(conditionA,levels=opts$conditions)]
# diff.dt[,conditionB:=strsplit(groupB,"_") %>% map_chr(2)] %>% .[,conditionB:=factor(conditionB,levels=opts$conditions)]

##########################
## Load gene signatures ##
##########################

source(here::here("rna/gene_signatures/load_gene_signatures.R"))

#########################
## GSEA using GO terms ##
#########################

gsea_go_list <- opts$groups %>% future_map(function(i) {
  
  gsea_list <- list()
  
  for (j in opts$groups) {
    if (i!=j) {
      gsea_list[[sprintf("%s_vs_%s",i,j)]] <- list()
      
      diff_filt.dt <- diff.dt[groupA==i & groupB==j] %>% .[transcript_type!="te"]
      
      genes.up <- diff_filt.dt[logFC>1 & sig==T,feature]
      genes.down <- diff_filt.dt[logFC<(-1) & sig==T,feature]
      background.genes <- diff_filt.dt[sig==F]
      all.genes <- diff_filt.dt$feature
      
      # Upregulated genes
      gsea_list[[sprintf("%s_vs_%s",i,j)]][["up"]] <- enrichGO(
        gene = genes.up,
        universe = all.genes,
        OrgDb       = org.Hs.eg.db,
        keyType     = 'SYMBOL',
        ont         = "BP", # MF, CC
        pAdjustMethod = "BH",
        qvalueCutoff = 0.10,
        minGSSize = 10,
        maxGSSize = 500
      )
      
      # Downregulated genes
      gsea_list[[sprintf("%s_vs_%s",i,j)]][["down"]] <- enrichGO(
        gene = genes.down,
        universe = all.genes,
        OrgDb       = org.Hs.eg.db,
        keyType     = 'SYMBOL',
        ont         = "BP", # MF, CC
        pAdjustMethod = "BH",
        qvalueCutoff = 0.10,
        minGSSize = 10,
        maxGSSize = 500
      )
      
    }
  }
  return(gsea_list)
}) %>% unlist(., recursive = FALSE)

# save
saveRDS(gsea_go_list, file.path(io$outdir,"gsea_go_list.rds"))


#############################
## GSEA using MSigDB terms ##
#############################

# i <- opts$groups[1]; j <- opts$groups[2]

# gsea_msigdb_list <- opts$groups %>% future_map(function(i) {
  
#   gsea_list <- list()
  
#   for (j in opts$groups) {
#     if (i!=j) {
#       gsea_list[[sprintf("%s_vs_%s",i,j)]] <- list()
      
#       diff_filt.dt <- diff.dt[groupA==i & groupB==j] %>% .[transcript_type!="te"]
      
#       genes.up <- diff_filt.dt[logFC>1 & sig==T,feature]
#       genes.down <- diff_filt.dt[logFC<(-1) & sig==T,feature]
#       background.genes <- diff_filt.dt[sig==F]
#       all.genes <- diff_filt.dt$feature
  
#       # Upregulated genes      
#       gsea_list[[sprintf("%s_vs_%s",i,j)]][["up"]] <- enricher(
#         gene = genes.up,
#         universe = all.genes,
#         TERM2GENE = msigdbr_hallmark_t2g.dt,
#         pvalueCutoff = 0.01,
#         pAdjustMethod = "BH",
#         minGSSize = 10,
#         maxGSSize = 500,
#         qvalueCutoff = 0.1
#       )
      
#       # Downregulated genes      
#       gsea_list[[sprintf("%s_vs_%s",i,j)]][["down"]] <- enricher(
#         gene = genes.down,
#         universe = all.genes,
#         TERM2GENE = msigdbr_hallmark_t2g.dt,
#         pvalueCutoff = 0.01,
#         pAdjustMethod = "BH",
#         minGSSize = 10,
#         maxGSSize = 500,
#         qvalueCutoff = 0.1
#       )
#     }
#   }
#   return(gsea_list)
# }) %>% unlist(., recursive = FALSE)

# # save
# saveRDS(gsea_msigdb_list, file.path(io$outdir,"gsea_mgsidb_list.rds"))
       

#############
## EnrichR ##
#############

###GeneOntology
#DE_2F = DE genes from wilcox test

dbs = listEnrichrDbs()
dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023")
enr_results <- enrichr(gene = head(rownames(DE_2F),500), database = dbs)

p1 = plotEnrich(enr_results[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
