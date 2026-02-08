source(here::here("settings.R"))
source(here::here("utils.R"))

library(enrichR)

#####################
## Define settings ##
#####################

# I/O
io$diff_results <- file.path(io$basedir,"results/rna/differential/metacells/sample/parsed/diff_expr_results.txt.gz")
io$outdir <- file.path(io$basedir,"results/rna/gene_signatures/GSEA_differential/enrichR"); dir.create(io$outdir, showWarnings = F, recursive=T)

# Options
opts$group_variable <- "sample"

opts$groups <- c( 
  "OS_8OFF",
  "OS_6OFF-2ON",
  "OS_4OFF-4ON",
  "OS_8ON",
  "OS_2OFF-2ON-4OFF",
  "OS_4ON-4OFF",
  "OSK_8OFF",
  "OSK_4OFF-4ON",
  "OSK_8ON",
  "OSK_4ON-4OFF",
  "OSKM_8OFF",
  "OSKM_6OFF-2ON",
  "OSKM_4OFF-4ON",
  "OSKM_8ON",
  "OSKM_2OFF-2ON-4OFF",
  "OSKM_4ON-4OFF"
)

opts$databases <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023", "Reactome_2022", "MSigDB_Hallmark_2020", "BioCarta_2016", "KEGG_2021_Human")

library(furrr)
plan(multisession, workers = 6)

########################
## Load gene metadata ##
########################

TEs <- fread(io$TEs, header = F)[[1]]

TFs <- fread(io$TFs)[[1]]

surface_markers <- unique(fread(io$cell_surface_markers)[[1]])

##########################################
## Load differential expression results ##
##########################################

diff.dt <- fread(io$diff_results) %>%
  setnames("gene","feature") %>%
  .[,sig:=padj_fdr<0.01 & abs(logFC)>=1] %>%
  .[,transcript_type:="Genes"] %>%
  .[feature%in%TEs,transcript_type:="Transposable Elements"] %>%
  .[feature%in%TFs,transcript_type:="Transcription Factors"] %>%
  .[feature%in%surface_markers,transcript_type:="Surface markers"] %>%
  .[,transcript_type:=factor(transcript_type, levels=c("Genes","Transcription Factors","Surface markers","Transposable Elements"))]

# diff.dt[,cell_lineA:=strsplit(groupA,"_") %>% map_chr(1)] %>% .[,cell_lineA:=factor(cell_lineA,levels=opts$cell_lines)]
# diff.dt[,cell_lineB:=strsplit(groupB,"_") %>% map_chr(1)] %>% .[,cell_lineB:=factor(cell_lineB,levels=opts$cell_lines)]
# 
# diff.dt[,conditionA:=strsplit(groupA,"_") %>% map_chr(2)] %>% .[,conditionA:=factor(conditionA,levels=opts$conditions)]
# diff.dt[,conditionB:=strsplit(groupB,"_") %>% map_chr(2)] %>% .[,conditionB:=factor(conditionB,levels=opts$conditions)]

diff.dt <- diff.dt[!transcript_type%in%c("Transposable Elements")]
diff.dt <- diff.dt[!feature%in%c("exogenous_UTR")]

##########################
## Load gene signatures ##
##########################

source(here::here("rna/gene_signatures/load_gene_signatures.R"))

########################
## GSEA using enrichR ##
########################

# dbs = listEnrichrDbs()


# opts$groups  <- opts$groups[1:2]
gsea_list <- list()

# i <- opts$groups[1]; j <- opts$groups[2]; k <- opts$databases[1]
for (k in opts$databases) {
  print(k)
  
  gsea_list[[k]] <- list()

  for (i in opts$groups) {
    for (j in opts$groups) {
      if (i!=j) {
        gsea_list[[k]][[sprintf("%s_vs_%s",i,j)]] <- list()
        
        diff_filt.dt <- diff.dt[groupA==i & groupB==j]# %>% .[transcript_type!="te"]
        
        genes.up <- diff_filt.dt[logFC>0 & sig==T,feature]
        genes.down <- diff_filt.dt[logFC<0 & sig==T,feature]
        # background.genes <- diff_filt.dt[sig==F]
        # all.genes <- diff_filt.dt$feature
        
        # Upregulated genes
        gsea_list[[k]][[sprintf("%s_vs_%s",i,j)]][["up"]] <- enrichr(
          genes = genes.up,
          database = k
        )[[1]] %>% as.data.table %>% .[,c("Term","Overlap","Adjusted.P.value","Odds.Ratio","Genes")] %>%
          setnames(c("gene_signature","overlap","pval_adj","odds_ratio","genes"))
        
        # Downregulated genes
        gsea_list[[k]][[sprintf("%s_vs_%s",i,j)]][["down"]] <- enrichr(
          gene = genes.down,
          database = k
        )[[1]] %>% as.data.table %>% .[,c("Term","Overlap","Adjusted.P.value","Odds.Ratio","Genes")] %>%
          setnames(c("gene_signature","overlap","pval_adj","odds_ratio","genes"))
      }
    }
  }
}

# save
saveRDS(gsea_list, file.path(io$outdir,"gsea_enrichr_list.rds"))

#############
## TEST ##
#############

# head(gsea_list$GO_Molecular_Function_2023$`OS_8OFF_vs_OS_6OFF-2ON`[["up"]][,c(1,2,3)])
# head(gsea_list$GO_Molecular_Function_2023$`OS_8OFF_vs_OS_6OFF-2ON`[["down"]][,c(1,2,3)])
# 
# head(gsea_list$GO_Biological_Process_2023$`OS_8OFF_vs_OS_6OFF-2ON`[["up"]][,c(1,2,3)])
# head(gsea_list$GO_Biological_Process_2023$`OS_8OFF_vs_OS_6OFF-2ON`[["down"]][,c(1,2,3)])
# 
# head(gsea_list$GO_Cellular_Component_2023$`OS_8OFF_vs_OS_6OFF-2ON`[["up"]][,c(1,2,3)])
# head(gsea_list$GO_Cellular_Component_2023$`OS_8OFF_vs_OS_6OFF-2ON`[["down"]][,c(1,2,3)])
# 
# head(gsea_list$MSigDB_Hallmark_2020$`OS_8OFF_vs_OS_6OFF-2ON`[["up"]][,c(1,2,3)])
# head(gsea_list$MSigDB_Hallmark_2020$`OS_8OFF_vs_OS_6OFF-2ON`[["down"]][,c(1,2,3)])