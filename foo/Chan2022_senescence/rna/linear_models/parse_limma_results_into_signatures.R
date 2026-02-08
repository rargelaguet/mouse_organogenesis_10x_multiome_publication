source(here::here("settings.R"))
source(here::here("utils.R"))

# TO DO: SEPARATE QUIESCENCE VS SENESCENCE SIGNATURE
# TO DO: GET RID OF CELL CYCLE EFFECT

#####################
## Define settings ##
#####################

# I/O
io$limma_results_cell_density <- file.path(io$basedir,"rna/results/limma/limma_results_cell_density.txt.gz")
io$limma_results_hTERT <- file.path(io$basedir,"rna/results/limma/limma_results_hTERT.txt.gz")
io$limma_results_radiation_senescence <- file.path(io$basedir,"rna/results/limma/limma_results_radiation_senescence.txt.gz")
io$limma_results_replicative_senescence <- file.path(io$basedir,"rna/results/limma/limma_results_replicative_senescence.txt.gz")
io$outdir <- file.path(io$basedir,"rna/results/limma/signatures"); dir.create(io$outdir, showWarnings = F, recursive = T)

# options
opts$pval <- 1e-4
opts$coef <- 0.1

########################
## Load limma results ##
########################

limma_list <- list()
limma_list[["cell_density"]] <- fread(io$limma_results_cell_density) %>% .[is.na(symbol),symbol:=ens_id] %>% setnames("symbol","feature") %>% setorder(pval)
limma_list[["hTERT"]] <- fread(io$limma_results_hTERT) %>% .[is.na(symbol),symbol:=ens_id] %>% setnames("symbol","feature") %>% setorder(pval)
limma_list[["radiation_senescence"]] <- fread(io$limma_results_radiation_senescence) %>% .[is.na(symbol),symbol:=ens_id] %>% setnames("symbol","feature") %>% setorder(pval)
limma_list[["replicative_senescence"]] <- fread(io$limma_results_replicative_senescence) %>% .[is.na(symbol),symbol:=ens_id] %>% setnames("symbol","feature") %>% setorder(pval)

##########################
## Load gene signatures ##
##########################

gene_sets <- list()
gene_sets[["Transcription Factors"]] <- fread(io$TFs, header=F)[[1]]
gene_sets[["Surface markers"]] <- unique(fread(io$cell_surface_markers, header = F)[[1]])
gene_sets[["Epigenetic factors"]] <- unique(fread(io$epigenetic_factors, header=F)[[1]])
gene_sets[["Ribosomal"]] <-  fread(io$ribo_genes, header=F)[[1]]
gene_sets[["Mitochondrial"]] <- fread(io$mito_genes, header=F)[[1]]
gene_sets[["DNA repair"]] <- fread(io$dna_repair_genes, header=F)[[1]]
gene_sets[["Cell cycle"]] <- fread(io$cell_cycle_g2m_genes, header=F)[[1]]
gene_sets[["ECM"]] <- fread(io$ecm_genes, header=F)[[1]]


###########
## Parse ##
###########

contrasts <- names(limma_list)

tmp <- list()
# i <- contrasts[1]
for (i in contrasts) {
  
  tmp[[i]] <- copy(limma_list[[i]]) %>%
    .[is.na(coef),c("coef","pval_adj"):=list(0,1)] %>%
    .[,weight:=-log10(pval_adj+1e-100)] %>% .[,weight:=round(weight,3)] %>%
    .[,direction:=ifelse(coef>0,"Up","Down")] %>%
    .[,sig:=pval_adj<opts$pval & abs(coef)>=opts$coef] %>% .[sig==TRUE] %>%
    # .[pval_adj<opts$pval_max] %>%
    .[,signature:=sprintf("%s-%s",i,direction)] %>%
    .[,c("coef","pval_adj"):=NULL] %>%
    setorder(direction,-weight) %>%
    .[,c("signature","feature","direction","weight")] %>%
    .[,condition:=i]
    
}

#############################
## Save general signatures ##
#############################

signature.dt <- rbindlist(tmp)
signature.dt[,.N,by=c("condition","direction")]

fwrite(signature.dt, file.path(io$outdir,"Chan2022_senescence_signatures.txt.gz"), quote=F, na="NA",sep="\t")


##############################
## Save gene set signatures ##
##############################

outdir <- file.path(io$outdir,"per_gene_set"); dir.create(outdir, showWarnings = F)

gene.sets.to.save <- names(gene_sets)

for (i in gene.sets.to.save) {
  
  signature_i.dt <- signature.dt[feature%in%gene_sets[[i]]] %>%
    .[,gene_set:=i]
  
  print(signature_i.dt[,.N,by=c("condition","direction","gene_set")])
  
  fwrite(signature.dt, file.path(outdir,sprintf("Chan2022_senescence_signatures_%s.txt.gz",i)), quote=F, na="NA",sep="\t")
}

###################################
## Senescence-specific signature ##
###################################



