library(msigdbr)

##########################
## Load gene signatures ##
##########################

# KEGG metabolism
io$kegg_metabolism <- ifelse(grepl("^(BAI-|SDI-|CI-)",Sys.info()['nodename']),"/Users/rargelaguet/data/hg38_regulation/gene_lists/kegg_metabgenes.csv","/share/rargelaguet/data/hg38_regulation/gene_lists/kegg_metabgenes.csv")
kegg_metabolism.mtx <- fread(io$kegg_metabolism) %>% as.matrix
kegg_metabolism.mtx[kegg_metabolism.mtx==""] <- as.character(NA)
kegg_metabolism_t2g.dt <- kegg_metabolism.mtx %>% as.data.table(keep.rownames = F) %>%
  data.table::melt(measure.vars=colnames(kegg_metabolism.mtx), variable.name="pathway", value.name="gene", na.rm=T) %>%
  .[,pathway:=as.character(pathway)]

# Mitochondria
io$mitocarta <- ifelse(grepl("^(BAI-|SDI-|CI-)",Sys.info()['nodename']),"/Users/rargelaguet/data/hg38_regulation/mitocarta/MitoPathways3.0.gmx","/share/rargelaguet/data/hg38_regulation/mitocarta/MitoPathways3.0.gmx")
mito_gsea.mtx <- fread(io$mitocarta) %>% .[-1] %>% as.matrix
mito_gsea.mtx[mito_gsea.mtx==""] <- as.character(NA)
mito_t2g.dt <- mito_gsea.mtx %>% as.data.table(keep.rownames = F) %>%
  data.table::melt(measure.vars=colnames(mito_gsea.mtx), variable.name="pathway", value.name="gene", na.rm=T) %>%
  .[,pathway:=as.character(pathway)]

# MsigDB 
# genes_msigdb.dt <- msigdbr(species = "Homo sapiens", category = "C2") %>% as.data.table
# genes_msigdb.dt <- msigdbr(species = "Homo sapiens", category = "C5") %>% as.data.table
# msigdbr_t2g.dt <- genes_msigdb.dt %>% .[,c("gs_name", "human_gene_symbol")] %>% setnames(c("pathway","gene"))

# MsigDB hallmarks
genes_msigdb_hallmark.dt <- msigdbr(species = "Homo sapiens", category = "H") %>% as.data.table
msigdbr_hallmark_t2g.dt <- genes_msigdb_hallmark.dt %>% .[,c("gs_name", "human_gene_symbol")] %>% setnames(c("pathway","gene"))

# HR3 pathways
io$HR3_pathways <- ifelse(grepl("^(BAI-|SDI-|CI-)",Sys.info()['nodename']),"/Users/rargelaguet/data/hg38_regulation/HR3/HR3_t2g_curated-v2.txt","/share/rargelaguet/data/hg38_regulation/HR3/HR3_t2g_curated-v2.txt")
hr3_t2g.dt <- fread(io$HR3_pathways)[,c("gs_name","gene_symbol")] %>% setnames(c("pathway","gene"))

# DecoupleR
io$decoupler_pathways <- ifelse(grepl("^(BAI-|SDI-|CI-)",Sys.info()['nodename']),"/Users/rargelaguet/data/hg38_regulation/decoupleR/decoupler_tf2gene_positive.txt.gz","/share/rargelaguet/data/hg38_regulation/decoupleR/decoupler_tf2gene_positive.txt.gz")
decoupler_t2g.dt <- fread(io$decoupler_pathways)[,c("tf","target")] %>% setnames(c("pathway","gene"))

# Reactome
io$reactome_pathways <- ifelse(grepl("^(BAI-|SDI-|CI-)",Sys.info()['nodename']),"/Users/rargelaguet/data/hg38_regulation/Reactome/reactome_pathway2gene.txt.gz","/share/rargelaguet/data/hg38_regulation/Reactome/reactome_pathway2gene.txt.gz")
reactome_t2g.dt <- fread(io$reactome_pathways)[,c("description","symbol")] %>% setnames(c("pathway","gene"))

#######################
## Prepare gene sets ##
#######################

gene.sets <- list()

# kegg metabolism
for (i in unique(kegg_metabolism_t2g.dt$pathway)) {
  gene.sets[[paste0("kegg#",gsub(" ","_",i))]] <- kegg_metabolism_t2g.dt[pathway==i,gene]
}

# MSigDB
for (i in unique(msigdbr_hallmark_t2g.dt$pathway)) {
  gene.sets[[paste0("msigdb#",gsub(" ","_",i))]] <- msigdbr_hallmark_t2g.dt[pathway==i,gene]
}

# mito
for (i in unique(mito_t2g.dt$pathway)) {
  gene.sets[[paste0("mito#",i)]] <- mito_t2g.dt[pathway==i,gene]
}

# HR3
for (i in unique(hr3_t2g.dt$pathway)) {
  gene.sets[[paste0("hr3#",i)]] <- hr3_t2g.dt[pathway==i,gene]
}

# Decoupler
for (i in unique(decoupler_t2g.dt$pathway)) {
  gene.sets[[paste0("decoupler#",i)]] <- decoupler_t2g.dt[pathway==i,gene]
}

# Reactome
for (i in unique(reactome_t2g.dt$pathway)) {
  gene.sets[[paste0("reactome#",i)]] <- reactome_t2g.dt[pathway==i,gene]
}

# filter gene sets by size
tmp <- sapply(gene.sets,length)
gene.sets <- gene.sets[tmp>=10 & tmp<1000]
