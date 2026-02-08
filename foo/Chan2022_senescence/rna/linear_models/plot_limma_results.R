source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$limma_results_cell_density <- file.path(io$basedir,"rna/results/limma/limma_results_cell_density.txt.gz")
io$limma_results_hTERT <- file.path(io$basedir,"rna/results/limma/limma_results_hTERT.txt.gz")
io$limma_results_radiation_senescence <- file.path(io$basedir,"rna/results/limma/limma_results_radiation_senescence.txt.gz")
io$limma_results_replicative_senescence <- file.path(io$basedir,"rna/results/limma/limma_results_replicative_senescence.txt.gz")
io$outdir <- file.path(io$basedir,"rna/results/limma/pdf"); dir.create(io$outdir, showWarnings = F, recursive = T)

# Options
opts$pval <- 1e-4
opts$coef <- 0.05

########################
## Load limma results ##
########################

limma_list <- list()
limma_list[["cell_density"]] <- fread(io$limma_results_cell_density) %>% .[is.na(symbol),symbol:=ens_id] %>% setnames("symbol","feature") %>% setorder(pval)
limma_list[["hTERT"]] <- fread(io$limma_results_hTERT) %>% .[is.na(symbol),symbol:=ens_id] %>% setnames("symbol","feature") %>% setorder(pval)
limma_list[["radiation_senescence"]] <- fread(io$limma_results_radiation_senescence) %>% .[is.na(symbol),symbol:=ens_id] %>% setnames("symbol","feature") %>% setorder(pval)
limma_list[["replicative_senescence"]] <- fread(io$limma_results_replicative_senescence) %>% .[is.na(symbol),symbol:=ens_id] %>% setnames("symbol","feature") %>% setorder(pval)

########################
## Load gene metadata ##
########################

gene_metadata.dt <- fread(io$gene_metadata)

gene_sets <- list()
gene_sets[["Transcription Factors"]] <- fread(io$TFs, header=F)[[1]]
gene_sets[["Surface markers"]] <- unique(fread(io$cell_surface_markers, header = F)[[1]])
gene_sets[["Epigenetic factors"]] <- unique(fread(io$epigenetic_factors, header=F)[[1]])
gene_sets[["Ribosomal"]] <-  fread(io$ribo_genes, header=F)[[1]]
gene_sets[["Mitochondrial"]] <- fread(io$mito_genes, header=F)[[1]]
gene_sets[["DNA repair"]] <- fread(io$dna_repair_genes, header=F)[[1]]
gene_sets[["Cell cycle"]] <- fread(io$cell_cycle_g2m_genes, header=F)[[1]]
gene_sets[["ECM"]] <- fread(io$ecm_genes, header=F)[[1]]

#####################
## Filter features ##
#####################

# limma.dt <- limma.dt[!grepl("^RPL|^MT",feature)]

##################
## Volcano plot ##
##################

conditions.to.plot <- names(limma_list)

genes.to.highlight <- c()

# i <- conditions.to.plot[1]
for (i in conditions.to.plot) {
  
  to.plot <- copy(limma_list[[i]]) %>%
    # .[,pval_adj:=NULL] %>% setnames("pval","pval_adj") %>% # TEMPORARY
    .[is.na(coef),c("coef","pval_adj"):=list(0,1)] %>%
    .[,sig:=pval_adj<opts$pval & abs(coef)>=opts$coef] %>%
    .[,log_pval:=-log10(pval_adj+1e-100)] %>%
    .[,transcript_type:="Genes"] %>%
    .[feature%in%gene_sets[["Surface markers"]],transcript_type:="Surface markers"] %>%
    .[feature%in%gene_sets[["Epigenetic factors"]],transcript_type:="Epigenetic factors"] %>%
    .[feature%in%gene_sets[["Ribosomal"]],transcript_type:="Ribosomal"] %>%
    .[feature%in%gene_sets[["Mitochondrial"]],transcript_type:="Mitochondrial"] %>%
    .[feature%in%gene_sets[["DNA repair"]],transcript_type:="DNA repair"] %>%
    .[feature%in%gene_sets[["Cell cycle"]],transcript_type:="Cell cycle"] %>%
    .[feature%in%gene_sets[["ECM"]],transcript_type:="ECM"] %>%
    .[feature%in%gene_sets[["Transcription Factors"]],transcript_type:="Transcription Factors"] %>%
    .[,transcript_type:=factor(transcript_type, levels=c("Genes",names(gene_sets)))]
    
    negative_hits <- to.plot[sig==TRUE & coef<0,feature]
    positive_hits <- to.plot[sig==TRUE & coef>0,feature]
    all <- nrow(to.plot)
    
    to.plot <- to.plot[abs(coef)>0.001]
    
    xlim <- max(c(0.05,abs(to.plot$coef)), na.rm=T)
    # xlim <- 0.051
    ylim <- max(c(5,to.plot$log_pval), na.rm=T)
    
    tmp <- to.plot[,.(N=as.character(.N), Npos=as.character(sum(sig==TRUE & coef>0)), Nneg=as.character(sum(sig==TRUE & coef<0))),by="transcript_type"]
    tmp %>% .[,Npos:=sprintf("%s (+)",Npos)] %>% .[,Nneg:=sprintf("%s (-)",Nneg)]
    
    p <- ggplot(to.plot, aes(x=coef, y=log_pval)) +
      ggrastr::geom_point_rast(aes(color=sig, size=sig), data=to.plot[!feature%in%genes.to.highlight], raster.dpi=300) +
      geom_point(size=0.65, color="purple", data=to.plot[feature%in%genes.to.highlight]) +
      geom_segment(x=0, xend=0, y=0, yend=ylim, color="orange", linewidth=0.1) +
      scale_color_manual(values=c("gray30","red")) +
      scale_size_manual(values=c(0.20,0.30)) +
      scale_x_continuous(limits=c(-xlim,xlim)) +
      scale_y_continuous(limits=c(0,ylim+0.5)) +
      labs(x="limma coefficient", y=expression(paste("-log"[10],"(p.value)")), title="") +
      geom_text_repel(data=to.plot[feature%in%genes.to.highlight], aes(x=coef, y=log_pval, label=feature), color="purple", size=3.5, segment.size=0.1, max.overlaps = Inf) +
      geom_text(aes(label=N), x=0, y=ylim+0.5, size=3, data=tmp) +
      geom_text(aes(label=Nneg), x=-xlim+(xlim/3), y=ylim+0.5, size=3, data=tmp) +
      geom_text(aes(label=Npos), x=xlim-(xlim/3), y=ylim+0.5, size=3, data=tmp) +
      annotate("text", x=-xlim+xlim/3, y=0, size=3.5, label="Down") +
      annotate("text", x=xlim-xlim/3, y=0, size=3.5, label="Up") +
      facet_wrap(~transcript_type, nrow=3) +
      theme_classic() +
      theme(
        strip.text = element_text(size=rel(1), color='black'),
        axis.text = element_text(size=rel(0.90), color='black'),
        axis.title = element_text(size=rel(1.15), color='black'),
        legend.position="none"
      )
    
    for (j in unique(to.plot$transcript_type)) {
      p <- p + geom_text_repel(data=head(to.plot[transcript_type==j & coef>0 & sig==TRUE],n=10), aes(x=coef, y=log_pval, label=feature), segment.size=0.05, size=2.50, max.overlaps = Inf)
      p <- p + geom_text_repel(data=head(to.plot[transcript_type==j & coef<0 & sig==TRUE],n=10), aes(x=coef, y=log_pval, label=feature), segment.size=0.05, size=2.50, max.overlaps = Inf)
    }
    
  pdf(file.path(io$outdir,sprintf("limma_expr_%s_volcano_plot_chan2022.pdf",i)), width=12, height=10)
  print(p)
  dev.off()
}


################################
## Volcano plots per gene set ##
################################

outdir <- file.path(io$basedir,"rna/results/limma/pdf/per_gene_set"); dir.create(outdir, showWarnings = F, recursive = T)

conditions.to.plot <- names(limma_list)

# i <- conditions.to.plot[1]; j <- names(gene_sets)[1]
for (i in conditions.to.plot) {
  for (j in names(gene_sets)) {
    
    to.plot <- limma_list[[i]] %>%
      .[feature%in%gene_sets[[j]]] %>%
      .[is.na(coef),c("coef","pval_adj"):=list(0,1)] %>%
      setorder(pval_adj) %>%
      .[,sig:=pval_adj<opts$pval & abs(coef)>=opts$coef] %>%
      .[,log_pval:=-log10(pval_adj+1e-100)]
    
    negative_hits <- to.plot[sig==TRUE & coef<0,feature]
    positive_hits <- to.plot[sig==TRUE & coef>0,feature]
    all <- nrow(to.plot)
    
    to.plot <- to.plot[abs(coef)>0.001]
    
    xlim <- max(c(0.05,abs(to.plot$coef)), na.rm=T)
    ylim <- max(c(5,to.plot$log_pval), na.rm=T)
    
    tmp <- to.plot[,.(N=as.character(.N), Npos=as.character(sum(sig==TRUE & coef>0)), Nneg=as.character(sum(sig==TRUE & coef<0)))]
    tmp %>% .[,Npos:=sprintf("%s (+)",Npos)] %>% .[,Nneg:=sprintf("%s (-)",Nneg)]
    
    p <- ggplot(to.plot, aes(x=coef, y=log_pval)) +
      geom_point(aes(color=sig, size=sig)) +
      geom_segment(x=0, xend=0, y=0, yend=ylim, color="orange", linewidth=0.1) +
      scale_color_manual(values=c("black","red")) +
      scale_size_manual(values=c(0.50,1)) +
      scale_x_continuous(limits=c(-xlim,xlim)) +
      scale_y_continuous(limits=c(0,ylim+1.5)) +
      labs(x="limma coefficient", y=expression(paste("-log"[10],"(p.value)")), title=i) +
      geom_text(aes(label=N), x=0, y=ylim+1, size=4, data=tmp) +
      geom_text(aes(label=Nneg), x=-xlim+(xlim/3), y=ylim+1, size=4, data=tmp) +
      geom_text(aes(label=Npos), x=xlim-(xlim/3), y=ylim+1, size=4, data=tmp) +
      annotate("text", x=-xlim+xlim/3, y=0, size=3, label="Up in young") +
      annotate("text", x=xlim-xlim/3, y=0, size=3, label="Up in old") +
      geom_text_repel(data=head(to.plot[sig==TRUE & coef>0],n=15), aes(x=coef, y=log_pval, label=feature), segment.size=0.05, size=2.75, max.overlaps = Inf) +
      geom_text_repel(data=head(to.plot[sig==TRUE & coef<0],n=15), aes(x=coef, y=log_pval, label=feature), segment.size=0.05, size=2.75, max.overlaps = Inf) +
      theme_classic() +
      theme(
        plot.title = element_text(hjust=0.5),
        strip.text = element_text(size=rel(1.05), color='black'),
        axis.text = element_text(size=rel(0.75), color='black'),
        axis.title = element_text(size=rel(1.0), color='black'),
        legend.position="none"
      )
    
    pdf(file.path(outdir,sprintf("%s_%s_limma_rna_volcano_plot_chan2022.pdf",i,j)), width=7, height=6)
    print(p)
    dev.off()
  }
}


###############################
## Volcano plots per pathway ##
###############################

outdir <- file.path(io$basedir,"rna/results/limma/pdf/per_pathway"); dir.create(outdir, showWarnings = F, recursive = T)

source(here::here("rna/gene_signatures/load_gene_signatures.R"))

# yamanaka_regulons <- readRDS("/Users/rargelaguet/data/Rand2021_oskm_combinations_microarray/results/regulons/yamanaka_regulons_100_Rand2017.rds")
# gene_sets <- yamanaka_regulons
# names(gene_sets) <- paste0(names(gene_sets),"_regulon")


gene.sets.to.plot <- gene.sets[grep("Senescence",names(gene.sets),value=T,ignore.case=T)]
gene.sets.to.plot <- gene.sets[c("hr3#ISR","hr3#UPR","reactome#Unfolded Protein Response (UPR)")]
gene.sets.to.plot[[c("hr3#ISR")]] <- c(gene.sets.to.plot[[c("hr3#ISR")]],"EIF2B1", "EIF2B2", "EIF2B3", "EIF2B4", "EIF2B5", "EEF2", "EEF2K")

# i <- conditions.to.plot[1]; j <- names(gene.sets.to.plot)[[1]]
for (i in conditions.to.plot) {
  for (j in names(gene.sets.to.plot)) {
    
    to.plot <- limma_list[[i]] %>%
      .[feature%in%gene.sets.to.plot[[j]]] %>%
      .[is.na(coef),c("coef","pval_adj"):=list(0,1)] %>%
      .[,sig:=pval_adj<opts$pval & abs(coef)>=opts$coef] %>%
      .[,log_pval:=-log10(pval_adj+1e-100)]
    
    negative_hits <- to.plot[sig==TRUE & coef<0,feature]
    positive_hits <- to.plot[sig==TRUE & coef>0,feature]
    all <- nrow(to.plot)
    
    # to.plot <- to.plot[abs(coef)>0.0001]
    
    xlim <- max(c(0.02,abs(to.plot$coef)), na.rm=T)
    ylim <- max(c(5,to.plot$log_pval), na.rm=T)
    
    p <- ggplot(to.plot, aes(x=coef, y=log_pval)) +
      geom_point(aes(color=sig, size=sig)) +
      geom_segment(x=0, xend=0, y=0, yend=ylim, color="orange", linewidth=0.1) +
      scale_size_manual(values=c(0.50,1)) +
      # scale_x_continuous(limits=c(-xlim,xlim+1.5)) +
      # scale_y_continuous(limits=c(0,ylim+2)) +
      scale_color_manual(values=c("black","red")) +
      labs(x="limma coefficient", y=expression(paste("-log"[10],"(p.value)")), title=sprintf("%s; %s",i,j)) +
      annotate("text", x=0, y=ylim+1, size=4, label=sprintf("(%d)", all)) +
      annotate("text", x=-xlim, y=ylim+1, size=4, label=sprintf("%d (-)",length(negative_hits))) +
      annotate("text", x=xlim, y=ylim+1, size=4, label=sprintf("%d (+)",length(positive_hits))) +
      annotate("text", x=-xlim+(xlim/15), y=0, size=4, label="Down") +
      annotate("text", x=xlim-(xlim/15), y=0, size=4, label="Up") +
      # geom_text_repel(data=head(to.plot,n=50), aes(x=coef, y=log_pval, label=feature), size=2.5, segment.size=0.05, max.overlaps = Inf) +
      geom_text_repel(data=head(to.plot[coef>0],n=15), aes(x=coef, y=log_pval, label=feature), segment.size=0.05, size=2.50, max.overlaps = Inf) +
      geom_text_repel(data=head(to.plot[coef<0],n=15), aes(x=coef, y=log_pval, label=feature), segment.size=0.05, size=2.50, max.overlaps = Inf) +
      theme_classic() +
      theme(
        plot.title = element_text(hjust=0.5),
        axis.text = element_text(size=rel(1.0), color='black'),
        axis.title = element_text(size=rel(1.10), color='black'),
        legend.position="none"
      )
    
    pdf(file.path(outdir,sprintf("%s_%s_limma_rna_volcano_plot_chan2022.pdf",i,gsub("[ /#]","_",j))), width=8, height=7)
    print(p)
    dev.off()
  }
}


########################################
## Scatter plots comparing conditions ##
########################################

outdir <- file.path(io$basedir,"rna/results/limma/pdf/scatter_plots"); dir.create(outdir, showWarnings = F, recursive = T)

conditions.to.plot <- names(limma_list)

# i <- conditions.to.plot[1]; j <- conditions.to.plot[2]
for (i in conditions.to.plot) {
  for (j in conditions.to.plot) {
    if (i!=j) {
      
      foo <- copy(limma_list[[i]]) %>%
        .[is.na(coef),c("coef","pval_adj"):=list(0,1)] %>%
        .[,log_pval:=-log10(pval_adj+1e-100)]
      
      bar <- copy(limma_list[[j]]) %>%
        .[is.na(coef),c("coef","pval_adj"):=list(0,1)] %>%
        .[,log_pval:=-log10(pval_adj+1e-100)]
      
      to.plot <- merge(foo[,c("ens_id","feature","coef")],bar[,c("ens_id","feature","coef")], by=c("ens_id","feature"), suffixes = c("_A", "_B")) %>%
        .[,transcript_type:="Genes"] %>%
        .[feature%in%gene_sets[["Surface markers"]],transcript_type:="Surface markers"] %>%
        .[feature%in%gene_sets[["Epigenetic factors"]],transcript_type:="Epigenetic factors"] %>%
        .[feature%in%gene_sets[["Ribosomal"]],transcript_type:="Ribosomal"] %>%
        .[feature%in%gene_sets[["Mitochondrial"]],transcript_type:="Mitochondrial"] %>%
        .[feature%in%gene_sets[["DNA repair"]],transcript_type:="DNA repair"] %>%
        .[feature%in%gene_sets[["Cell cycle"]],transcript_type:="Cell cycle"] %>%
        .[feature%in%gene_sets[["ECM"]],transcript_type:="ECM"] %>%
        .[feature%in%gene_sets[["Transcription Factors"]],transcript_type:="Transcription Factors"] %>%
        .[,transcript_type:=factor(transcript_type, levels=c("Genes",names(gene_sets)))]
      
      to.plot <- to.plot %>%
        .[,tmp:=abs(coef_A)*abs(coef_B)] %>% setorder(-tmp) %>% .[,tmp:=NULL] %>% 
        .[abs(coef_A)>0.001 | abs(coef_B)>0.001]
        
      xlim.max <- max(abs(to.plot$coef_A))
      ylim.max <- max(abs(to.plot$coef_B))
      
      p <- ggplot(to.plot, aes(x=coef_A, y=coef_B)) +
        geom_point(size=0.15, color="gray30", data=to.plot[!(abs(coef_B)>=0.01 | abs(coef_A)>=0.01)][!feature%in%genes.to.highlight]) +
        geom_point(size=0.45, color="gray30", data=to.plot[abs(coef_B)>=0.01 | abs(coef_A)>=0.01][!feature%in%genes.to.highlight]) +
        geom_point(size=0.45, color="purple", data=to.plot[feature%in%genes.to.highlight]) +
        stat_smooth(method="lm") +
        stat_cor(aes(label = after_stat(r.label)),  method = "pearson") +
        geom_segment(x=-xlim.max, xend=xlim.max, y=0, yend=0, color="orange", linewidth=0.1) +
        geom_segment(x=0, xend=0, y=-ylim.max, yend=ylim.max, color="orange", linewidth=0.1) +
        labs(x=sprintf("limma coefficient (%s)",i), y=sprintf("limma coefficient (%s)",j), title="") +
        geom_text_repel(data=to.plot[feature%in%genes.to.highlight], aes(x=coef_A, y=coef_B, label=feature), color="purple", size=2.50, segment.size=0.1, max.overlaps = Inf) +
        facet_wrap(~transcript_type, nrow=2) +
        theme_classic() +
        theme(
          plot.title = element_text(hjust=0.5, size=rel(1.0)),
          axis.text = element_text(size=rel(0.75), color='black'),
          axis.title = element_text(size=rel(1.0), color='black'),
          legend.position="none"
        )
      
      for (k in unique(to.plot$transcript_type)) {
        p <- p +
          geom_text_repel(data=head(to.plot[transcript_type==k & coef_B>=0.01 & coef_A>=0.01 & !feature%in%genes.to.highlight],n=15), aes(x=coef_A, y=coef_B, label=feature), color="black", size=2.50, segment.size=0.1, max.overlaps = Inf) +
          geom_text_repel(data=head(to.plot[transcript_type==k & coef_B<(-0.01) & coef_A<=(-0.01) & !feature%in%genes.to.highlight],n=15), aes(x=coef_A, y=coef_B, label=feature), color="black", size=2.50, segment.size=0.1, max.overlaps = Inf)
      }
      
      pdf(file.path(outdir,sprintf("limma_expr_%s_vs_%s_scatter_plot_chan2022.pdf",i,j)), width=13, height=10)
      print(p)
      dev.off()
    }
  }
}
