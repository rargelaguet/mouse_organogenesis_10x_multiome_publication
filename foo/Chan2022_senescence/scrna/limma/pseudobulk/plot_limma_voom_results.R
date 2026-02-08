source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$limma_results <- file.path(io$basedir,"results/limma_voom/limma_results_passage.txt.gz")
# io$limma_results <- file.path(io$basedir,"results/limma_voom/limma_results_passage_proliferation_corrected.txt.gz")
io$outdir <- file.path(io$basedir,"results/limma_voom/pdf"); dir.create(io$outdir, showWarnings = F, recursive = T)

# Options
opts$coef <- 0.01
opts$pvalue <- 0.001

########################
## Load limma results ##
########################

limma.list <- list("passage"=fread(io$limma_results))
# limma.dt <- fread(io$limma_results)

########################
## Load gene metadata ##
########################

source("load_gene_sets.R")

###################
## Volcano plots ##
###################

contrasts <- names(limma.list)

# i <- contrasts[1]
for (i in contrasts) {
  
  to.plot <- copy(limma.list[[i]]) %>%
    .[is.na(coef),c("coef","pval_adj"):=list(0,1)] %>%
    setorder(pval_adj) %>%
    .[,sig:=pval_adj<opts$pvalue & abs(coef)>=opts$coef] %>%
    .[,log_pval:=-log10(pval_adj+1e-100)] %>%
    .[,transcript_type:="Other genes"]
  
  for (j in names(gene_sets)) {
    to.plot <- to.plot %>%
    .[feature%in%gene_sets[[j]],transcript_type:=j]
  }
  to.plot <- to.plot %>% .[,transcript_type:=factor(transcript_type, levels=c("Other genes",names(gene_sets)))]
  
  negative_hits <- to.plot[sig==TRUE & coef<0,feature]
  positive_hits <- to.plot[sig==TRUE & coef>0,feature]
  all <- nrow(to.plot)
  
  to.plot <- to.plot[abs(coef)>0.001]
  
  xlim <- max(c(0.05,abs(to.plot$coef)), na.rm=T)
  ylim <- max(c(5,to.plot$log_pval), na.rm=T)
  
  tmp <- to.plot[,.(N=as.character(.N), Npos=as.character(sum(sig==TRUE & coef>0)), Nneg=as.character(sum(sig==TRUE & coef<0))),by="transcript_type"]
  tmp %>% .[,Npos:=sprintf("%s (+)",Npos)] %>% .[,Nneg:=sprintf("%s (-)",Nneg)]
  
  p <- ggplot(to.plot, aes(x=coef, y=log_pval)) +
    geom_point(aes(color=sig, size=sig)) +
    geom_segment(x=0, xend=0, y=0, yend=ylim, color="orange", linewidth=0.1) +
    scale_color_manual(values=c("black","red")) +
    scale_size_manual(values=c(0.20,0.40)) +
    # scale_x_continuous(limits=c(-xlim,xlim+1.5)) +
    scale_y_continuous(limits=c(0,ylim+1.5)) +
    labs(x="limma coefficient", y=expression(paste("-log"[10],"(p.value)")), title=i) +
    geom_text(aes(label=N), x=0, y=ylim+1, size=3, data=tmp) +
    geom_text(aes(label=Nneg), x=-xlim+(xlim/3), y=ylim+1, size=3, data=tmp) +
    geom_text(aes(label=Npos), x=xlim-(xlim/3), y=ylim+1, size=3, data=tmp) +
    annotate("text", x=-xlim+xlim/1.75, y=0, size=3, label="Down passage") +
    annotate("text", x=xlim-xlim/2.5, y=0, size=3, label="Up passage") +
    facet_wrap(~transcript_type, nrow=3) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust=0.5),
      strip.text = element_text(size=rel(1.05), color='black'),
      axis.text = element_text(size=rel(0.75), color='black'),
      axis.title = element_text(size=rel(1.0), color='black'),
      legend.position="none"
    )
  
  for (j in unique(to.plot$transcript_type)) {
    p <- p + geom_text_repel(data=head(to.plot[transcript_type==j & sig==TRUE],n=25), aes(x=coef, y=log_pval, label=feature), segment.size=0.05, size=2.50, max.overlaps = Inf)
  }
  
  pdf(file.path(io$outdir,sprintf("%s_limma_rna_volcano_plot_Sturm2022.pdf",i)), width=14, height=10)
  print(p)
  dev.off()
}

################################
## Volcano plots per gene set ##
################################

outdir <- file.path(io$outdir,"per_gene_set"); dir.create(outdir, showWarnings = F, recursive = T)

contrasts <- names(limma.list)[2]

# i <- contrasts[1]
for (i in contrasts) {
  for (j in names(gene_sets)) {
    
    to.plot <- limma.list[[i]] %>%
      .[feature%in%gene_sets[[j]]] %>%
      .[is.na(coef),c("coef","pval_adj"):=list(0,1)] %>%
      setorder(pval_adj) %>%
      .[,sig:=pval_adj<opts$pval & abs(coef)>=opts$coef] %>%
      .[,log_pval:=-log10(pval_adj+1e-300)]
    
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
    
    pdf(file.path(outdir,sprintf("%s_%s_limma_rna_volcano_plot_Sturm2022.pdf",i,j)), width=7, height=6)
    print(p)
    dev.off()
  }
}

######################################
## Volcano plots per gene signature ##
######################################

outdir <- file.path(io$outdir,"gene_signatures"); dir.create(outdir, showWarnings = F, recursive = T)

source(here::here("gene_signatures/load_manual_gene_signatures.R"))
source(here::here("gene_signatures/load_gene_signatures.R"))

# gene_sets_to_plot <- names(gene_signatures_list)
gene_sets_to_plot <- grep("HALLMARK_INTERFERON",names(gene_signatures_list),ignore.case = T, value=T)

contrast <- "passage"

# i <- "hr3#ISR"
for (i in gene_sets_to_plot) {
  
  to.plot <- limma.list[[contrast]] %>%
    .[feature%in%gene_signatures_list[[i]]] %>%
    .[is.na(coef),c("coef","pval_adj"):=list(0,1)] %>%
    .[,sig:=pval_adj<opts$pvalue & abs(coef)>=opts$coef] %>%
    .[,log_pval:=-log10(pval_adj+1e-100)]
  
  
  negative_hits <- to.plot[sig==TRUE & coef<0,feature]
  positive_hits <- to.plot[sig==TRUE & coef>0,feature]
  all <- nrow(to.plot)
  
  xlim <- max(c(0.01,abs(to.plot$coef)), na.rm=T)
  ylim <- max(c(5,to.plot$log_pval), na.rm=T)
  
  p <- ggplot(to.plot, aes(x=coef, y=log_pval)) +
    geom_point(aes(color=sig, size=sig)) +
    geom_segment(x=0, xend=0, y=0, yend=ylim, color="orange", linewidth=0.1) +
    scale_size_manual(values=c(0.50,1)) +
    # scale_x_continuous(limits=c(-xlim,xlim+1.5)) +
    # scale_y_continuous(limits=c(0,ylim+2)) +
    scale_color_manual(values=c("black","red")) +
    labs(x="Coefficient", y=expression(paste("-log"[10],"(p.value)")), title=i) +
    annotate("text", x=0, y=ylim+1, size=4, label=sprintf("(%d)", all)) +
    annotate("text", x=-xlim, y=ylim+1, size=4, label=sprintf("%d (-)",length(negative_hits))) +
    annotate("text", x=xlim, y=ylim+1, size=4, label=sprintf("%d (+)",length(positive_hits))) +
    annotate("text", x=-xlim+(xlim/15), y=0, size=3, label="Passage down") +
    annotate("text", x=xlim-(xlim/15), y=0, size=3, label="Passage up") +
    geom_text_repel(data=head(to.plot,n=20), aes(label=feature), segment.size=0.05, size=2.50, max.overlaps = Inf) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust=0.5),
      axis.text = element_text(size=rel(0.75), color='black'),
      axis.title = element_text(size=rel(1.0), color='black'),
      legend.position="none"
    )
  
  pdf(file.path(outdir,sprintf("%s_limma_rna_volcano_plot_Sturm2022.pdf",i)), width=8, height=7)
  print(p)
  dev.off()
}
