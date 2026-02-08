source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# Options
opts$group_variable <- "Treatment"

# I/O
io$diff_results <- file.path(io$basedir,sprintf("results/rna/differential/%s/diff_Control_vs_TG.txt.gz",opts$group_variable))
io$outdir <- file.path(io$basedir,sprintf("results/rna/differential/%s/pdf",opts$group_variable)); dir.create(io$outdir, showWarnings = F)
io$outdir.volcano <- file.path(io$basedir,sprintf("results/rna/differential/%s/pdf/volcano_plots",opts$group_variable)); dir.create(io$outdir.volcano, showWarnings = F, recursive = T)
io$outdir.scatter <- file.path(io$basedir,sprintf("results/rna/differential/%s/pdf/scatter_plots",opts$group_variable)); dir.create(io$outdir.scatter, showWarnings = F, recursive = T)

########################
## Load gene metadata ##
########################

gene_metadata.dt <- fread(io$gene_metadata)

# io$TEs <- "/Users/rargelaguet/data/hg38_regulation/transposable_elements/TEtranscripts/TEtranscripts_TE_metadata.txt.gz"
TEs <- fread(io$TEs, header = F)[[1]]

TFs <- fread(io$TFs)[[1]]

surface_markers <- unique(fread(io$cell_surface_markers)[[1]])

##########################################
## Load differential expression results ##
##########################################

diff.dt <- fread(io$diff_results) %>%
  setnames("gene","feature") %>%
  # .[groupA%in%opts$groups & groupB%in%opts$groups] %>%
  .[is.na(logFC),c("logFC","padj_fdr"):=list(0,1)] %>%
  .[,transcript_type:="Genes"] %>%
  .[feature%in%TEs,transcript_type:="Transposable Elements"] %>%
  .[feature%in%TFs,transcript_type:="Transcription Factors"] %>%
  .[feature%in%surface_markers,transcript_type:="Surface markers"] %>%
  .[,transcript_type:=factor(transcript_type, levels=c("Genes","Transcription Factors","Surface markers","Transposable Elements"))]

table(diff.dt$transcript_type)

##########################
## Filter diff results ##
##########################

genes.to.remove <- c(grep("^MT-",gene_metadata.dt$symbol, value=T),grep("^RPS|^RPL",gene_metadata.dt$symbol, value=T))

diff.dt <- diff.dt[!feature%in%genes.to.remove]

###################
## Volcano plots ##
###################

i <- "Control"; j <- "TG"
# to.plot <- diff.dt[groupA==i & groupB==j] %>% 
to.plot <- diff.dt %>% 
  .[abs(logFC)>0.05 & padj_fdr<1] %>%
  .[,sig:=padj_fdr<0.01 & abs(logFC)>=1]

negative_hits <- to.plot[sig==TRUE & logFC<0,feature]
positive_hits <- to.plot[sig==TRUE & logFC>0,feature]
all <- nrow(to.plot)

xlim <- max(c(3,abs(to.plot$logFC)), na.rm=T)
ylim <- max(c(10,-log10(to.plot$padj_fdr+1e-100)), na.rm=T)

p <- ggplot(to.plot, aes(x=logFC, y=-log10(padj_fdr+1e-100))) +
  geom_point(aes(size=sig), color="gray40") +
  geom_segment(x=0, xend=0, y=0, yend=ylim, color="orange", linewidth=0.15) +
  # scale_color_manual(values=c("black","red")) +
  scale_size_manual(values=c(0.25,0.50)) +
  scale_x_continuous(limits=c(-xlim-1.5,xlim+1.5)) +
  scale_y_continuous(limits=c(0,ylim+2)) +
  labs(x="Log fold change", y=expression(paste("-log"[10],"(p.value)")), title=sprintf("%s vs %s",i,j)) +
  # annotate("text", x=0, y=ylim+1, size=4, label=sprintf("(%d)", all)) +
  # annotate("text", x=-xlim, y=ylim+1, size=4, label=sprintf("%d (-)",length(negative_hits))) +
  # annotate("text", x=xlim, y=ylim+1, size=4, label=sprintf("%d (+)",length(positive_hits))) +
  annotate("text", x=-xlim+(xlim/5), y=0, size=2, label=sprintf("Up in %s",i)) +
  annotate("text", x=xlim-(xlim/5), y=0, size=2, label=sprintf("Up in %s",j)) +
  ggrepel::geom_text_repel(data=head(to.plot[transcript_type=="Genes"],n=30), aes(x=logFC, y=-log10(padj_fdr+1e-100), label=feature), segment.size=0.05, size=2.50, max.overlaps = Inf) +
  ggrepel::geom_text_repel(data=head(to.plot[transcript_type=="Transposable Elements"],n=6), aes(x=logFC, y=-log10(padj_fdr+1e-100), label=feature), segment.size=0.05, size=2.50, max.overlaps = Inf) +
  ggrepel::geom_text_repel(data=head(to.plot[transcript_type=="Transcription Factors"],n=12), aes(x=logFC, y=-log10(padj_fdr+1e-100), label=feature), segment.size=0.05, size=2.50, max.overlaps = Inf) +
  ggrepel::geom_text_repel(data=head(to.plot[transcript_type=="Surface markers"],n=12), aes(x=logFC, y=-log10(padj_fdr+1e-100), label=feature), segment.size=0.05, size=2.50, max.overlaps = Inf) +
  facet_wrap(~transcript_type, nrow=1) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust=0.5, size=rel(1.0)),
    axis.text = element_text(size=rel(0.75), color='black'),
    axis.title = element_text(size=rel(1.0), color='black'),
    legend.position="none"
  )

pdf(file.path(io$outdir.volcano,sprintf("volcano_plot_%s_vs_%s_BJ_TG.pdf",i,j)), width=12, height=6)
print(p)
dev.off()
