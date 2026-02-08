# here::i_am("differential/differential_gene_expr.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$TEs <- "/Users/rargelaguet/data/hg38_regulation/transposable_elements/scTE/scTE_TEs.txt"
io$diff_results <- file.path(io$basedir,"results/rna/differential/metacells/alias/parsed/diff_expr_results.txt.gz")
io$outdir <- file.path(io$basedir,"results/rna/differential/metacells/alias/pdf"); dir.create(io$outdir, showWarnings = F)

# Options
opts$groups <- opts$aliases

opts$group_variable <- "alias"

###############
## Load data ##
###############

diff.dt <- fread(io$diff_results) %>%
  setnames("gene","feature") %>%
  .[,sig:=padj_fdr<0.01 & abs(logFC)>=2]

TEs <- fread(io$TEs, header = F)[[1]]
diff.dt %>% .[,transcript_type:=c("gene","te")[as.numeric(feature%in%TEs)+1]]

###################
## Volcano plots ##
###################

# i <- opts$groups[1]; j <- opts$groups[2]
for (i in opts$groups) {
  for (j in opts$groups) {
    if (i!=j) {
      
      groupA <- i; groupB <- j
      
      # i <- "ERV1"
      to.plot <- diff.dt[groupA==i & groupB==j] %>%
        .[,transcript_type_plot:=ifelse(transcript_type=="gene","Genes","Transposable Elements")]
      
      negative_hits <- to.plot[sig==TRUE & logFC<0,feature]
      positive_hits <- to.plot[sig==TRUE & logFC>0,feature]
      all <- nrow(to.plot)
      
      xlim <- max(c(3,abs(to.plot$logFC)), na.rm=T)
      ylim <- max(c(10,-log10(to.plot$padj_fdr+1e-100)), na.rm=T)
      
      p <- ggplot(to.plot, aes(x=logFC, y=-log10(padj_fdr+1e-100))) +
        geom_point(aes(color=sig, size=sig)) +
        geom_segment(x=0, xend=0, y=0, yend=ylim, color="orange", linewidth=0.15) +
        scale_color_manual(values=c("black","red")) +
        scale_size_manual(values=c(0.30,0.90)) +
        scale_x_continuous(limits=c(-xlim-1.5,xlim+1.5)) +
        scale_y_continuous(limits=c(0,ylim+2)) +
        labs(x="Log fold change", y=expression(paste("-log"[10],"(p.value)")), title=sprintf("%s vs %s",i,j)) +
        # annotate("text", x=0, y=ylim+1, size=4, label=sprintf("(%d)", all)) +
        # annotate("text", x=-xlim, y=ylim+1, size=4, label=sprintf("%d (-)",length(negative_hits))) +
        # annotate("text", x=xlim, y=ylim+1, size=4, label=sprintf("%d (+)",length(positive_hits))) +
        annotate("text", x=-xlim+(xlim/5), y=0, size=2, label=sprintf("Up in %s",groupA)) +
        annotate("text", x=xlim-(xlim/5), y=0, size=2, label=sprintf("Up in %s",groupB)) +
        ggrepel::geom_text_repel(data=head(to.plot[sig==T & transcript_type=="gene"],n=50), aes(x=logFC, y=-log10(padj_fdr+1e-100), label=feature), segment.size=0.05, size=2.25, max.overlaps = Inf) +
        ggrepel::geom_text_repel(data=head(to.plot[transcript_type=="te"],n=10), aes(x=logFC, y=-log10(padj_fdr+1e-100), label=feature), segment.size=0.05, size=2.25, max.overlaps = Inf) +
        facet_wrap(~transcript_type_plot, nrow=1) +
        theme_classic() +
        theme(
          plot.title = element_text(hjust=0.5, size=rel(1.0)),
          axis.text = element_text(size=rel(0.75), color='black'),
          axis.title = element_text(size=rel(1.0), color='black'),
          legend.position="none"
        )
      
      pdf(file.path(io$outdir,sprintf("volcano_plot_%s_vs_%s_Xing2020.pdf",i,j)), width=10, height=6)
      print(p)
      dev.off()
    }
  }
}



#######################################
## Scatterplots pairwise comparisons ##
#######################################


###########################
## Compare to BJ RNA-seq ##
###########################

diff_bj_24h.dt <- fread("/Users/rargelaguet/data/BJ_fibroblasts_reprogramming/results/differential/genes/48h_diff.txt.gz")
# diff_bj_2D.dt <- fread("/Users/rargelaguet/data/BJ_fibroblasts_reprogramming/results/differential/genes/48h_diff.txt.gz")

#########################################
## Compare to BJ reprogramming (CHAMP) ##
#########################################

diff_bj_champ_d2_vs_d0.dt <- fread("/Users/rargelaguet/data/BJ_fibroblasts_reprogramming_scrna_altos/results/differential/metacells/group/diff_group_D02_Dox_vs_D00_Dox.txt.gz") %>%
  setnames("gene","feature")

features <- intersect(diff_bj_champ_d2_vs_d0.dt$feature, unique(diff.dt$feature))
to.plot <- rbind(
  diff.dt[feature%in%features & groupA=="D2" & groupB=="D0",c("feature","transcript_type","logFC")] %>% .[,dataset:="Xing2020"],
  diff_bj_champ_d2_vs_d0.dt[feature%in%features,c("feature","transcript_type","logFC")] %>% .[,dataset:="CHAMP"]
) %>%  .[is.na(logFC),logFC:=0] %>%
  dcast(feature+transcript_type~dataset, value.var="logFC")

xlim.max <- max(abs(to.plot$CHAMP), na.rm=T)
ylim.max <- max(abs(to.plot$Xing2020), na.rm=T)

p <- ggplot(to.plot, aes(x=CHAMP, y=Xing2020)) +
  # geom_point(size=0.30, width=0.01, height=0.01) +
  geom_point(size=1, color="black", data=to.plot[abs(Xing2020)>=1.5 & abs(CHAMP)>=1.5]) +
  ggrastr::geom_point_rast(size=0.15, color="gray70", data=to.plot[!(abs(Xing2020)>=1.5 & abs(CHAMP)>=1.5)]) +
  # stat_smooth(method="lm", alpha=0.5, color="gray30") +
  stat_cor(method = "pearson") +
  geom_segment(x=-xlim.max, xend=xlim.max, y=0, yend=0, color="orange", size=0.1) +
  geom_segment(x=0, xend=0, y=-ylim.max, yend=ylim.max, color="orange", size=0.1) +
  ggrepel::geom_text_repel(data=to.plot[abs(Xing2020)>=1.5 & abs(CHAMP)>=1.5], aes(x=CHAMP, y=Xing2020, label=feature), color="black", size=2, segment.size=0.1, max.overlaps = Inf) +
  labs(x="logFC (D2 vs D0, CHAMP)", y="logFC (D2 vs D0, Xing2020)", title="") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(size=rel(1.0), color="black"),
    axis.text.y = element_text(size=rel(1.0), color="black"),
    legend.position = "none",
    legend.title = element_blank()
  )

pdf(file.path(io$outdir,"scatterplot_Xing2020_vs_CHAMP_D0_vs_D2.pdf"), width=8, height=6)
print(p)
dev.off()




