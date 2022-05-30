#####################
## Define settings ##
#####################

# Load default settings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/utils.R")
  source("/Users/ricard/gastrulation_multiome_10x/atac/archR/differential/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/utils.R")
  source("/homes/ricard/gastrulation_multiome_10x/atac/archR/differential/utils.R")
} else {
  stop("Computer not recognised")
}

# I/O
io$outdir <- paste0(io$basedir,"/results/atac/archR/differential/PeakMatrix/pdf"); dir.create(io$outdir, showWarnings = F)

# Options
# opts$groups <- strsplit(list.files(io$diff.dir, pattern="*.gz"),"_vs_") %>% map(~ .[[1]]) %>% unlist %>% unique
opts$celltypes <- c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  # "PGC",
  # "Anterior_Primitive_Streak",
  "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  # "Mixed_mesoderm",
  "Intermediate_mesoderm",
  "Caudal_Mesoderm",
  "Paraxial_mesoderm",
  "Somitic_mesoderm",
  "Pharyngeal_mesoderm",
  "Cardiomyocytes",
  "Allantois",
  "ExE_mesoderm",
  "Mesenchyme",
  "Haematoendothelial_progenitors",
  "Endothelium",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  # "Erythroid1",
  # "Erythroid2",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
  "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm"
  # "Visceral_endoderm",
  # "ExE_endoderm",
  # "ExE_ectoderm",
  # "Parietal_endoderm"
)

opts$celltypes <- c("Rostral_neurectoderm", "Forebrain_Midbrain_Hindbrain", "Neural_crest")

opts$min.MeanDiff <- 0.15
opts$fdr <- 0.01
opts$atac.matrix <- "PeakMatrix"

##################
## Load results ##
##################

# NOTE: POSITIVE HITS ARE MORE UPREGULATED IN CELLTYPE A
# NOTE: NEGATIVE HITS ARE MORE UPREGULATED IN CELLTYPE B

dt <- opts$celltypes %>% map(function(i) { opts$celltypes %>% map(function(j) {
  file <- sprintf("%s/%s_%s_vs_%s.txt.gz", io$archR.peak.differential.dir,opts$atac.matrix,i,j)
  if (file.exists(file)) {
    fread(file, select = c(1,2,3)) %>% 
      .[,.(nhits_upregulated=sum(MeanDiff>=opts$min.MeanDiff & FDR<=opts$fdr,na.rm=T), 
           nhits_downregulated=sum(MeanDiff<=-(opts$min.MeanDiff) & FDR<=opts$fdr,na.rm=T))] %>%
      .[,nhits:=nhits_upregulated+nhits_downregulated] %>%
      .[,c("celltypeA","celltypeB"):=list(i,j)] %>%
      return
  }
}) %>% rbindlist }) %>% rbindlist


dt <- dt %>% rbind(
  dt %>% copy %>% setnames(c("nhits_downregulated","nhits_upregulated","nhits","celltypeB","celltypeA")) %>% .[,c("nhits_upregulated","nhits_downregulated","nhits","celltypeA","celltypeB")]
)

##################
## Circle plots ##
##################

to.plot <- dt %>% copy %>% 
  .[,fraction_upregulated_hits:=nhits_upregulated/nhits] %>%
  .[nhits>=8e3,nhits:=8e3] %>%
  .[fraction_upregulated_hits>=0.90,fraction_upregulated_hits:=0.90] %>%
  .[fraction_upregulated_hits<=0.10,fraction_upregulated_hits:=0.10] %>%
  .[,nhits_scaled:=minmax.normalisation(nhits)]

# TFs.to.plot <- cor_dt[,.N,by="TF"] %>% setorder(-N) %>% head(n=100) %>% .$TF
# to.plot <- to.plot[TF%in%TFs.to.plot] %>% .[,TF:=factor(TF,levels=TFs.to.plot)]

to.plot[,celltypeA:=factor(celltypeA,levels=rev(opts$celltypes))]
to.plot[,celltypeB:=factor(celltypeB,levels=opts$celltypes)]

p <- ggplot(to.plot, aes_string(x="celltypeB", y="celltypeA", color="fraction_upregulated_hits", size="nhits_scaled")) +
  geom_point()  +    ## geom_point for circle illusion
  # scale_colour_gradientn(colours = terrain.colors(10)) +
  scale_color_gradient2(low = "blue", mid="gray90", high = "red", midpoint=0.5) +
  guides(x = guide_axis(angle = 90)) +
  labs(x="", y="") +
  theme_classic() +
  guides(size=F) +
  theme(
    axis.text = element_text(colour="black",size=rel(0.75)),
    axis.ticks.x = element_blank(),
    legend.position = "top"
  )

pdf(sprintf("%s/atac_fraction_upregulated_peaks_circles.pdf",io$outdir), width=9, height=8)
print(p)
dev.off()

###############
## Box plots ##
###############

to.plot2 <- to.plot %>% .[celltypeA!=celltypeB & nhits>=100]
celltype.order <- to.plot2[,mean(fraction_upregulated_hits),by="celltypeA"] %>% setorder(-V1) %>% .$celltypeA
to.plot2[,celltypeA:=factor(celltypeA,levels=celltype.order)]

p <- ggplot(to.plot2, aes(x=celltypeA, y=fraction_upregulated_hits)) +
  geom_point(aes(fill = celltypeA), shape=21, size=1, alpha=0.5) +
  geom_boxplot(aes(fill = celltypeA), alpha=0.5, outlier.shape=NA, coef=1) +
  coord_flip(ylim=c(0,1)) +
  geom_hline(yintercept=0.50, linetype="dashed", size=0.5) +
  scale_fill_manual(values=opts$celltype.colors, drop=F) +
  theme_classic() +
  labs(y="Fraction of upregulated peaks", x="") +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  )


pdf(sprintf("%s/atac_fraction_upregulated_peaks_boxplots.pdf",io$outdir), width=5, height=7)
print(p)
dev.off()


###################
## Volcano plots ##
###################

dt <- opts$celltypes %>% map(function(i) { opts$celltypes %>% map(function(j) {
  file <- sprintf("%s/PeakMatrix_%s_vs_%s.txt.gz", io$archR.peak.differential.dir,i,j)
  if (file.exists(file)) {
    fread(file) %>% .[,c("celltypeA","celltypeB"):=list(as.factor(i),as.factor(j))] %>%
      return
  }
}) %>% rbindlist }) %>% rbindlist %>%
  .[,sig:=FALSE] %>% .[abs(MeanDiff)>opts$min.MeanDiff & FDR<opts$fdr,sig:=TRUE] %>%
  .[,direction:=c("up","down")[as.numeric(MeanDiff<0)+1]]  # up = higher accessibility in celltype A

for (i in unique(dt$celltypeA)) {
  for (j in unique(dt[celltypeA==i,celltypeB])) {
    to.plot <- dt[celltypeA==i & celltypeB==j] %>% .[!is.na(sig)] %>% .[,MeanDiff:=100*MeanDiff]

    p <- gg_volcano_plot(to.plot, top_genes=0, label_groups = c(i,j))

    # pdf(sprintf("%s/%s_%s_%s_volcano.pdf",io$outdir,i,opts$groupA,opts$groupB), width=9, height=5)
    png(sprintf("%s/%s_%s_DA_peaks_volcano.png",io$outdir,i,j), width=900, height=500)
    print(p)
    dev.off()
  }
}

