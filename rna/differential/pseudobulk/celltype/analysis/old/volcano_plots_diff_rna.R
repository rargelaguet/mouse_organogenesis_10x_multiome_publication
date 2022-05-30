#####################
## Define settings ##
#####################

# Load default settings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/utils.R")
  source("/Users/ricard/gastrulation_multiome_10x/rna/differential/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/utils.R")
  source("/homes/ricard/gastrulation_multiome_10x/rna/differential/utils.R")
} else {
  stop("Computer not recognised")
}

# I/O
io$outdir <- paste0(io$basedir,"/results/rna/differential/pdf"); dir.create(io$outdir, showWarnings = F)

# Options
# opts$groups <- strsplit(list.files(io$diff.dir, pattern="*.gz"),"_vs_") %>% map(~ .[[1]]) %>% unlist %>% unique
opts$celltypes <- c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  "PGC",
  "Anterior_Primitive_Streak",
  "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  "Mixed_mesoderm",
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
  "Erythroid1",
  "Erythroid2",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
  "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm",
  "Visceral_endoderm",
  "ExE_endoderm",
  "ExE_ectoderm",
  "Parietal_endoderm"
)# %>% head(n=4)

# opts$celltypes <- c("Rostral_neurectoderm", "Forebrain_Midbrain_Hindbrain", "Neural_crest")

opts$min.logFC <- 1
opts$fdr <- 0.01

##################
## Load results ##
##################

# NOTE: POSITIVE logFC ARE MORE EXPRESSED IN CELLTYPE B
# NOTE: NEGATIVE logFC ARE MORE EXPRESSED IN CELLTYPE A

io$rna.file <- paste0(io$rna.differential,"/rna_diff_precomputed.txt.gz")

if (file.exists(io$rna.file)) {
  rna.dt <- fread(io$rna.file)
} else {
  rna.dt <- opts$celltypes %>% map(function(i) { opts$celltypes %>% map(function(j) {
    file <- sprintf("%s/%s_vs_%s.txt.gz", io$rna.differential,i,j)
    if (file.exists(file)) {
      fread(file, select = c(1,2)) %>% 
        .[,.(nhits_positive=sum(logFC>=opts$rna.min_logFC,na.rm=T), 
             nhits_negative=sum(logFC<=-(opts$rna.min_logFC),na.rm=T))] %>%
        .[,nhits:=nhits_positive+nhits_negative] %>%
        .[,c("celltypeA","celltypeB"):=list(i,j)] %>%
        return
    }
  }) %>% rbindlist }) %>% rbindlist %>%
    .[,celltypeA:=factor(celltypeA,levels=opts$celltypes)] %>%
    .[,celltypeB:=factor(celltypeB,levels=opts$celltypes)]
  
  rna.dt <- rna.dt %>% rbind(
    rna.dt %>% copy %>% setnames(c("nhits_negative","nhits_positive","nhits","celltypeB","celltypeA")) %>% .[,c("nhits_positive","nhits_negative","nhits","celltypeA","celltypeB")]
  )
  fwrite(rna.dt, io$rna.file)
}


##################
## Circle plots ##
##################

to.plot <- rna.dt %>% copy %>% 
  .[,fraction_positive_hits:=nhits_positive/nhits] %>%
  # .[nhits>=8e3,nhits:=8e3] %>%
  # .[fraction_positive_hits>=0.90,fraction_positive_hits:=0.90] %>%
  # .[fraction_positive_hits<=0.10,fraction_positive_hits:=0.10] %>%
  .[,nhits_scaled:=minmax.normalisation(nhits)]

# TFs.to.plot <- cor_dt[,.N,by="TF"] %>% setorder(-N) %>% head(n=100) %>% .$TF
# to.plot <- to.plot[TF%in%TFs.to.plot] %>% .[,TF:=factor(TF,levels=TFs.to.plot)]

to.plot[,celltypeA:=factor(celltypeA,levels=rev(opts$celltypes))]
to.plot[,celltypeB:=factor(celltypeB,levels=opts$celltypes)]

p <- ggplot(to.plot, aes_string(x="celltypeB", y="celltypeA", color="fraction_positive_hits", size="nhits_scaled")) +
  geom_point() +    ## geom_point for circle illusion
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

pdf(sprintf("%s/rna_fraction_positive_genes_circles.pdf",io$outdir), width=9, height=8)
print(p)
dev.off()

###############
## Box plots ##
###############

to.plot2 <- to.plot %>% .[celltypeA!=celltypeB & nhits>=100]
celltype.order <- to.plot2[,mean(fraction_positive_hits),by="celltypeA"] %>% setorder(-V1) %>% .$celltypeA
to.plot2[,celltypeA:=factor(celltypeA,levels=celltype.order)]

p <- ggplot(to.plot2, aes(x=celltypeA, y=fraction_positive_hits)) +
  geom_point(aes(fill = celltypeA), shape=21, size=1, alpha=0.5) +
  geom_boxplot(aes(fill = celltypeA), alpha=0.5, outlier.shape=NA, coef=1) +
  coord_flip(ylim=c(0,1)) +
  geom_hline(yintercept=0.50, linetype="dashed", size=0.5) +
  scale_fill_manual(values=opts$celltype.colors, drop=F) +
  theme_classic() +
  labs(y="Fraction of positive genes", x="") +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  )


pdf(sprintf("%s/rna_fraction_positive_genes_boxplots.pdf",io$outdir), width=5, height=7)
print(p)
dev.off()


###################
## Volcano plots ##
###################

opts$celltypes %>% walk(function(i) { opts$celltypes %>% walk(function(j) {
  print(sprintf("%s vs %s",i,j))
  file <- sprintf("%s/%s_vs_%s.txt.gz", io$rna.differential,i,j)
  if (file.exists(file)) {
    to.plot <- fread(file, select=c(1,2,4)) %>% 
      .[,c("celltypeA","celltypeB"):=list(as.factor(i),as.factor(j))] %>%
      .[,sig:=FALSE] %>% .[abs(logFC)>opts$min.logFC & padj_fdr<opts$fdr,sig:=TRUE] %>% 
      .[!is.na(sig)]
      
      p <- gg_volcano_plot(to.plot, top_genes=15, label_groups = c(i,j))
      
      # pdf(sprintf("%s/%s_%s_%s_volcano.pdf",io$outdir,i,opts$groupA,opts$groupB), width=9, height=5)
      png(sprintf("%s/%s_%s_DE_genes_volcano.png",io$outdir,i,j), width=700, height=400)
      print(p)
      dev.off()
    }
  })
})

