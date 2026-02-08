source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

io$var_results <- file.path(io$basedir,"results/rna/variance_partition/variance_partition.rds")
io$outdir <- paste0(io$basedir,"/results/rna/variance_partition"); dir.create(io$outdir, showWarnings = F, recursive = T)

#####################################
## Load variance partition results ##
#####################################

varPart <- readRDS(io$var_results)

###########################################################################
## Boxplots of variance explained by each covariate, including residuals ##
###########################################################################

to.plot <- as.matrix(varPart) %>% as.data.table(keep.rownames = T) %>%
  setnames("rn","gene") %>%
  melt(id.vars="gene") %>%
  .[,value:=100*value]

p <- ggplot(to.plot, aes(x=variable, y=value, fill=variable)) +
  geom_violin(width=1, alpha=0.8) +
  geom_boxplot(width=0.25, outlier.shape=NA, alpha=0.9) +
  theme_classic() +
  labs(x="" , y="Variance explained (%)") +
  theme(
    axis.text = element_text(colour="black",size=rel(1.0)),
    axis.title.y = element_text(colour="black",size=rel(1.2))
  )

print(p)


#############################################
## Barplots of variance explained per gene ##
#############################################

genes.to.plot <- rownames(varPart) %>% head(n=25)
to.plot.bar <- to.plot[gene%in%genes.to.plot]

# to.plot.bar <- to.plot[gene%in%head(unique(to.plot$gene),n=25)]# %>%
# .[,gene:=factor(gene,levels=gene[order(gene)])

p <- ggplot(to.plot.bar, aes(x=gene, y=value, fill=variable)) +
  geom_bar(stat="identity", color="black") +
  coord_flip() +
  theme_classic() +
  labs(x="" , y="Variance explained (%)") +
  theme(
    axis.text.x = element_text(colour="black",size=rel(1.0)),
    axis.ticks.y = element_line(size=rel(0.5)),
    axis.text.y = element_text(colour="black",size=rel(0.75)),
    legend.position="right",
    legend.title = element_blank(),
    legend.text = element_text(size=rel(1.0))
  )

# pdf(sprintf("%s/%s.pdf",io$outdir,i), width=5, height=3.5, useDingbats = F)
print(p)
# dev.off()