# https://www.ArchRProject.com/bookdown/identifying-marker-peaks-with-archr.html

########################
## Load ArchR Project ##
########################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
} else {
  stop("Computer not recognised")
}

#####################
## Define settings ##
#####################

# I/O
# io$metadata <- paste0(io$basedir,"/processed/atac/archR/sample_metadata_after_archR.txt.gz")
io$outdir <- paste0(io$basedir,"/results/atac/archR/marker_peaks/all_peaks")

# Options
opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E8.5_rep1",
  "E8.5_rep2"
)

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE] %>%
  .[sample%in%opts$samples]
stopifnot(sample_metadata$archR_cell %in% rownames(ArchRProject))

# subset celltypes with sufficient number of cells
opts$min.cells <- 100

sample_metadata <- sample_metadata %>%
  .[,N:=.N,by=c("celltype.predicted")] %>% .[N>opts$min.cells] %>% .[,N:=NULL]
table(sample_metadata$celltype.predicted)

##################
## Subset ArchR ##
##################

ArchRProject.filt <- ArchRProject[sample_metadata$archR_cell,]

############################################
## Load results from marker gene analysis ##
############################################

io$marker.genes <- paste0(io$basedir,"/results/atac/archR/marker_peaks/all_peaks/markersPeaks.tsv.gz")
dt <- fread(io$marker.genes) %>% setorder(FDR)

# Define significant hits
dt[,sig:=abs(Log2FC)>=1 & FDR<0.1]

#####################################################
## Barplot of number of marker genes per cell type ##
#####################################################

to.plot <- dt %>%
  .[,.(Nmarkers=sum(sig==T & Log2FC>0)),by="celltype"] 

p <- ggbarplot(to.plot, x="celltype", y="Nmarkers", fill="celltype") +
  labs(x="", y="Number of marker peaks") +
  scale_fill_manual(values=opts$celltype.colors) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(colour="black",size=rel(0.8)),  
    axis.text.x = element_text(colour="black",size=rel(0.65), angle=45, hjust=1, vjust=1),  
  )

pdf(sprintf("%s/barplot_number_gene_markers_per_celltype.pdf",io$outdir), width = 9, height = 5)
print(p)
dev.off()

################################################################################
## Scatterplot of number of cells versus number of marker genes per cell type ##
################################################################################

celltype_numbers <- sample_metadata %>% 
  .[,.(Ncells=.N),by=c("celltype.predicted")] %>%
  setnames("celltype.predicted","celltype")

to.plot <- dt %>% 
  .[,.(Nmarkers=sum(sig==T & Log2FC>0)),by="celltype"] %>% 
  merge(celltype_numbers,by="celltype")

p <- ggscatter(to.plot, x="Nmarkers", y="Ncells", fill="celltype", shape=21, size=3, 
               add="reg.line", add.params = list(color="blue", fill="lightgray"), conf.int=TRUE) +
  labs(x="Number of marker genes", y="Number of cells") +
  scale_fill_manual(values=opts$celltype.colors) +
  theme(
    axis.text = element_text(colour="black",size=rel(0.65)),  
    legend.position = "none"
  )
p

###################
## Volcano plots ##
###################

for (i in unique(dt$celltype)) {
  
  to.plot <- dt[celltype==i] %>% .[,id:=sprintf("%s:%s-%s",chr,start,end)]
  negative_hits <- to.plot[sig==TRUE & Log2FC<0,id]
  positive_hits <- to.plot[sig==TRUE & Log2FC>0,id]
  all <- nrow(to.plot)
  
  xlim <- max(abs(to.plot$Log2FC), na.rm=T)
  ylim <- max(-log10(to.plot$FDR+1e-100), na.rm=T)
  
  p <- ggplot(to.plot, aes(x=Log2FC, y=-log10(FDR+1e-100))) +
    # geom_hline(yintercept = -log10(opts$threshold_fdr), color="blue") +
    geom_segment(aes(x=0, xend=0, y=0, yend=ylim-1), color="orange", size=0.5) +
    ggrastr::geom_point_rast(aes(color=sig, size=sig)) +
    # ggrepel::geom_text_repel(data=head(to.plot[sig==T],n=15), aes(x=Log2FC, y=-log10(FDR+1e-100), label=gene), size=3) +
    scale_color_manual(values=c("black","red")) +
    scale_size_manual(values=c(0.75,1.25)) +
    scale_x_continuous(limits=c(-xlim-0.5,xlim+0.5)) +
    scale_y_continuous(limits=c(0,ylim+6)) +
    annotate("text", x=0, y=ylim+6, size=4, label=sprintf("(%d)", all)) +
    annotate("text", x=-xlim-0.5, y=ylim+6, size=4, label=sprintf("%d (-)",length(negative_hits))) +
    annotate("text", x=xlim+0.5, y=ylim+6, size=4, label=sprintf("%d (+)",length(positive_hits))) +
    annotate("text", x=xlim-0.5, y=0, size=3, label=sprintf("Up in %s (N=%s)",i,celltype_numbers[celltype==i,Ncells])) +
    labs(x="Log fold change", y=expression(paste("-log"[10],"(p.value)"))) +
    theme_classic() +
    theme(
      axis.text = element_text(size=rel(0.75), color='black'),
      axis.title = element_text(size=rel(1.0), color='black'),
      legend.position="none"
    )
  
  # pdf(sprintf("%s/volcano_marker_genes_%s.pdf",io$outdir,i), width = 9, height = 6)
  png(sprintf("%s/volcano_marker_peaks_%s.png",io$outdir,i), width = 800, height = 500)
  print(p)
  dev.off()
}

####################
## Browser tracks ##
####################

gene <- dt[1,gene]

gr <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)
p <- plotBrowserTrack(ArchRProject.filt, 
                      groupBy = "celltype.predicted", 
                      geneSymbol = gene,
                      # features =  gr,
                      useMatrix = "GeneScoreMatrix",
                      upstream = 50000,
                      downstream = 50000
)

grid::grid.newpage()
grid::grid.draw(p$Tbx)

# plotPDF(p, name = "Plot-Tracks-With-Features", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)


#############
## Explore ##
#############

dt[,idx:=sprintf("%s:%s-%s",chr,start,end)]
dt[celltype=="Spinal_cord" & sig==T & MeanDiff>0.25] %>%  View
