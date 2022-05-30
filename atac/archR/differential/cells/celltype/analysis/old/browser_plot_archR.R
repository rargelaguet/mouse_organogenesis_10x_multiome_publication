####################
## Browser tracks ##
####################

# i <- "chr2:39483639-39484239"

opts$extend.upstream <- 2500
opts$extend.downstream <- 2500
opts$tileSize <- 50

# Ugly hack
# rename <- paste(1:length(opts$celltypes),opts$celltypes,sep="_")
# names(rename) <- opts$celltypes
# ArchRProject.filt$celltype.predicted <- stringr::str_replace_all(ArchRProject.filt$celltype.predicted,rename)

for (i in unique(markers_peaks.dt$idx) %>% head(n=5)) {
  
  # Fetch GRanges
  to.plot <- peakset.gr[peakset.gr$idx==i]
  start(to.plot) <- start(to.plot) - opts$extend.upstream
  end(to.plot) <- end(to.plot) + opts$extend.downstream
  
  # Plot
  p <- plotBrowserTrack(
    ArchRProj = ArchRProject.filt, 
    region = to.plot,
    groupBy = "celltype.predicted", 
    tileSize = opts$tileSize,
    pal = opts$celltype.colors,
    plotSummary = c("bulkTrack", "featureTrack", "geneTrack"),
    sizes = c(10, 1.5, 1.5),
  )
  
  # grid::grid.newpage()
  
  pdf(sprintf("%s/%s_BrowserTrack.pdf",io$outdir,gsub(":","-",i)), width = 9, height = 5)
  grid::grid.draw(p)
  dev.off()
}