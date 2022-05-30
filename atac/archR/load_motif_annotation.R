# opts$motif_annotation <- "Motif_cisbp" # "Motif_JASPAR2020" 

motif2gene_file <- sprintf("%s/Annotations/%s_TFs.txt.gz",io$archR.directory,opts$motif_annotation)

if (file.exists(motif2gene_file)) {

  motif2gene.dt <- fread(motif2gene_file)

} else {

  peakAnnotation <- readRDS(sprintf("%s/Annotations/peakAnnotation.rds",io$archR.directory))
  stopifnot(opts$motif_annotation%in%names(peakAnnotation))
  motif2gene.dt <- peakAnnotation[[opts$motif_annotation]]$motifSummary %>%
    as.data.table(keep.rownames = T) %>% setnames("rn","motif") %>% .[,strand:=NULL] %>% setnames("symbol","gene")

  # Rename genes
  if (grepl("cisbp",opts$motif_annotation, ignore.case = T)) {
    
    tf2gene_rename <- c(
      "Tcfe"="Tfe", "Nkx1"="Nkx1-", "Nkx2"="Nkx2-", "Nkx3"="Nkx3-", "Nkx4"="Nkx4-", "Nkx5"="Nkx5-", "Nkx6"="Nkx6-", "Foxf1a"="Foxf1",
      "Hmga1rs1"="rs1", "Mycl1$"="Mycl", "Dux$"="Duxf3", "Duxbl$"="Duxbl1", "Pit1$"="Prop1",
      "ENSMUSG00000079994"="Sox1", "Tcfap"="Tfap"
    )
    
    motif2gene.dt[,gene:=stringr::str_replace_all(gene,tf2gene_rename)]
    
  } else if (grepl("JASPAR",opts$motif_annotation, ignore.case = T)) {
    
    # conflictive motifs: fusion proteins (UN::JUNB) and versions (TFAP2A(var.2))
    # stop("To-do")
    tf2gene_rename <- c("TBXT"="T")
    motif2gene.dt[,gene:=stringr::str_replace_all(gene,tf2gene_rename)]
    
    # for JASPAR motifs
    motif2gene.dt[,motif:=str_replace(motif,"\\.VAR\\.","\\.var\\."),]
    
  } else {
    stop("Motif annotation not recognised")
  }

  motif2gene.dt[,c("motif","gene"):=list(toupper(motif),toupper(gene))]

  # Save
  fwrite(motif2gene.dt, motif2gene_file, sep="\t", quote=F)
}


