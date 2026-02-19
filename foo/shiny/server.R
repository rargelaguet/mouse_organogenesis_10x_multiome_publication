###############
## Shiny app ##
###############

server <- function(input, output, session) {
 
  #####################
  ## Gene expression ##
  #####################

  updateSelectizeInput(session = session, inputId = 'gene_expr_gene', choices = genes, server = TRUE, selected = "CD70")
  
  plot_gene_expr <- reactive({
    
    ## START TEST ##
    # input <- list()
    # input$gene_expr_gene <- "CD70"
    # input$gene_expr_timepoints <- timepoints
    # input$gene_expr_donors <- donors
    # input$gene_expr_split <- TRUE
    # input$gene_expr_norm <- "cpm"
    ## END TEST ##
    
    # Subset sample metadata
    sample_metadata_subset.dt <- sample_metadata.dt %>%
      .[cell_line%in%input$gene_expr_donors & time%in%input$gene_expr_timepoints]
    
    # Fetch gene expression data
    expr.dt <- data.table(
      sample = colnames(expr.sce),
      expr = assay(expr.sce,input$gene_expr_norm)[input$gene_expr_gene,]
    ) %>% merge(sample_metadata_subset.dt, by="sample")
    
    # Plot all donors

    if (input$gene_expr_norm=="cpm") {
      ylim.max <- max(15,max(expr.dt$expr))
    } else if (input$gene_expr_norm=="log_cpm") {  
      ylim.max <- max(5,max(expr.dt$expr))
    }

    to.plot.means <- expr.dt[,.(expr=mean(expr), sd=sd(expr)), by=c("condition","time","surface_marker","condition_surface_marker","condition_surface_marker_time")]
    
    p1 <- ggplot(expr.dt, aes(x=time, y=expr, fill=time)) +
      geom_bar_interactive(aes(tooltip=condition_surface_marker_time, data_id=condition_surface_marker_time), stat="identity", alpha=0.75, color="black", data=to.plot.means) +
      # geom_point_interactive(aes(tooltip=cell_line, data_id=cell_line), size=2.75, shape=21, stroke=0.2, color="black", alpha=0.90) +
      geom_point(size=2.75, shape=21, stroke=0.2, color="black", alpha=0.90) +
      geom_errorbar(aes(ymin=expr-sd, ymax=expr+sd), color="black", width=0.4, alpha=0.9, linewidth=0.5, data=to.plot.means) +
      scale_fill_manual(values=timepoints_colours) +
      facet_grid(~condition_surface_marker, scales = "free_x", space = "free_x") +
      guides(x = guide_axis(angle = 90)) +
      coord_cartesian(ylim=c(0,ylim.max)) +
      labs(x="", y=sprintf("Gene expression (%s)",input$gene_expr_norm), title=sprintf("%s (all donors)",input$gene_expr_gene)) +
      theme_classic() +
      theme(
        strip.text = element_text(size=rel(0.90)),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(colour="black",size=rel(1.05)),
        axis.text.y = element_text(colour="black",size=rel(1.10)),
        axis.title.y = element_text(colour="black",size=rel(1.15)),
        legend.position = "none"
      )
    
    # Plot per donor
    if (input$gene_expr_split) {
      
      p_list <- list()
      for (j in input$gene_expr_donors) {
        
        to.plot <- expr.dt %>% .[cell_line==j]
        
        to.plot.means <- to.plot[,.(expr=mean(expr), sd=sd(expr)), by=c("condition","time","surface_marker","condition_surface_marker","condition_surface_marker_time")]
        
        p_list[[j]] <- ggplot(to.plot, aes(x=time, y=expr, fill=time)) +
          geom_bar_interactive(aes(tooltip=condition_surface_marker_time, data_id=condition_surface_marker_time), stat="identity", alpha=0.75, color="black", data=to.plot.means) +
          geom_point(size=2.25, shape=21, stroke=0.2, color="black", alpha=0.90) +
          # geom_point_interactive(aes(tooltip=cell_line, data_id=cell_line), size=2.25, shape=21, stroke=0.2, color="black", alpha=0.90) +
          geom_point(size=3, shape=21, stroke=0.2, color="black") +
          scale_fill_manual(values=timepoints_colours) +
          facet_grid(~condition_surface_marker, scales = "free_x", space = "free_x") +
          guides(x = guide_axis(angle = 90)) +
          coord_cartesian(ylim=c(0,ylim.max)) +
          # labs(x="", y=input$gene_expr_norm, title=sprintf("%s (%s)",input$gene_expr_gene,j)) +
          labs(x="", y="", title=sprintf("%s (%s)",input$gene_expr_gene,j)) +
          theme_classic() +
          theme(
            strip.text = element_text(size=rel(0.70)),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(colour="black",size=rel(0.80)),
            axis.text.y = element_text(colour="black",size=rel(0.90)),
            legend.position = "none"
        )
      }
      p2 <- cowplot::plot_grid(plotlist=p_list, ncol=1)
      p <- cowplot::plot_grid(plotlist=list(p1,p2), rel_widths = c(0.6,0.4), ncol=2)
      width_svg <- 14
    } else {
      p <- p1
      width_svg <- 9
    }
  
    girafe(
      code = print(p),
      width_svg = width_svg, height_svg = 7.5,
      options = list(
        opts_selection(type = "single", css = ""),
        opts_hover_inv(css = "opacity:0.25"),
        opts_hover(css = "")
      )
    ) %>% return(.)
    
  })
  
  output$gene_expr_plot = renderGirafe({
    shiny::validate(need(input$gene_expr_gene%in%genes,"Gene not found"))
    plot_gene_expr()
  })
  
  
  ##############################
  ## Dimensionality reduction ##
  ##############################
  
  updateSelectizeInput(session = session, inputId = 'dimred_gene', choices = genes, server = TRUE, selected = "CD24")
  
  plot_dimensionality_reduction <- reactive({
    
    ## START TEST ##
    # input <- list()
    # input$dimred_timepoints <- timepoints
    # input$dimred_donors <- donors
    # input$dimred_layout <- "pca"
    # input$dimred_colourby <- "gene_expression"
    # # input$dimred_colourby <- "donor"
    # input$dimred_gene <- "KLF4"
    # input$dimred_layout_dim1 <- "PC1"
    # input$dimred_layout_dim2 <- "PC2"
    ## END TEST ##
    
    # Subset samples
    sample_metadata_subset.dt <- sample_metadata.dt %>%
      .[cell_line%in%input$dimred_donors & time%in%input$dimred_timepoints]
    sce.filt <- expr.sce[,sample_metadata_subset.dt$sample]
    
    # Feature selection
    tmp <- apply(assay(sce.filt),1,var)
    sce.filt <- sce.filt[tmp>0.25]
    
    # PCA
    set.seed(42)
    sce.filt <- runPCA(sce.filt, ncomponents=5, exprs_values="log_cpm")
    
    if (input$dimred_layout == "pca") {
      dimred.dt <- reducedDim(sce.filt,"PCA")[,c(input$dimred_layout_dim1,input$dimred_layout_dim2)] %>% 
        as.data.table(keep.rownames = T) %>% setnames(c("sample","x","y")) %>% setkey(sample)
      var.explained <- attr(reducedDims(sce.filt)[["PCA"]],"percentVar")
      var.explained <- 100*(var.explained / sum(var.explained))
      labs <- c(
        sprintf("%s (%.2f%%)",input$dimred_layout_dim1,var.explained[as.numeric(gsub("\\D", "", input$dimred_layout_dim1))]), 
        sprintf("%s (%.2f%%)",input$dimred_layout_dim2,var.explained[as.numeric(gsub("\\D", "", input$dimred_layout_dim2))])
      )
    } else {
      stop("Not implemented")
    }
    
    to.plot <- dimred.dt %>% merge(sample_metadata.dt, by="sample")
    
    if (input$dimred_colourby == "gene_expression") {
      to.plot %>% .[,color:=assay(sce.filt,"log_cpm")[input$dimred_gene,sample]] %>% setorder(color)
    } else {
      to.plot$color <- to.plot[[input$dimred_colourby]]
    }
    
    # Plot dimensionality reduction
    p <- ggplot(to.plot, aes(x=x, y=y, fill=color)) +
      geom_point_interactive(aes(tooltip=condition_surface_marker_time, data_id=condition_surface_marker_time), size=7.5, shape=21, alpha=0.9) +
      geom_text_repel(data=to.plot, aes(label=condition_surface_marker_time), size=3.15, segment.size=0.1, max.overlaps = Inf) +
      labs(x=labs[1], y=labs[2]) +
      theme_classic(base_family = "Helvetica") +
      theme(
        axis.title = element_text(color="black", size=rel(1.35)),
        axis.text = element_text(color="black", size=rel(1.10)),
        legend.text = element_text(color="black", size=rel(1.1)),
        legend.title = element_blank()
      )
    
    # Modify legends
    if (input$dimred_colourby=="time") {
      p <- p + scale_fill_manual(values=timepoints_colours)
    } else if (input$dimred_colourby=="donor") {
      # p <- p + scale_color_manual(values=donors_colours)
    } else if (input$dimred_colourby=="cell_cycle_phase") {
      # p <- p + 
      #   # guides(fill = guide_legend(override.aes = list(size=5, alpha = 1))) +
      #   theme(
      #     legend.position = "none",
      #   )
    } else if (input$dimred_colourby=="gene_expression") {
      p <- p + 
        # gene_expr_palette_fill +
        scale_fill_distiller(palette = "YlOrRd", direction=1) +
        labs(title=input$dimred_gene) + 
        theme(
          plot.title = element_text(hjust=0.5)
        )
    }
    
    # Define palette
    # palette <- switch(input$dimred_colourby, "factor" = factor_palette_fill, "gene_expression" = gene_expr_palette_fill)
    # p <- p + palette
    
    # girafe    
    girafe(
      code = print(p),
      width_svg = 10, height_svg = 8.5,
      options = list( 
        opts_sizing(rescale = FALSE),
        opts_selection(type = "single", css = ""),
        opts_hover(css = "cursor:pointer;fill:magenta;color:magenta")
      )
    ) %>% return(.)
    
  })
  
  output$dimred_plot = renderGirafe({
    shiny::validate(need(input$dimred_gene%in%genes, "Select gene..." ))
    plot_dimensionality_reduction()
  })
  
  ##############################
  ## gene expression heatmaps ##
  ##############################

  plot_gene_expr_heatmap <- reactive({

    ## START TEST ##
    # input$gene_expr_heatmap_samples <- "Rich media"
    # input$gene_expr_heatmap_genes <- "ENPP5,DPP4,TXNIP,TNXB,FEN1,CD274,H3C6,COL10A1"
    # input$gene_expr_heatmap_include_hgps <- FALSE
    # input$gene_expr_heatmap_cluster_rows <- F
    # input$gene_expr_heatmap_cluster_cols <- F
    # input$gene_expr_heatmap_fontsize_row <- 10
    # input$gene_expr_heatmap_display_numbers <- FALSE
    # input$gene_expr_heatmap_scale_rows <- TRUE
    # input$gene_expr_heatmap_group_samples <- TRUE
    ## END TEST ##

    # Process gene names
    genes.to.plot <- gsub("[^[:alnum:],]", "", input$gene_expr_heatmap_genes)
    genes.to.plot <- strsplit(genes.to.plot,split = ",")[[1]]
    
    if (input$gene_expr_heatmap_group_samples) {
      if (input$gene_expr_heatmap_samples=="Rich media") {
        rna.mtx <- assay(readRDS("data/aggregate_samples/SingleCellExperiment_rich.rds"),"log_cpm")[genes.to.plot,]
      } else if (input$gene_expr_heatmap_samples=="Base media") {
        rna.mtx <- assay(readRDS("data/aggregate_samples/SingleCellExperiment_base.rds"),"log_cpm")[genes.to.plot,]
      }
      angle_col <- 0
    } else {
      
      # Subset sample metadata
      sample_metadata_subset.dt <- sample_metadata.dt[media==input$gene_expr_heatmap_samples]
      
      # if (!input$gene_expr_heatmap_include_hgps) {
      #   sample_metadata_subset.dt <- sample_metadata_subset.dt[hgps==FALSE]
      # }
      # if (!input$gene_expr_include_longitudinal) {
      #   sample_metadata_subset.dt <- sample_metadata_subset.dt[longitudinal==FALSE]
      # }
      
      # Fetch gene expression
      samples.order <- sample_metadata_subset.dt %>% setorder(-horvath2_clock) %>% .$sample
      rna.mtx <- assay(expr.sce[genes.to.plot,samples.order],"log_cpm")
      colnames(rna.mtx) <- sample_metadata_subset.dt$donor
      
      angle_col <- 90
    }
    
    if (input$gene_expr_heatmap_scale_rows) {
      
      pheatmap(
        mat = rna.mtx,
        color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
        cluster_cols = input$gene_expr_heatmap_cluster_cols, cluster_rows = input$gene_expr_heatmap_cluster_rows,
        show_colnames = T, show_rownames = T,
        legend = T, annotation_legend = F,
        fontsize_row = input$gene_expr_heatmap_fontsize_row, fontsize_col = 12,
        treeheight_row = 20, treeheight_col = 20,
        display_numbers = input$gene_expr_heatmap_display_numbers,
        scale = "row", angle_col=angle_col
      )
      
    } else {
      
      breaks <- seq(from=min(abs(rna.mtx)), to=max(abs(rna.mtx)), by=0.5)
      
      pheatmap(
        mat = rna.mtx,
        breaks = breaks,
        color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaks)),
        cluster_cols = input$gene_expr_heatmap_cluster_cols, cluster_rows = input$gene_expr_heatmap_cluster_rows,
        show_colnames = T, show_rownames = T,
        legend = T, annotation_legend = F,
        fontsize_row = input$gene_expr_heatmap_fontsize_row, fontsize_col = 12,
        treeheight_row = 20, treeheight_col = 20,
        angle_col=angle_col,
        display_numbers = input$gene_expr_heatmap_display_numbers
      )
      
    }
  })

  output$gene_expr_heatmap = renderPlot({
    tmp <- gsub("[^[:alnum:],]", "", input$gene_expr_heatmap_genes); tmp <- strsplit(tmp,split = ",")[[1]]
    genes.not.found <- tmp[!tmp%in%genes]
    shiny::validate(need(length(genes.not.found)==0, sprintf("The following genes are not found:\n%s", paste(genes.not.found,collapse="\n"))))
    shiny::validate(need(length(tmp)>=2, "A heatmap requires at least 2 genes"))
    plot_gene_expr_heatmap()
  })
  

  ##################################
  ## Gene Set Enrichment Analysis ##
  ##################################
  
  plot_gsea <- reactive({

    ## START TEST ##
    # input$gsea_db <- "EnrichR_Reactome_2022"
    # input$gsea_groupA <- "Dox_4ON_CD13"
    # input$gsea_groupB <- "Dox_4ON_CD70"
    # input$gsea_num_pathways <- 15
    ## END TEST ##

    # Fetch data
    if (input$gsea_db=="GO") {
      gsea_list <- readRDS(file.path(data_folder,sprintf("gsea/gsea_GO_%s_vs_%s.rds",input$gsea_groupA,input$gsea_groupB)))
    } else if (input$gsea_db=="MSigDB_Hallmarks") {
      gsea_list <- readRDS(file.path(data_folder,sprintf("gsea/gsea_msigdb_%s_vs_%s.rds",input$gsea_groupA,input$gsea_groupB)))
    } else if (input$gsea_db=="MitoCarta") {
      gsea_list <- readRDS(file.path(data_folder,sprintf("gsea/gsea_mitocarta_%s_vs_%s.rds",input$gsea_groupA,input$gsea_groupB)))
    # } else if (input$gsea_db=="Reactome") {
    #   gsea_list <- readRDS(file.path(data_folder,"gsea/gsea_reactome_list.rds"))
    } else if (input$gsea_db=="EnrichR_GO_Molecular_Function_2023") {
      gsea_list <- readRDS(file.path(data_folder,sprintf("gsea/gsea_enrichr_GO_Molecular_Function_2023_%s_vs_%s.rds",input$gsea_groupA,input$gsea_groupB)))
    } else if (input$gsea_db=="EnrichR_GO_Cellular_Component_2023") {
      gsea_list <- readRDS(file.path(data_folder,sprintf("gsea/gsea_enrichr_GO_Cellular_Component_2023_%s_vs_%s.rds",input$gsea_groupA,input$gsea_groupB)))
    } else if (input$gsea_db=="EnrichR_GO_Biological_Process_2023") {
      gsea_list <- readRDS(file.path(data_folder,sprintf("gsea/gsea_enrichr_GO_Biological_Process_2023_%s_vs_%s.rds",input$gsea_groupA,input$gsea_groupB)))
    } else if (input$gsea_db=="EnrichR_Reactome_2022") {
      gsea_list <- readRDS(file.path(data_folder,sprintf("gsea/gsea_enrichr_Reactome_2022_%s_vs_%s.rds",input$gsea_groupA,input$gsea_groupB)))
    } else if (input$gsea_db=="EnrichR_MSigDB_Hallmark_2020") {
      gsea_list <- readRDS(file.path(data_folder,sprintf("gsea/gsea_enrichr_MSigDB_Hallmark_2020_%s_vs_%s.rds",input$gsea_groupA,input$gsea_groupB)))
    }

    # Plot
    p_list <- list()

    for (k in c("down","up")) {

      tmp <- gsea_list[[k]]
      if (is.null(tmp)) {
        to.plot <- data.table(Description="NA", Count=1, qvalue=1)
      } else {
        to.plot <- as.data.table(tmp) %>% head(n=input$gsea_num_pathways) %>% setorder(-qvalue)
        
        to.plot %>%
          .[,log_pval:=-log10(qvalue)] %>%
          .[,Description:=substr(Description,1,70)] %>%
          .[,Description:=factor(Description, levels=Description)]
      }

      p_list[[k]] <- ggplot(to.plot, aes(x=log_pval, y=Description)) +
        geom_point(aes(fill=log_pval), size=4, shape=21) +
        scale_fill_gradient(low = "gray80", high = "red") +
        coord_cartesian(xlim=c(0,max(10,max(to.plot$log_pval)))) +
        labs(x="Enrichment score (log10 q-value)", y="", title=sprintf("%s-regulated in %s\n(background: %s)",k,input$gsea_groupB,input$gsea_groupA)) +
        theme_bw() +
        theme(
          plot.title = element_text(hjust=0.5, size=rel(1.20)),
          axis.text.x = element_text(color="black", size=rel(1.1)),
          axis.text.y = element_text(color="black", size=rel(1.0)),
          axis.title = element_text(color="black", size=rel(1.15)),
          legend.position = "none"
        )

    }
    p <- cowplot::plot_grid(plotlist=p_list, nrow=1)

    return(p)
  })

  output$gsea_plot = renderPlot({
    shiny::validate(need(input$gsea_db%in%databases_gsea, "GSEA Database not recognised"))
    plot_gsea()
  })
  
  # output$gsea_download <- downloadHandler(
  #   filename = function() {
  #     sprintf("gsea_%s_list.rds",input$gsea_db)
  #   },
  #   content = function(file) {
  #     diff.dt <- readRDS(file.path(data_folder,sprintf("limma/diff_results_%s.rds",input$diff_samples)))[[input$diff_covariate]] %>%
  #       .[feature%in%genes] %>% setnames("pval","pvalue")
  #     fwrite(diff.dt, file, quote = F, na="NA", sep="\t")
  #   },
  #   contentType = "text/plain"
  # )
  
  ##################################
  ## Differential gene expression ##
  ##################################
  
  updateSelectizeInput(session = session, inputId = 'diff_gene', choices = genes, server = TRUE, selected = "CD70")
  
  # # Update gene upon click
  # observeEvent(input$diff_plot_selected, {
  #   updateSelectizeInput(session = session, inputId = 'diff_gene', choices = genes, server = TRUE, selected = input$diff_plot_selected)
  # })
  # 
  # observeEvent(input$diff_dt_row_last_clicked, {
  #   file <- file.path(data_folder,sprintf("limma/pairwise_all/diff_%s_vs_%s.csv.gz",input$diff_groupA,input$diff_groupB))
  #   diff.dt <- fread(file) %>% .[gene%in%genes]
  #   updateSelectizeInput(session = session, inputId = 'diff_gene', choices = genes, server = TRUE, selected = diff.dt[input$diff_dt_row_last_clicked]$gene)
  # })

  load_diff_data <- reactive({
    
    ## START TEST ##
    # input <- list()
    # input$diff_groupA <- "Dox_4ON_CD13"
    # input$diff_groupB <- "Dox_4ON_CD70"
    # input$diff_gene_class <- "All"
    # input$diff_gene <- "CD70"
    # input$diff_range <- c(-10,10)
    # input$diff_label_pos <- 15
    # input$diff_label_neg <- 15
    ## END TEST ##
  
    req(input$diff_groupA, input$diff_groupB, input$diff_gene)
    
    ## Volcano plot ##
    
    # Fetch data
    file <- file.path(data_folder,sprintf("differential/diff_%s_vs_%s.txt.gz",input$diff_groupA,input$diff_groupB))
    
    diff.dt <- fread(file) %>% .[feature%in%genes] %>% 
      setnames("p.value","pvalue")
    
    diff.dt %>% .[logFC>10,logFC:=10] %>% .[logFC<(-10),logFC:=(-10)]
    
    # Parse data
    to.plot <- diff.dt %>% 
      .[logFC>=input$diff_range[1] & logFC<=input$diff_range[2]] %>% 
      .[is.na(logFC),logFC:=0] %>%
      .[is.na(pvalue),pvalue:=1] %>%
      .[,log_pval:=-log10(pvalue+1e-100)]
    
    # Subset features
    if (input$diff_gene_class=="Surface markers") {
      to.plot <- to.plot[feature%in%gene_sets[["surface_markers"]]]
    } else if (input$diff_gene_class=="Transcription Factors") {
      to.plot <- to.plot[feature%in%gene_sets[["transcription_factors"]]]
    } else if (input$diff_gene_class=="DNA repair") {
      to.plot <- to.plot[feature%in%gene_sets[["dna_repair"]]]
    } else if (input$diff_gene_class=="Epigenetic regulators") {
      to.plot <- to.plot[feature%in%gene_sets[["epigenetic_regulators"]]]
    } else if (input$diff_gene_class=="Cell cycle genes") {
      to.plot <- to.plot[feature%in%gene_sets[["cell_cycle"]]]
    } else if (input$diff_gene_class=="ECM genes") {
      to.plot <- to.plot[feature%in%gene_sets[["ecm"]]]
    } else if (input$diff_gene_class=="Mitochondrial") {
      to.plot <- to.plot[feature%in%gene_sets[["mitochondrial"]]]
    } else if (input$diff_gene_class=="Ribosomal") {
      to.plot <- to.plot[feature%in%gene_sets[["ribosomal"]]]
    }
    return(to.plot)
  })
  
  plot_diff <- reactive({
    
    # Fetch diff output
    to.plot <- load_diff_data()
    
    # Define plot margins
    xlim_min <- input$diff_range[1]; xlim_max <- input$diff_range[2]; margin_text <- 2; margin_dots <- 0.1
    to.plot[logFC<=xlim_min,logFC:=xlim_min]; to.plot[logFC>=xlim_max,logFC:=xlim_max]
    
    # to.plot.subset <- to.plot[log_pval>=input$min_log_pval]
    to.plot.subset <- to.plot[log_pval>0]
    
    if (!input$diff_gene%in%to.plot.subset$feature) {
      to.plot.subset <- rbind(to.plot.subset,to.plot[feature==input$diff_gene])
    }
    
    ylim_min <- min(to.plot.subset$log_pval,na.rm=T); ylim_max <- max(to.plot$log_pval,na.rm=T)
    
    p.volcano <- ggplot(to.plot.subset, aes(x=logFC, y=log_pval)) +
      # geom_segment(x=0, xend=0, y=input$min_log_pval, yend=ylim_max, color="orange", linewidth=0.25) +
      geom_segment(x=0, xend=0, y=0, yend=ylim_max, color="orange", linewidth=0.25) +
      geom_point_interactive(aes(color=log_pval, tooltip=feature, data_id=feature, alpha=log_pval, size=log_pval, onclick=feature)) + 
      geom_point_interactive(aes(tooltip=feature, data_id=feature, onclick=feature), size=2.5, fill="green", shape=21, data=to.plot.subset[feature==input$diff_gene]) + 
      geom_hline(yintercept=input$min_log_pval, linetype="dashed") +
      ggrepel::geom_text_repel(aes(label=feature), size=4, data=to.plot[feature==input$diff_gene]) +
      ggrepel::geom_text_repel(aes(label=feature), size=3, data=head(to.plot[logFC<0 & feature!=input$diff_gene],n=input$diff_label_neg), segment.size=0.02, max.overlaps = Inf) +
      ggrepel::geom_text_repel(aes(label=feature), size=3, data=head(to.plot[logFC>0 & feature!=input$diff_gene],n=input$diff_label_pos), segment.size=0.02, max.overlaps = Inf) +
      scale_color_gradient(low = "gray80", high = "black") +
      scale_alpha_continuous(range=c(0.25,1)) +
      scale_size_continuous(range=c(0.15,2)) +
      scale_x_continuous(limits=c(xlim_min-margin_dots,xlim_max+margin_dots)) +
      scale_y_continuous(limits=c(ylim_min,ylim_max+1)) +
      annotate("text", x=0, y=ylim_max+ylim_max/12, size=5, label=sprintf("(%d)", nrow(to.plot.subset))) +
      # annotate("text", x=-xlim_max+(xlim_max/2.5), y=ylim_max+ylim_max/12, size=5, label=sprintf("(%d)", )) +
      coord_cartesian(ylim=c(0,ylim_max+ylim_max/10)) +
      # annotate("text", x=xlim_min-0.05, y=ylim_max+2, size=4, label=sprintf("%d (-)",length(negative_hits))) +
      # annotate("text", x=xlim_max+0.05, y=ylim_max+2, size=4, label=sprintf("%d (+)",length(positive_hits))) +
      annotate("text", x=xlim_min+margin_text, y=0.1, size=4.5, label=sprintf("Higher in %s",input$diff_groupA)) +
      annotate("text", x=xlim_max-margin_text, y=0.1, size=4.5, label=sprintf("Higher in %s",input$diff_groupB)) +
      labs(x="Limma coefficient", y=expression(paste("-log"[10],"(p.value)"))) +
      theme_classic(base_family = "Helvetica") +
      theme(
        axis.text = element_text(size=rel(1.3), color='black'),
        axis.title = element_text(size=rel(1.5), color='black'),
        legend.position="none"
      )
    
    # girafe call
    girafe(
      code = print(p.volcano),
      width_svg = 11, height_svg = 9,
      options = list(
        opts_zoom(min = 1, max = 5),
        opts_sizing(rescale = TRUE),
        opts_selection(type = "single", css = ""),
        opts_hover(css = "cursor:pointer")
      )
    ) %>% return(.)
    
  })
  
  output$diff_plot <- renderGirafe({
    # shiny::validate(need(input$diff_groupA%in%factors,"Please select factor A"))
    # shiny::validate(need(input$diff_groupB%in%factors,"Please select factor B"))
    # shiny::validate(need(input$diff_groupA!=input$diff_groupB,"group A and B must be different"))
    shiny::validate(need(input$diff_gene%in%genes,"Please select a gene"))
    plot_diff()
  })
  
  output$diff_download <- downloadHandler(
    filename = function() {
      "diff_output.tsv"
    },
    content = function(file) {
      diff.dt <- load_diff_data()
      fwrite(diff.dt, file, quote = F, na="NA", sep="\t")
    },
    contentType = "text/plain"
  )


  ###################################
  ## Gene signatures (precomputed) ##
  ###################################
  
  updateSelectizeInput(session = session, inputId = 'gs_gene_signature', choices = gene_signatures_names, server = TRUE, selected = "OSKM_regulon_positive")
  
  plot_gene_signatures_precomputed <- reactive({
    
    ## START TEST ##
    # input <- list()
    # input$gs_gene_signature <- "OSKM_regulon_positive"
    # input$gs_timepoints <- timepoints
    # input$gs_donors <- donors
    # input$gs_split <- TRUE
    ## END TEST ##
    
    # Subset sample metadata
    sample_metadata_subset.dt <- sample_metadata.dt %>%
      .[cell_line%in%input$gs_donors & time%in%input$gs_timepoints]
    
    # Fetch UCell data
    ucell.dt <- data.table(
      sample = colnames(ucell.sce),
      expr = assay(ucell.sce)[input$gs_gene_signature,]
    ) %>% merge(sample_metadata_subset.dt, by="sample")
    
    # Plot all donors
    to.plot.means <- ucell.dt[,.(expr=mean(expr), sd=sd(expr)), by=c("condition","time","surface_marker","condition_surface_marker","condition_surface_marker_time")]
    
    p1 <- ggplot(ucell.dt, aes(x=time, y=expr, fill=time)) +
      geom_bar_interactive(aes(tooltip=condition_surface_marker_time, data_id=condition_surface_marker_time), stat="identity", alpha=0.75, color="black", data=to.plot.means) +
      # geom_point_interactive(aes(tooltip=cell_line, data_id=cell_line), size=2.75, shape=21, stroke=0.2, color="black", alpha=0.90) +
      geom_point(size=2.75, shape=21, stroke=0.2, color="black", alpha=0.90) +
      geom_errorbar(aes(ymin=expr-sd, ymax=expr+sd), color="black", width=0.4, alpha=0.9, linewidth=0.5, data=to.plot.means) +
      scale_fill_manual(values=timepoints_colours) +
      facet_grid(~condition_surface_marker, scales = "free_x", space = "free_x") +
      guides(x = guide_axis(angle = 90)) +
      # coord_cartesian(ylim=c(0,ylim.max)) +
      labs(x="", y="Gene Signature (UCell)", title="All donors") +
      theme_classic() +
      theme(
        strip.text = element_text(size=rel(0.90)),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(colour="black",size=rel(1.05)),
        axis.text.y = element_text(colour="black",size=rel(1.10)),
        axis.title.y = element_text(colour="black",size=rel(1.15)),
        legend.position = "none"
      )
    
    # Plot per donor
    if (input$gs_split) {
      
      p_list <- list()
      for (j in input$gs_donors) {
        
        to.plot <- ucell.dt %>% .[cell_line==j]
        
        to.plot.means <- to.plot[,.(expr=mean(expr), sd=sd(expr)), by=c("condition","time","surface_marker","condition_surface_marker","condition_surface_marker_time")]
        
        p_list[[j]] <- ggplot(to.plot, aes(x=time, y=expr, fill=time)) +
          geom_bar_interactive(aes(tooltip=condition_surface_marker_time, data_id=condition_surface_marker_time), stat="identity", alpha=0.75, color="black", data=to.plot.means) +
          geom_point(size=2.25, shape=21, stroke=0.2, color="black", alpha=0.90) +
          # geom_point_interactive(aes(tooltip=cell_line, data_id=cell_line), size=2.25, shape=21, stroke=0.2, color="black", alpha=0.90) +
          geom_point(size=3, shape=21, stroke=0.2, color="black") +
          scale_fill_manual(values=timepoints_colours) +
          facet_grid(~condition_surface_marker, scales = "free_x", space = "free_x") +
          guides(x = guide_axis(angle = 90)) +
          # coord_cartesian(ylim=c(0,ylim.max)) +
          # labs(x="", y=input$gs_norm, title=sprintf("%s (%s)",input$gs_gene,j)) +
          labs(x="", y="", title=sprintf("%s",j)) +
          theme_classic() +
          theme(
            strip.text = element_text(size=rel(0.70)),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(colour="black",size=rel(0.80)),
            axis.text.y = element_text(colour="black",size=rel(0.90)),
            legend.position = "none"
          )
      }
      p2 <- cowplot::plot_grid(plotlist=p_list, ncol=1)
      p <- cowplot::plot_grid(plotlist=list(p1,p2), rel_widths = c(0.6,0.4), ncol=2)
      width_svg <- 14
    } else {
      p <- p1
      width_svg <- 9
    }
    
    girafe(
      code = print(p),
      width_svg = width_svg, height_svg = 7.5,
      options = list(
        opts_selection(type = "single", css = ""),
        opts_hover_inv(css = "opacity:0.25"),
        opts_hover(css = "")
      )
    ) %>% return(.)
  })
    
  output$plot_gs_precomputed = renderGirafe({
    shiny::validate(need(input$gs_gene_signature%in%gene_signatures_names,"Precomputed gene signature not found"))
    plot_gene_signatures_precomputed()
  })
  

#############################
## De novo gene signatures ##
#############################

  # updateSelectizeInput(session = session, inputId = 'gsn_genes', choices = gene_signatures_names, server = TRUE, selected = "msigdb_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
  # 
  # plot_gene_signatures_denovo <- reactive({
  #   
  #   ## START TEST ##
  #   input <- list()
  #   input$gsn_genes <- "CD177,PI16,LY6D,DPP4"
  #   input$gsn_x_axis <- "epigenetic_age"
  #   input$gsn_colourby <- "epigenetic_age"
  #   input$gsn_include_hgps <- TRUE
  #   input$gsn_include_longitudinal <- TRUE
  #   input$gsn_split_media <- TRUE
  #   input$gsn_label_samples <- TRUE
  #   ## END TEST ##
  #   
  #   # Process gene names
  #   genes.to.use <- gsub("[^[:alnum:],]", "", input$gene_sign_genes)
  #   genes.to.use <- strsplit(genes.to.use,split = ",")[[1]]
  #   
  #   # Subset sample metadata
  #   sample_metadata_subset.dt <- sample_metadata.dt
  #   
  #   if (!input$gsn_include_hgps) {
  #     sample_metadata_subset.dt <- sample_metadata_subset.dt[hgps==FALSE]
  #   }
  #   if (!input$gsn_include_longitudinal) {
  #     sample_metadata_subset.dt <- sample_metadata_subset.dt[longitudinal==FALSE]
  #   }
  #   
  #   # Calculate gene signature
  #   sce <- ScoreSignatures_UCell(sce, features=list("tmp"=genes.to.use), maxRank=1e4)
  # 
  #   # Fetch gene signature scores
  #   to.plot <- data.table(
  #     sample = colnames(expr.sce),
  #     expr = assay(ucell.sce)[input$gsn_gene_signature,]
  #   ) %>% merge(sample_metadata_subset.dt, by="sample")
  #   
  #   if (input$gsn_x_axis=="chronological_age") {
  #     x.axis <- "Chronological age"
  #     to.plot$value <- to.plot$age
  #   } else if (input$gsn_x_axis=="epigenetic_age") {
  #     x.axis <- "Epigenetic age (Horvath2)"
  #     to.plot$value <- to.plot$horvath2_clock
  #   }
  #   
  #   # Plot
  #   p1 <- ggplot(to.plot, aes(x=media, y=expr, fill=media)) +
  #     geom_boxplot(color="black") +
  #     geom_point(color="black", shape=21, alpha=0.5, size=2) +
  #     scale_fill_brewer(palette="Dark2") +
  #     labs(x="", y="UCell signature", title="") +
  #     theme_classic() +
  #     theme(
  #       axis.text.y = element_text(colour="black",size=rel(1.25)),
  #       axis.text.x = element_text(colour="black",size=rel(1.25)),
  #       axis.title.y = element_text(colour="black",size=rel(1.20)),
  #       # axis.text.y = element_text(colour="black",size=rel(1.25)),
  #       legend.position = "none"
  #     )
  #   
  #   p <- ggplot(to.plot, aes(x=value, y=expr)) +
  #     stat_smooth(method="lm", color="gray30", alpha=0.05) +
  #     stat_cor(aes(label = after_stat(r.label)),  method = "pearson") +
  #     labs(x=x.axis, y="UCell signature", title=input$gsn_gene_signature) +
  #     theme_classic() +
  #     theme(
  #       plot.title = element_text(hjust = 0.5),
  #       axis.text.y = element_text(colour="black",size=rel(1.25)),
  #       axis.title.x = element_text(colour="black",size=rel(1.20)),
  #       axis.title.y = element_text(colour="black",size=rel(1.20)),
  #       strip.text = element_text(colour="black",size=rel(1.15)),
  #       legend.position = "right",
  #       legend.title = element_blank()
  #     )
  #   
  #   if (input$gsn_colourby%in%c("epigenetic_age","chronological_age")) {
  #     p <- p + geom_point_interactive(aes(fill=value, tooltip=donor, data_id=donor), color="black", alpha=0.90, size=4, shape=21, stroke=0.2) +
  #       scale_fill_distiller(palette = "Blues", direction=1)
  #   } else if (input$gsn_colourby=="hgps") {
  #     p <- p + geom_point_interactive(aes(fill=hgps, tooltip=donor, data_id=donor), color="black", alpha=0.90, size=4, shape=21, stroke=0.2) +
  #       scale_fill_brewer(palette="Dark2")
  #   }
  #   
  #   
  #   if (input$gsn_split_media) { 
  #     p <- p + facet_wrap(~media, nrow=1)
  #   }
  #   
  #   if (input$gsn_label_samples) {
  #     p <- p + ggrepel::geom_text_repel(aes(label=donor), size=2, segment.size=0.1, max.overlaps = Inf)
  #   }
  #   
  #   girafe(
  #     # code = print(p),
  #     code = print(p1+p + plot_layout(ncol=2, widths=c(1.25,4))),
  #     width_svg = 13, height_svg = 5.5,
  #     options = list(
  #       opts_selection(type = "single", css = ""),
  #       opts_hover_inv(css = "opacity:0.25"),
  #       opts_hover(css = "")
  #     )
  #   ) %>% return(.)
  #   
  #   })
  #   
  # output$plot_gene_signature_denovo = renderGirafe({
  #   
  #   tmp <- gsub("[^[:alnum:],]", "", input$gene_sign_genes); tmp <- strsplit(tmp,split = ",")[[1]]
  #   genes.not.found <- tmp[!tmp%in%genes]
  #   shiny::validate(need(length(genes.not.found)==0, sprintf("The following genes are not found:\n%s", paste(genes.not.found,collapse="\n"))))
  #   shiny::validate(need(length(tmp)>=2, "Calculation of gene signatures require at least 2 genes"))
  #   plot_gene_signatures()
  #   
  #   # shiny::validate(need(input$gsn_gene_signature%in%gene_signatures_names,"Precomputed gene signature not found"))
  #   plot_gene_signatures_denovo()
  # })
  
  
  #####################
  ## Gene expression ##
  #####################
  
  # updateSelectizeInput(session = session, inputId = 'gene_expr_gene', choices = genes, server = TRUE, selected = "CD70")
  
  plot_transcriptional_clock <- reactive({
    
    ## START TEST ##
    # input <- list()
    # input$tclock_clock <- "fibrosis"
    # input$tclock_timepoints <- timepoints
    # input$tclock_donors <- donors
    # input$tclock_split <- TRUE
    ## END TEST ##
    
    # Subset sample metadata
    sample_metadata_subset.dt <- sample_metadata.dt %>%
      .[cell_line%in%input$tclock_donors & time%in%input$tclock_timepoints]
    
    # Fetch gene expression data
    tmp <- transcriptional_clocks.dt[,c("sample",input$tclock_clock),with=F] %>%
      setnames(c("sample","value")) %>%
      merge(sample_metadata_subset.dt, by="sample")
    
    # Plot per donor
    dt_list <- list()
    p_list <- list()
    for (j in cell_lines) {
      
      to.plot <- tmp %>% .[cell_line==j]
      
      dt_list[[j]] <- to.plot %>%
        .[,ctrl_value:=to.plot[cell_line==j & group=="Ctrl_OFF",value]] %>%
        .[,value_diff:=value-ctrl_value] %>%
        .[!group=="Ctrl_OFF"]
      
      # to.plot.means <- dt_list[[j]][,.(value_diff=mean(value_diff)), by=c("condition","time","surface_marker","condition_surface_marker")]
      to.plot.means <- to.plot[,.(value=mean(value)), by=c("condition","time","surface_marker","condition_surface_marker","condition_surface_marker_time")]
      
      p_list[[j]] <- ggplot(to.plot.means, aes(x=condition_surface_marker, y=value, fill=condition_surface_marker)) +
        # geom_bar(stat="identity", alpha=0.75, color="black", width=0.80) +
        geom_bar_interactive(aes(tooltip=condition_surface_marker_time, data_id=condition_surface_marker_time), stat="identity", alpha=0.75, color="black", data=to.plot.means) +
        scale_fill_manual(values=conditions_colours) +
        facet_grid(~time, scales = "free_x", space = "free_x") +
        guides(x = guide_axis(angle = 90)) +
        labs(x="", y="Transcriptome age", title=j) +
        theme_classic() +
        theme(
          strip.text = element_text(size=rel(0.90)),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(colour="black",size=rel(0.90)),
          axis.text.y = element_text(colour="black",size=rel(1)),
          axis.title.y = element_text(colour="black",size=rel(0.95)),
          legend.position = "none"
        )
    }
    p2 <- cowplot::plot_grid(plotlist=p_list, ncol=1)
    
    # Plot all donors
    to.plot.all <- rbindlist(dt_list)
    to.plot.means <- to.plot.all %>% .[,.(value_diff=mean(value_diff), sd=sd(value_diff)), by=c("condition","time","surface_marker","condition_surface_marker")]
    
    p1 <- ggplot(to.plot.all, aes(x=time, y=value_diff, fill=time)) +
      # geom_point(size=4, shape=21, stroke=0.2, color="black", alpha=0.90) +
      geom_point_interactive(aes(tooltip=cell_line, data_id=cell_line), size=4, shape=21, stroke=0.2, color="black", alpha=0.90) +
      geom_bar(stat="identity", alpha=0.75, color="black", data=to.plot.means) +
      # geom_bar_interactive(aes(tooltip=condition_surface_marker_time, data_id=condition_surface_marker_time), stat="identity", alpha=0.75, color="black", data=to.plot.means) +
      geom_errorbar(aes(ymin=value_diff-sd, ymax=value_diff+sd), color="black", width=0.4, alpha=0.9, linewidth=0.5, data=to.plot.means) +
      scale_fill_manual(values=timepoints_colours) +
      facet_grid(~condition_surface_marker, scales = "free_x", space = "free_x") +
      guides(x = guide_axis(angle = 90)) +
      # coord_cartesian(ylim=c(0,ylim.max)) +
      labs(x="", y="Transcriptomic age difference (relative control)", title="") +
      theme_classic() +
      theme(
        strip.text = element_text(size=rel(1.05)),
        plot.title = element_text(hjust=0.5, size=rel(1.25), color="black"),
        axis.text.x = element_text(colour="black",size=rel(1)),
        axis.text.y = element_text(colour="black",size=rel(1.15)),
        axis.title.y = element_text(colour="black",size=rel(1.15)),
        legend.position = "none"
      )
    
    if (input$tclock_split) {
      p <- p1
    } else {
      p <- p2
    }
      
    girafe(
      code = print(p),
      width_svg = 9.5, height_svg = 7.5,
      options = list(
        opts_selection(type = "single", css = ""),
        opts_hover_inv(css = "opacity:0.25"),
        opts_hover(css = "")
      )
    ) %>% return(.)
    
  })
  
  output$plot_transcriptional_clock = renderGirafe({
    # shiny::validate(need(input$tclock_gene%in%genes,"Gene not found"))
    plot_transcriptional_clock()
  })
  
  
  
}