
##################
## Shiny app UI ##
##################

ui <- shinyUI(
  fluidPage(
    navbarPage(
        title = "Gene expression of transient partial OSKM reprogramming of Human Dermal Fibroblasts (Reik Lab)",
        theme = shinytheme("spacelab"),
        
        tabPanel(
          title = "Overview", id = "overview",
          br(),
          includeMarkdown("overview.md")
        ),


        #####################
        ## Gene expression ##
        #####################
        
        tabPanel(
          title = "Gene expression", id = "gene_expr",
          sidebarPanel(width=3,
                       selectizeInput("gene_expr_gene", "Select gene", choices = NULL),
                       selectInput("gene_expr_timepoints", "Select timepoints", choices=timepoints, selected=timepoints, multiple=TRUE),
                       selectInput("gene_expr_donors", "Select donors", choices=donors, selected=donors, multiple=TRUE),
                       checkboxInput("gene_expr_split", "Split donors?", value=TRUE),
                       selectizeInput("gene_expr_norm", "Select Normalisation", choices = c("CPM"="cpm", "log CPM"="log_cpm"), selected="cpm")
          ),
          mainPanel(
            girafeOutput("gene_expr_plot", width = "1100px", height = "800px")
          )
        ),
        
        ##############################
        ## Dimensionality reduction ##
        ##############################
        
        tabPanel(
          title = "Dimensionality Reduction", id = "dimred",
          sidebarPanel(width=3,
                       selectizeInput("dimred_timepoints", "Select timepoints", choices=timepoints, selected=timepoints, multiple=TRUE),
                       selectizeInput("dimred_donors", "Select donors", choices=donors, selected=c("63F"), multiple=TRUE),
                       selectInput(inputId = "dimred_layout", label = "Layout", choices = c("PCA" = "pca"), selected = "pca"),
                       conditionalPanel(
                         condition = "input.dimred_layout == 'pca'",

                         selectizeInput("dimred_layout_dim1", "Select x-axis:", choices = paste0("PC",1:5), selected = "PC1"),
                         selectizeInput("dimred_layout_dim2", "Select y-axis:", choices = paste0("PC",1:5), selected = "PC3")
                       ),

                       # selectInput(inputId = "dimred_sample", label = "Sample", choices = samples.choices, selected = samples, multiple = TRUE),
                       selectInput(inputId = "dimred_colourby", label = "Plot colour", choices = c("Donor"="donor", "Timepoint"="time", "Cell cycle phase"="cell_cycle_phase" ,"Gene expression"="gene_expression"), selected = "time"),
                       # sliderInput("dimred_dot_size", label = "Minimum log pvalue", min=0, max=5, step=0.5, value = 3),
                       conditionalPanel(
                         condition = "input.dimred_colourby == 'gene_expression'",
                         selectizeInput("dimred_gene", "Select gene", choices = NULL, selected = "CD70")
                       ),
          ),
          mainPanel(
            girafeOutput("dimred_plot")
          )
        ),
        
        ##############################
        ## Gene expression heatmaps ##
        ##############################

        # tabPanel(
        #   title = "Gene expression heatmaps", id = "gene_expr_heatmap",
        #   sidebarPanel2(width=3,
        #                selectInput(inputId = "gene_expr_heatmap_samples", label = "Select samples", choices = c("Base media","Rich media"), selected = "Rich media"),
        #                textInput("gene_expr_heatmap_genes", label = "Comma-separated list of genes", value = "ENPP5,DPP4,TXNIP,TNXB,FEN1,CD274,H3C6,COL10A1"),
        #                # checkboxInput("gene_expr_heatmap_include_hgps", "Include HGPS samples?", value = FALSE),
        #                checkboxInput("gene_expr_heatmap_cluster_rows", "Cluster rows?", value = TRUE),
        #                checkboxInput("gene_expr_heatmap_cluster_cols", "Cluster columns?", value = FALSE),
        #                checkboxInput("gene_expr_heatmap_scale_rows", "Scale rows (z-scores)?", value = TRUE),
        #                checkboxInput("gene_expr_heatmap_display_numbers", "Display numbers?", value = FALSE),
        #                checkboxInput("gene_expr_heatmap_group_samples", "Group samples in discrete age groups?", value = FALSE),
        #                sliderInput("gene_expr_heatmap_fontsize_row", label = "Font size (rows)", min=5, max=25, step=1, value = 12),
        #     out=tags$p("Samples are sorted by epigenetic age (Horvath2)")
        #   ),
        #   mainPanel(
        #     plotOutput("gene_expr_heatmap", width = "900px", height = "600px")
        #   )
        # ),
        
        #############################
        ## Differential expression ##
        #############################
        
        tabPanel(
          title = "Differential expression", id = "diff",
          sidebarPanel(width=3,
            selectInput(inputId = "diff_groupA", label = "Select group A", choices = groups, selected = "Dox_4ON_CD13"),
            selectInput(inputId = "diff_groupB", label = "Select group B", choices = groups, selected = "Dox_4ON_CD70"),
            selectInput(inputId = "diff_gene_class", label = "Subset gene class", choices = c("All",gene_sets_names), selected = "All"),
            selectizeInput("diff_gene", "Select gene", choices = NULL),
            # sliderInput("min_log_pval", label = "Minimum log pvalue", min=0, max=10, step=0.25, value = 0),
            sliderInput("diff_range", label = "Range of coefficients", min=-10, max=10, step=0.5, value = c(-10,10)),
            sliderInput("diff_label_pos", label = "Label top hits (+)", min=10, max=50, step=1, value = 25),
            sliderInput("diff_label_neg", label = "Label top hits (-)", min=10, max=50, step=1, value = 25),
            downloadButton("diff_download", label = "Download file")
          ),
          mainPanel(
            girafeOutput("diff_plot", width = "1000px", height = "650px")
            # br(),
            # DT::dataTableOutput('limma_dt', width = "100%", height = "auto", fill = FALSE)
          )
        ),
        
        ##################################
        ## Gene Set Enrichment Analysis ##
        ##################################
        
        tabPanel(
          title = "Gene Set Enrichment", id = "gsea",
          sidebarPanel2(width=3,
                        selectInput(inputId = "gsea_groupA", label = "Select group A", choices = groups_no_ctrl, selected = "Dox_4ON_CD13"),
                        selectInput(inputId = "gsea_groupB", label = "Select group B", choices = groups_no_ctrl, selected = "Dox_4ON_CD70"),
                        selectInput(inputId = "gsea_db", label = "Select database", choices = databases_gsea, selected = "MSigDB_Hallmarks"),
                        sliderInput("gsea_num_pathways", label = "Number of pathways", min=5, max=50, step=1, value = 15),
                        out = h5("GSEA pathway scores are derived from the coefficients of the limma model of gene expression versus epigenetic age")
          ),
          mainPanel(
            plotOutput("gsea_plot", width = "1000px", height = "500px"),
          )
        ),


        #################################
        ## Precomputed gene signatures ##
        #################################
        
        tabPanel(
          title = "Precomputed gene signatures", id = "gs_precomputed",
          sidebarPanel2(width=3,
                        selectizeInput("gs_gene_signature", "Select gene signature", choices = NULL),
                        selectInput("gs_timepoints", "Select timepoints", choices=timepoints, selected=timepoints, multiple=TRUE),
                        selectInput("gs_donors", "Select donors", choices=donors, selected=donors, multiple=TRUE),
                        checkboxInput("gs_split", "Split donors?", value=TRUE),
                        out = h5("Gene signatures are calculated using the UCell package. Scores range from 0 to 1, as a non-rigurous rule of thumb any difference <0.1 is likely not significant.")
          ),
          mainPanel(
            # girafeOutput("plot_gene_signature_precomputed", width = "1000px", height = "400px")
            girafeOutput("plot_gs_precomputed", width = "1100px", height = "800px")
          )
        ),

        ############################
        ## Transcriptional clocks ##
        ############################
        
        tabPanel(
          title = "Transcriptional clocks", id = "transcriptional_clock",
          sidebarPanel2(width=3,
                        selectizeInput("tclock_clock", "Select clock", choices = transcriptional_clocks_names, selected="Multi_tissue_RNA_clock"),
                        selectInput("tclock_timepoints", "Select timepoints", choices=timepoints, selected=timepoints, multiple=TRUE),
                        selectInput("tclock_donors", "Select donors", choices=donors, selected=donors, multiple=TRUE),
                        checkboxInput("tclock_split", "Plot difference vs control?", value=TRUE),
                        out = h5("Transcriptional clocks are described here: TO-DO")
          ),
          mainPanel(
            girafeOutput("plot_transcriptional_clock", width = "750px", height = "700px")
          )
        ),
        
        #############################
        ## De novo Gene signatures ##
        #############################
        
        # tabPanel(
        #   title = "De novo gene signatures", id = "gs_denovo",
        #   sidebarPanel2(width=3,
        #                 selectizeInput("gsn_genes", "Add genes", choices = NULL),
        #                 selectizeInput("gsn_x_axis", "Select x-axis", choices = c("Epigenetic age (Horvath2)"="epigenetic_age", "Chronological age"="chronological_age"), selected="chronological_age"),
        #                 selectizeInput("gsn_colourby", "Colour by", choices = c("Epigenetic age (Horvath2)"="epigenetic_age", "Chronological age"="chronological_age", "HGPS"="hgps"), selected="Epigenetic age (Horvath2)"),
        #                 checkboxInput("gsn_include_hgps", "Include HGPS samples?", value = TRUE),
        #                 checkboxInput("gsn_include_longitudinal", "Include longitudinal samples?", value = TRUE),
        #                 checkboxInput("gsn_split_media", "Split by media condition?", value = TRUE),
        #                 checkboxInput("gsn_label_samples", "Label samples?", value = FALSE),
        #                 # selectInput(inputId = "gene_sign_classes", label = "Classes", choices = classes.choices, selected = classes, multiple = TRUE)
        #                 out = h5("Gene signatures are calculated using the UCell package. Scores range from 0 to 1, as a non-rigurous rule of thumb any difference <0.1 is likely insignificant.")
        #   ),
        #   mainPanel(
        #     # girafeOutput("plot_gene_signature_denovo", width = "1000px", height = "400px")
        #     girafeOutput("plot_gene_signature_denovo")
        #   )
        # ),
        
        shiny::tags$style(HTML(".navbar-header { width:100% } .navbar-brand { width: 100%; text-align: center; font-size: 27px }")) 
      )
    )
  )
