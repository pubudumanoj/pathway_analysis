#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(data.table)
library(MAGeCKFlute)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(DT)
library(shinycssloaders)


# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {

  find_pathways <- function(){
    filelistgene <- input$file$datapath
    if(is.null(filelistgene) & input$genelist == ""){
      showModal(modalDialog(
        title = "Warning",
        "Please upload or paste a gene list",
        easyClose = TRUE
      ))
    
    } else{
      
    if(input$genelist != ""){
      gene_list <- strsplit(input$genelist, "\n") %>% unlist()
      gene_list <- data.frame(V1=gene_list)
    } else{
      gene_list <- fread(filelistgene, header = F)
      
    }
    
    
    fdr =input$fdr
    
    i=1
    
    gene_list <- unique(gene_list)
    
    if(input$genetype!="Entrez"){
      Genes_GeneHancer <- lapply(X = gene_list, FUN = function(x) {x <- TransGeneID(x, fromType =input$genetype, toType = "Entrez", organism = input$org, ensemblHost = "www.ensembl.org")})
      rm(gene_list)
    } 
    
    
    
    ##############################################################################################################
    ###########                               Kegg pathway analysis                                    ###########
    ##############################################################################################################
    
    if("Kegg" %in% input$dbs){
      KEGG_GeneHancer <- lapply(X = Genes_GeneHancer, FUN = function(x) {x <- enrichKEGG(x, organism=input$org, pAdjustMethod=input$keggfdr, qvalueCutoff=as.numeric(input$keggqvalue), pvalueCutoff=as.numeric(input$keggpvalue), minGSSize = input$keggGsize)})
      KEGG_GeneHancer <- lapply(X = KEGG_GeneHancer, FUN = function(x) {x <- setReadable(x = x, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")})
      Manhattan_KEGG <- KEGG_GeneHancer[i][[1]]@result
      rm(KEGG_GeneHancer)
      Manhattan_KEGG$Approach <- "KEGG"
    }
   
    
    
    ##############################################################################################################
    ###########                                  Reactome  analyses                                    ###########
    ##############################################################################################################
    
    if("Reactome" %in% input$dbs){
      org="hsa"
      if(input$org=="hsa"){
        org="human"
      }
      Reactome_GeneHancer <- lapply(X = Genes_GeneHancer, FUN = function(x) {x <- enrichPathway(x, organism=org, pAdjustMethod=input$reactomefdr, qvalueCutoff=as.numeric(input$reactomeqvalue), pvalueCutoff=as.numeric(input$reactomepvalue), minGSSize = input$reactomeGsize, readable = TRUE)})
      Manhattan_Reactome <- Reactome_GeneHancer[i][[1]]@result
      rm(Reactome_GeneHancer)
      Manhattan_Reactome$Approach <- "Reactome"
    }
    ##############################################################################################################
    ###########                                   GO BP analyses                                       ###########
    ##############################################################################################################
    
    if("GO" %in% input$dbs){
      GO_BP_GeneHancer <- lapply(X = Genes_GeneHancer, FUN = function(x) {x <- enrichGO(x, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod=input$gofdr, qvalueCutoff = as.numeric(input$goqvalue), pvalueCutoff = as.numeric(input$gopvalue), minGSSize = input$goGsize, readable = TRUE)})
      Manhattan_GO <- GO_BP_GeneHancer[i][[1]]@result
      rm(GO_BP_GeneHancer)
      Manhattan_GO$Approach <- "GO"
      }
    
    #############################################################################################
    #################             Manhattan plot for all tested databases      ##################
    #############################################################################################
    
     
      Manhattan_data <- ""
     
      manhat <- list()

    
    if("Kegg" %in% input$dbs){
      manhat <- c(manhat, list(Manhattan_KEGG))
      rm(Manhattan_KEGG)
    } 
    if ("Reactome" %in% input$dbs){
      manhat <- c(manhat, list(Manhattan_Reactome))
      rm(Manhattan_Reactome)
    }
    if ("GO" %in% input$dbs){
      manhat <- c(manhat, list(Manhattan_GO))
      rm(Manhattan_GO)
    }
      Manhattan_data <- bind_rows(manhat)
      rm(manhat)
      
      
      Manhattan_data_table <- Manhattan_data
      rownames(Manhattan_data_table) <- NULL
      Manhattan_data_table <- Manhattan_data_table %>% mutate(pvalue = round(pvalue, 3), p.adjust = round(p.adjust, 3), qvalue = round(qvalue, 3))
      Manhattan_data_table <- select(Manhattan_data_table, -geneID)
      Manhattan_data <- Manhattan_data[order(Manhattan_data$ID),]
      Manhattan_data <- subset(Manhattan_data, Count >= 5)
      Manhattan_data$BH_combined = p.adjust(Manhattan_data$pvalue, method = input$fdrmethod)
      Manhattan_data <- Manhattan_data[order(Manhattan_data$BH_combined),]
      
      # Manhattan plot 
      #check whether there are more significant pathways. if there are no significant pathways. there will be no lables
      if(dim(Manhattan_data)[1]>0){
        Selected_terms <- c(if (dim(subset(Manhattan_data, BH_combined < fdr & Approach == "GO"))[1] > input$Gsize) subset(Manhattan_data, Approach == "GO")[1:10,"ID"] else subset(Manhattan_data, Approach == "GO")[1:dim(subset(Manhattan_data, BH_combined < fdr & Approach == "GO"))[1],"ID"], 
                            if (dim(subset(Manhattan_data, BH_combined < fdr & Approach == "KEGG"))[1] > input$Gsize) subset(Manhattan_data, Approach == "KEGG")[1:5,"ID"] else subset(Manhattan_data, Approach == "KEGG")[1:dim(subset(Manhattan_data, BH_combined < fdr & Approach == "KEGG"))[1],"ID"],
                            if (dim(subset(Manhattan_data, BH_combined < fdr & Approach == "Reactome"))[1] > input$Gsize) subset(Manhattan_data, Approach == "Reactome")[1:5,"ID"] else subset(Manhattan_data, Approach == "Reactome")[1:dim(subset(Manhattan_data, BH_combined < fdr & Approach == "Reactome"))[1],"ID"]
        )
        
        
        
        p <- ggplot(Manhattan_data, aes(x = ID, y = -log10(BH_combined), color = Approach, label = Description)) +
          geom_point(aes(size = Count), alpha = 0.5) +
          ylab("-log10 (FDR)") +
          scale_size(range = c(1,10)) +
          geom_text_repel(data = subset(Manhattan_data, ID %in% Selected_terms), aes(label = Description),
                          segment.color = "grey75", color = "grey25", nudge_y = 2.5, angle = 0, size = 4,
                          show.legend = FALSE, force = 10) +
          geom_hline(yintercept = 1.3, color="black", linetype="dotted", size = 1) +
          ylim(0,5) +
          theme_bw() +
          ggtitle(input$ftitle)+
          theme(axis.title = element_text(size = 14, face = "bold"), axis.text.y = element_text(size = 12),
                axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank())
        
      }
      else{
        print("No enoguh data")
      }

    
    return(list(name=p, table=Manhattan_data_table))
  }
  }
  
  
  observeEvent(input$showpathways, {

    
    output$plot <- renderPlot({
      finaloutput <- find_pathways()
      output$table <- DT::renderDataTable({
        
        DT::datatable(finaloutput$table, options = list(lengthMenu = c(5, 30, 50), pageLength = 5))
      })
      finaloutput$name
      
    })
     
    
    output$spinner <- renderUI({
      withSpinner(plotOutput("plot"))
    })
    
    
  })
  observeEvent(input$testgenes,{
    
    testegenes <- paste(readLines("data/testgenes.txt"), collapse = "\n")
   
    updateTextAreaInput(session, "genelist", value = testegenes)
    
  })
  
  observeEvent(input$cleargenes,{
    
    
    
    updateTextAreaInput(session, "genelist", value = "")
    
  })


  observe({
    if(length(input$dbs) < 1){
      updateCheckboxGroupInput(session, "dbs", selected= "Kegg")
    }
    
  })

})
