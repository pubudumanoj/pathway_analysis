#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
    
    
    title = "Pathway Analysis",
    
    plotOutput('plot'),
    
    hr(),
    titlePanel( h4("Pathway Enrichment Analysis", align = "center")),
    sidebarLayout(
      sidebarPanel(
        fileInput(inputId = "file", label = "Insert gene list"),
        tags$h4("or", style="text-align: center;"),
        textAreaInput(inputId = "genelist", label = "Paste Gene list here"),
        selectInput("org", "Select organism",
                    c("hsa")
                        ),
        selectInput("genetype", "Input gene ID type", selected = "Symbol",
                    c("Symbol","Entrez",  "hgnc", "ensembl", "fullname", "uniprotswissprot")
        ),
        checkboxGroupInput("dbs", "Choose Databases:",
                           choiceNames =
                             c("Kegg", "Reactome", "GO"),
                           choiceValues =
                             c("Kegg", "Reactome", "GO"),
                           selected = c("Kegg", "Reactome", "GO")
        ),
        
        fluidRow(column(12, align="center", offset = 0, div(actionButton(inputId = "showpathways", "Show Pathways")),
                        
        )),
       
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
        tags$h4("Change options"),
        tabsetPanel(
          tabPanel("Kegg", br(),
                   selectInput("keggfdr", "P Adjusted method", selected = "BH",
                               c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
                   ),
                   textInput("keggqvalue", "q value", value = 0.05),
                   textInput("keggpvalue", "p value", value = 0.05),
                   sliderInput("keggGsize", "Minimum gene count:",
                               min = 1, max = 30,
                               value = 10)
                   ),
          tabPanel("Reactome", br(),
                   selectInput("reactomefdr", "P Adjusted method", selected = "BH",
                               c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
                   ),
                   textInput("reactomeqvalue", "q value", value = 0.05),
                   textInput("reactomepvalue", "p value", value = 0.05),
                   sliderInput("reactomeGsize", "Minimum gene count:",
                               min = 1, max = 30,
                               value = 10)
                   ),
          tabPanel("Go", br(), 
                   selectInput("gofdr", "P Adjusted method", selected = "BH",
                               c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
                   ),
                   textInput("goqvalue", "q value", value = 0.05),
                   textInput("gopvalue", "p value", value = 0.05),
                   sliderInput("goGsize", "Minimum gene count:",
                               min = 1, max = 30,
                               value = 10),
                   tableOutput("table")),
          tabPanel("Final Figure", br(), 
                   selectInput("fdrmethod", "P Adjusted method", selected = "BH",
                               c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
                   ),
                   textInput("fdr", "FDR", value = 0.05),
                   textInput("ftitle", "Figure Title", value = "Pathway Enrichment Analysis"),
                   sliderInput("Gsize", "Minimum gene count:",
                               min = 1, max = 30,
                               value = 10),
                   )
        ),
       
        
      )
    )
    ))
