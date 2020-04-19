library(shiny)


# LOAD GENE NAMES
geneIDs=as.character(read.csv("Inputs/ReducedGenes.csv")$x)
  
  
# Define UI for application that draws a histogram
shinyUI(fluidPage(
  navbarPage("BRCA Methyl Shiny",
             # META DATA PLOT TAB
             tabPanel("Probe Annotation Plots",
                      fluidRow(sidebarLayout(
                        sidebarPanel(width=2,
                          radioButtons("metaDataPlotType", "Plot",c("Probe Mapping Categories"="1", "Probes Mapping to Unique Genes"="2"))
                        ),
                        mainPanel(width = 9,align="center",
                          plotOutput("metaDataPlot")
                        )
                      ))
             ),
             # GENE DATA PLOT TAB
             tabPanel("Gene Methylation Plots",
                      fluidRow(sidebarLayout(
                        sidebarPanel(width=2,
                          selectInput("geneID", "Select Gene",selected = "CDH1",
                                       c(geneIDs)),
                          selectInput("comparyBy", "Select Clinical Variable",selected = "ER Status",
                                         c("ER Status","Histology","PAM50")),
                          radioButtons("genePlotType", "Plot",c("Heatmap"="1", "TSS Scatter Plots"="2", "Gene Body Scatter Plots"="3")),
                            checkboxInput("DoColSplit", "Split By Selected Clinical Variable",value = F)
                          ),
                        mainPanel(width = 9,align="center",
                          plotOutput("genePlot")
                        )
                      ))
             ),
             # DIFFERENTIALLY METHYLATED REGION TAB
             tabPanel("About",
                     fluidRow(column(4,
                                     includeMarkdown("https://raw.githubusercontent.com/osamashiraz/BRCA_MethylShiny/master/README.md")
                     ),
                       column(5,
                              img(class="img-polaroid",
                                  src=paste0("https://www.cancer.gov/sites/g/files/xnrzdm211/files/styles/cgov_featured/public/cgov_image/media_image/100/900/5/files/TCGA%20people%20and%20layers%20of%20data%20425x319.jpg?h=982f41e1&itok=ZGxIcnmW")),
                              tags$small(
                                "Source: Cancer Gov ",
                                a(href="https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga",
                                  "Cancer.gov")
                              )
                       )
                     )
                        
             )
  )
))
