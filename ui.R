library(shiny)

shinyUI(fluidPage(
  # Title
  titlePanel("enrichPFAM"),
  
  #Sidebar with a slider input for the number of bins
  sidebarLayout(
    sidebarPanel(
      
      # PFAM ID box
      textInput("pfamid", label = h4("PFAM ID"),
                value = "Enter PFAM ID"),
      helpText("E.g. PF00230"),
      
      # Organism for doing enrichment analysisis
      selectInput("org", label = h4("Select organism"), 
                  choices = c("Anopheles" = "org.Ag.eg.db", "Arabidopsis" = "org.At.tair.db", 
                              "Bovine" = "org.Bt.eg.db", "Canine" = "org.Cf.eg.db", 
                              "Chimp" = "org.Pt.eg.db", "E. coli strain k12" = "org.EcK12.eg.db", 
                              "Fly" = "org.Dm.eg.db", 
                              "Malaria" = "org.Pf.plasmo.db", "Mouse" = "org.Mm.eg.db", 
                              "Pig" = "org.Ss.eg.db","Rat" = "org.Rn.eg.db",
                              "Rhesus" = "org.Mmu.eg.db", "Streptomyces coelicolor" = "org.Sco.eg.db",
                              "Toxoplasma gondii" = "org.Tgondii.eg.db", "Xenopus" = "org.Xl.eg.db",
                              "Yeast" = "org.Sc.sgd.db", "Zebrafish" = "org.Dr.eg.db"), 
                  selected = "org.Rn.eg.db"),
      
      # Action button for submiting request.
      actionButton("submit", "Submit")
    ),
    
    
    # RESULTS PANEL
    mainPanel(
      tabsetPanel(
        
        # Help and Usage
        tabPanel("Help",
                 h1("Welcome to enrichPFAM"),
                 p("enrichPFAM is a web tool for performing Gene Ontology (GO) enrichment analysis.
                   It uses a PFAM ID as input, selects all human genes in the family and performs an enrichment test
                   against a selected organism."), 
                 p("The results are provided in three different ways: (1) Histograms for p-values,
                   (2) a summary table for GO terms and their scores, (3) GO graphs with significative nodes and their parents."),
                 h3("Usage"),
                 p("This are the main steps you need to do in order to run this app:"),
                 tags$ol(
                   tags$li("Write your", tags$b("PFAM ID"), "of interest in the text box."), 
                   tags$li("Select the", tags$b("organism"), "in the list with which you want to compare human genes in the PFAM family."), 
                   tags$li("Click", tags$b("'Sumbit'"),  "for sending your request."),
                   tags$li("Select the", tags$b("ontology"), "tab in which you are interested.")
                   ),
                 h6(tags$b("Warning"), ": Program will not start calculations until you select an ontology tab. 
                    Each time you select a tab, calculations will restart for the selected ontology 
                    with the same data you submitted unless you submit other data. 
                    Be patient, it may take a while."),
                 h3("Download"),
                 p("If you are interested in the code for this app, you can download it freely", a("here", 
                                                                                                   href="https://github.com/Mirthelle/enrichPFAM",
                                                                                                   target="_blank"),
                   ".")
                 ),
        
        # Ontology: Molecular Munction
        tabPanel("Molecular Function",
                 tabsetPanel(
                    tabPanel("Histogram",
                              plotOutput("hist1MF", height = 600),
                              plotOutput("hist2MF", height = 600)
                    ),
                    tabPanel("Data table",
                             h2("Results Summary"),
                            dataTableOutput("pvalueMF"),
                            downloadLink("downloadTableMF", "Download Table")
                    ),
                    tabPanel("GO Graphs",
                             h3("GO graph with classic enrichment analysis"),
                            plotOutput("nodesallMF", height=1000),
#                             downloadLink('downloadData1', 'Download Graph'),
                            h3("GO graph with weight enrichment analysis"), 
                            plotOutput("nodesdefMF", height=1000)
#                             downloadLink('downloadData2', 'Download Graph')
                  )
                 )
        ),
        
        # Ontology: Biological Process
        tabPanel("Biological Process",
                 tabsetPanel(
                   tabPanel("Histogram",
                            plotOutput("hist1BP", height = 600),
                            plotOutput("hist2BP", height = 600)
                            ),
                   tabPanel("Data table",
                            h2("Results Summary"),
                            dataTableOutput("pvalueBP"),
                            downloadLink("downloadTableBP", "Download Table")
                            ),
                   tabPanel("GO Graphs",
                            h3("GO graph with classic enrichment analysis"),
                            plotOutput("nodesallBP", height=1000),
#                             downloadLink('downloadData3', 'Download Graph'),
                            h3("GO graph with weight enrichment analysis"),
                            plotOutput("nodesdefBP", height=1000)
#                             downloadLink('downloadData4', 'Download Graph')
                            )
                   )
                ),
        
        # Ontology: Celular Component
        tabPanel("Cellular Component",
                 tabsetPanel(
                   tabPanel("Histogram",
                            plotOutput("hist1CC"),
                            plotOutput("hist2CC")
                            ),
                   tabPanel("Data table",
                            h2("Results Summary"),
                            dataTableOutput("pvalueCC"),
                            downloadLink("downloadTableCC", "Download Table")
                            ),
                   tabPanel("GO Graphs",
                            h3("GO graph with classic enrichment analysis"),
                            plotOutput("nodesallCC", height=1000),
#                             downloadLink('downloadData5', 'Download Graph'),
                            h3("GO graph with weight enrichment analysis"),
                            plotOutput("nodesdefCC", height=1000)
#                             downloadLink('downloadData6', 'Download Graph')
                            )
                   )
        )
      )
    )
  )
))