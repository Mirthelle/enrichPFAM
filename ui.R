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
      helpText("E.g. PF00870 or PF00230"),
      
      # Organism for doing enrichment analysisis
      selectInput("org", label = h4("Select organism"), 
                  choices = c("Anopheles" = "org.Ag.eg.db", "Arabidopsis" = "org.At.tair.db", 
                              "Bovine" = "org.Bt.eg.db", "Canine" = "org.Cf.eg.db", 
                              "Chimp" = "org.Pt.eg.db", "E. coli strain k12" = "org.EcK12.eg.db", 
                              "E. coli strain Sakai" = "org.EcSakai.eg.db", "Fly" = "org.Dm.eg.db", 
                              "Malaria" = "org.Pf.plasmo.db", "Mouse" = "org.Mm.eg.db", 
                              "Pig" = "org.Ss.eg.db","Rat" = "org.Rn.eg.db",
                              "Rhesus" = "org.Mmu.eg.db", "Streptomyces coelicolor" = "org.Sco.eg.db",
                              "Toxoplasma gondii" = "org.Tgondii.eg.db", "Xenopus" = "org.Xl.eg.db",
                              "Yeast" = "org.Sc.sgd.db", "Zebrafish" = "org.Dr.eg.db"), 
                  selected = "org.Rn.eg.db"),
#       
#       # Ontologies in which doing the analysis:
#       checkboxGroupInput("ontology", 
#                          label = h4("Ontologies"), 
#                          choices = c("Molecular Function" = "MF", 
#                                       "Cellular Component" = "CC", 
#                                       "Biological Process" = "BP"),
#                                 selected = "MF"),
      actionButton("submit", "Submit")
    ),
    
    
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        
        tabPanel("Molecular Function",
                 tabsetPanel(
                    tabPanel("Histogram",
                              plotOutput("hist1MF", width = 400, height = 300),
                              plotOutput("hist2MF", width = 400, height = 300)
                    ),
                    tabPanel("Data table", 
                            dataTableOutput("pvalueMF")
                    ),
                    tabPanel("GO Graphs",
                            plotOutput("nodesallMF", width = 1000, height = 1500),
                            plotOutput("nodesdefMF", width = 1000, height = 1500)
                  )
                 )
        ),
        tabPanel("Biological Process",
                 tabsetPanel(
                   tabPanel("Histogram",
                            plotOutput("hist1BP", width = 400, height = 300),
                            plotOutput("hist2BP", width = 400, height = 300)
                            ),
                   tabPanel("Data table",
                            dataTableOutput("pvalueBP")
                            ),
                   tabPanel("GO Graphs",
                            plotOutput("nodesallBP"),
                            plotOutput("nodesdefBP")
                            )
                   )
                ),
        tabPanel("Cellular Component",
                 tabsetPanel(
                   tabPanel("Histogram",
                            plotOutput("hist1CC", width = 400, height = 300),
                            plotOutput("hist2CC", width = 400, height = 300)
                            ),
                   tabPanel("Data table",
                            dataTableOutput("pvalueCC")
                            ),
                   tabPanel("GO Graphs",
                            plotOutput("nodesallCC", width = 800, height = 1000),
                            plotOutput("nodesdefCC", width = 800, height = 1000)
                            )
                   )
        )
      )
    )
  )
))