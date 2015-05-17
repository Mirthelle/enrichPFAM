library("shiny")
source("fun_library.R")
library("topGO")
library("AnnotationDbi")
library("GO.db")
library("org.Hs.eg.db") # Annotated genome for Homo sapiens
library("org.Rn.eg.db") # Annotated genome for Rattus norvegicus
library("org.Ag.eg.db") #Anopheles
library("org.At.tair.db") #Arabidopsis
library("org.Bt.eg.db") # Bovine
library("org.Ce.eg.db") # Worm
library("org.Cf.eg.db") # Canine
library("org.Dm.eg.db") #  Fly
library("org.Dr.eg.db") # Zebrafish
library("org.EcK12.eg.db") # E coli strain k12
library("org.EcSakai.eg.db") # E coli strain Sakai
library("org.Mm.eg.db") # Mouse
library("org.Mmu.eg.db")  # Rhesus
library("org.Pf.plasmo.db") # Malaria
library("org.Pt.eg.db") # Chimp
library("org.Sc.sgd.db") # Yeast
library("org.Sco.eg.db") # Streptomyces coelicolor
library("org.Ss.eg.db") # Pig
library("org.Tgondii.eg.db") # Toxoplasma gondii
library("org.Xl.eg.db") # Xenopus

shinyServer(function(input, output) {
  
  ###############################################
  ##-------------------------------------------##
  ##   R E A C T I V E   C O N D U C T O R S   ##
  ##-------------------------------------------##
  ###############################################
  
  observe({
    if (input$submit == 0)
      return()
    
    isolate({
      ##--------------------##
      ## MOLECULAR FUNCTION ##
      ##--------------------##
    
      GOhumanMF <- reactive ({
        GOhuman <- humanAnnotation(getGeneList(input$pfamid), "MF")
      })
        
      resultFishMF  <- reactive ({
        resultFish <- classicEnrichTest(GOhumanMF(), orgAnnotation(input$org, "MF"))
      })
      
      resultWeightMF  <- reactive ({
        resultWeight <- weightEnrichTest(GOhumanMF(), orgAnnotation(input$org))
      })
      
      AllResMF <- reactive ({
        AllRes <- GenTable(GOhumanMF(), classic = resultFishMF(), weight = resultWeightMF(), 
                           orderBy = "weight", ranksOf = "classic", topNodes = 20)
      })
      
      ##--------------------##
      ## BIOLOGICAL PROCESS ##
      ##--------------------##
      
      GOhumanBP <- reactive ({
        GOhuman <- humanAnnotation(getGeneList(input$pfamid), "BP")
      })
      
      resultFishBP  <- reactive ({
        resultFish <- classicEnrichTest(GOhumanBP(), orgAnnotation(input$org, "BP"))
      })
      
      resultWeightBP  <- reactive ({
        resultWeight <- weightEnrichTest(GOhumanBP(), orgAnnotation(input$org))
      })
      
      AllResBP <- reactive ({
        AllRes <- GenTable(GOhumanBP(), classic = resultFishBP(), weight = resultWeightBP(), 
                           orderBy = "weight", ranksOf = "classic", topNodes = 20)
      })
      
      ##--------------------##
      ## CELLULAR COMPONENT ##
      ##--------------------##
      
      GOhumanCC <- reactive ({
        GOhuman <- humanAnnotation(getGeneList(input$pfamid), "CC")
      })
      
      resultFishCC  <- reactive ({
        resultFish <- classicEnrichTest(GOhumanCC(), orgAnnotation(input$org, "CC"))
      })
      
      resultWeightCC  <- reactive ({
        resultWeight <- weightEnrichTest(GOhumanCC(), orgAnnotation(input$org))
      })
      
      AllResCC <- reactive ({
        AllRes <- GenTable(GOhumanCC(), classic = resultFishCC(), weight = resultWeightCC(), 
                           orderBy = "weight", ranksOf = "classic", topNodes = 20)
      })
      
      #######################
      ##-------------------##
      ##   O U T P U T S   ##
      ##-------------------##
      #######################
      
      ##--------------------##
      ## MOLECULAR FUNCTION ##
      ##--------------------##
      
      output$hist1MF <- renderPlot ({
        hist(pValue(resultFishMF()), 50, xlab = "p-values")
      })
      
      output$hist2MF <- renderPlot ({
        hist(pValue(resultWeightMF()), 50, xlab = "p-values")
      })
      
       output$pvalueMF <- renderDataTable ({
         AllResMF()
       })
      
      output$nodesallMF <- renderPlot ({
         showSigOfNodes(GOhumanMF(), score(resultFishMF()), firstSigNodes = 5, useInfo = "all")
      })
      
      output$nodesdefMF <- renderPlot ({
         showSigOfNodes(GOhumanMF(), score(resultWeightMF()), firstSigNodes = 5, useInfo = "def")
      })
      
      ##--------------------##
      ## BIOLOGICAL PROCESS ##
      ##--------------------##
      
      output$hist1BP <- renderPlot ({
        hist(pValue(resultFishBP()), 50, xlab = "p-values")
      })
      
      output$hist2BP <- renderPlot ({
        hist(pValue(resultWeightBP()), 50, xlab = "p-values")
      })
      
      output$pvalueBP <- renderDataTable ({
        AllResBP()
      })
      
      output$nodesallBP <- renderPlot ({
        showSigOfNodes(GOhumanBP(), score(resultFishBP()), firstSigNodes = 5, useInfo = "all")
      })
      
      output$nodesdefBP <- renderPlot ({
        showSigOfNodes(GOhumanBP(), score(resultWeightBP()), firstSigNodes = 5, useInfo = "def")
      })
      
      ##--------------------##
      ## CELLULAR COMPONENT ##
      ##--------------------##
      
      output$hist1CC <- renderPlot ({
        hist(pValue(resultFishCC()), 50, xlab = "p-values")
      })
      
      output$hist2CC <- renderPlot ({
        hist(pValue(resultWeightCC()), 50, xlab = "p-values")
      })
      
      output$pvalueCC <- renderDataTable ({
        AllResCC()
      })
      
      output$nodesallCC <- renderPlot ({
        showSigOfNodes(GOhumanCC(), score(resultFishCC()), firstSigNodes = 5, useInfo = "all")
      })
      
      output$nodesdefCC <- renderPlot ({
        showSigOfNodes(GOhumanCC(), score(resultWeightCC()), firstSigNodes = 5, useInfo = "def")
      })
    })
  })
  
})