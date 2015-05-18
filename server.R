library("shiny")
source("fun_library.R")
source("packages_load.R")

shinyServer(function(input, output) {
  
  observe({
    if (input$submit == 0)
      return()
    
    isolate({
  
      ###############################################
      ##-------------------------------------------##
      ##   R E A C T I V E   C O N D U C T O R S   ##
      ##-------------------------------------------##
      ###############################################
  

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
        hist(pValue(resultFishMF()), 50, main = "p-values for Classic Fisher enrichment test", 
             xlab = "p-values", col = "skyblue", border = "white")
      })
      
      output$hist2MF <- renderPlot ({
        hist(pValue(resultWeightMF()), 50, main = "p-values for Weight Fisher enrichment test", 
             xlab = "p-values", col = "skyblue", border = "white")
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
      
#       output$downloadData <- downloadHandler(
#         printGraph(GOdataMF(), resultWeightMF(), firstSigNodes = 10, resultFisMF(), fn.prefix = "tGO", useInfo = "def")
#       )
      
      ##--------------------##
      ## BIOLOGICAL PROCESS ##
      ##--------------------##
      
      output$hist1BP <- renderPlot ({
        hist(pValue(resultFishBP()), 50, main = "p-values for Classic Fisher enrichment test", 
             xlab = "p-values", col = "skyblue", border = "white")
      })
      
      output$hist2BP <- renderPlot ({
        hist(pValue(resultWeightBP()), 50, main = "p-values for Weight Fisher enrichment test", 
             xlab = "p-values", col = "skyblue", border = "white")
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
        hist(pValue(resultFishCC()), 50, main = "p-values for Classic Fisher enrichment test", 
             xlab = "p-values", col = "skyblue", border = "white")
      })
      
      output$hist2CC <- renderPlot ({
        hist(pValue(resultWeightCC()), 50, main = "p-values for Weight Fisher enrichment test", 
             xlab = "p-values", col = "skyblue", border = "white")
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