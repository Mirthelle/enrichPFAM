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
                           orderBy = "weight", ranksOf = "classic", topNodes = 50)
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
         showSigOfNodes(GOhumanMF(), score(resultFishMF()), firstSigNodes = 15, useInfo = "all", swPlot = TRUE)
      })
      
      output$nodesdefMF <- renderPlot ({
         showSigOfNodes(GOhumanMF(), score(resultWeightMF()), firstSigNodes = 15, useInfo = "def")
      })
      
      output$downloadTableMF <- downloadHandler(
        filename = "dataTable_MF.csv",
        content <- function(file) {
          write.csv(AllResMF(), file)
        }
      )
      
#        output$downloadData1 <- downloadHandler(
#          filename = "nodes_all",
#          content  <- function (file) {
#            file.copy(printNodes(GOhumanMF(), resultWeightMF(), "tGO_MF", "all"), file)
#          }
#        )
#       
#       output$downloadData2 <- downloadHandler(
#         filename = "nodes_def",
#         content  <- function (file) {
#           printNodes(GOhumanMF(), resultWeightMF(), "def")
#         }
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
        showSigOfNodes(GOhumanBP(), score(resultFishBP()), firstSigNodes = 15, useInfo = "all")
      })
      
      output$nodesdefBP <- renderPlot ({
        showSigOfNodes(GOhumanBP(), score(resultWeightBP()), firstSigNodes = 15, useInfo = "def")
      })
      
      output$downloadTableBP <- downloadHandler(
        filename = "dataTable_BP.csv",
        content <- function(file) {
          write.csv(AllResBP(), file)
        }
      )
    
#       output$downloadData3 <- downloadHandler(
#         filename = "nodes_all",
#         content  <- function (file) {
#           printNodes(GOhumanBP(), resultFishBP(), "tGO_BP", "all")
#         }
#       )
#       
#       output$downloadData4 <- downloadHandler(
#         filename = "nodes_def",
#         content  <- function (file) {
#           printNodes(GOhumanBP(), resultWeightBP(), "tGO_BP", "def")
#         }
#       )
      
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
        showSigOfNodes(GOhumanCC(), score(resultFishCC()), firstSigNodes = 15, useInfo = "all")
      })
      
      output$nodesdefCC <- renderPlot ({
        showSigOfNodes(GOhumanCC(), score(resultWeightCC()), firstSigNodes = 15, useInfo = "def")
      })

      output$downloadTableCC <- downloadHandler(
        filename = "dataTable_CC.csv",
        content <- function(file) {
          write.csv(AllResCC(), file)
        }
      )
      
#       output$downloadData5 <- downloadHandler(
#         filename = "nodes_all",
#         content  <- function (file) {
#           printNodes(GOhumanMF(), resultFishMF(), "tGO_CC", "all")
#         }
#       )
#       
#       output$downloadData6 <- downloadHandler(
#         filename = "nodes_def",
#         content  <- function (file) {
#           printNodes(GOhumanCC(), resultWeightCC(), "tGO_CC", "def")
#         }
#       )
    })
  })
  
})