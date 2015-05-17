###########################
##-----------------------##
##   F U N C T I O N S   ##
##-----------------------##
###########################

##########################################################
getGeneList <- function(pfamid) {
  #--------------------------------------------------------#
  # Gets standard database names for human genes with the  #
  # selected PFAM ID using a perl script.                  #
  ##########################################################
  program.name <- "enrichPFAM.pl"
  perl.command  <- paste("perl", program.name, pfamid)
  x <- system(perl.command, intern=TRUE)
  genes_list <- strsplit(x, ",")[[1]]
  return (genes_list)
}

##########################################################
humanAnnotation <- function(gene_list, onto) {
  #--------------------------------------------------------#
  # Creates topGOdata object factoring all human genes in  #
  # a gene_list (getGeneList() output) and returns the     #
  # the GOdata object                                      #
  ##########################################################
  org1_term <- annFUN.org(onto, mapping = "org.Hs.eg.db", ID = "symbol")
  org1_allGenes <- unique(unlist(org1_term))
  org1_geneList <- factor(as.integer(org1_allGenes %in% gene_list))
  names(org1_geneList) <- org1_allGenes
  
  org1_GOdata <- new("topGOdata",
                     ontology = onto,
                     allGenes = org1_geneList,
                     annot = annFUN.org,
                     mapping = "org.Hs.eg.db",
                     ID = "symbol")
  return (org1_GOdata)
}

##########################################################
orgAnnotation <- function(org_db, onto) {
  #--------------------------------------------------------#
  # Using an organism annotation database, returns all     #
  # genes in its genome.                                   #
  ##########################################################
  org2_term <- annFUN.org(onto, mapping = org_db, ID = "symbol")
  org2_allGenes  <- unique(unlist(org2_term))
  return (org2_allGenes)
}

##########################################################
classicEnrichTest <- function(topGOdata, gene.universe) {
  #--------------------------------------------------------#
  test.stat <- new("classicCount",
                   testStatistic = GOFisherTest,
                   name = "Fisher test", 
                   allMembers = gene.universe)
  resultFisher <- getSigGroups(topGOdata, test.stat)
  return (resultFisher)
}

##########################################################
weightEnrichTest <- function(topGOdata, gene.universe) {
  #--------------------------------------------------------#
  test.stat <- new("weightCount",
                   testStatistic = GOFisherTest,
                   name = "Fisher test",
                   sigRatio = "ratio")
  
  resultWeight <- getSigGroups(topGOdata, test.stat)
  return (resultWeight)
}

##########################################################
pValue <- function(result) {
  #--------------------------------------------------------#
  pval <- score(result)
  return (pval)
}
