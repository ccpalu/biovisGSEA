###############################
###############################
# NSilico Life Sciences Ltda  #
# author: Cintia Palu         #
#                             #
# GSEA made from interface    #
# input                       #
# 05/05/2017                  #
# version: 0.1                #
###############################
###############################

# Input: genes list -> might be an R variable resulting from a 
#        differential expression analysis (from limma as example)
#        or the user will be able to input a simple gene list
# Gene IDs: the script recognises Entrez IDs, ...
# The script used Hypergeometric test to identify GO ontologies
# tht migth match the gene set role


######################################
######## Required Libraries###########
######################################


simplicity_gostats <- function(Gene.Ann, SP, Universe.Sets,Universe.Source,pvalue = 0.025){
  library(GOstats)
  #library(GSEABase)
  library(Category)
  ##############################
  ##########  Input  ###########
  ##############################
  
  ################################
  # ############################ #
  # #                          # #
  # #       GO  Analysis       # #
  # #     Custom Gene List     # #
  # #                          # #
  # ############################ #
  ################################
  
  
  #####################
  #   Gene Universe   #
  #####################
  
  universe <- switch(Universe.Source,
                   'current' = Gene.Ann$ENTREZID,
                   'new.source' = Universe.Gene.Ann$ENTREZID,
                   #NCBI full collection is default -> it takes a long time
                   mappedkeys(get(gsub('.db','ACCNUM',SP$DB)) )
                   
  )
  # Removing any duplicates
  if(any(duplicated(universe))){
    universe <- universe[-which(duplicated(universe))]
  }
  #####################
  #      GO GSEA      #
  #####################
  # Removing any duplicates
  if(any(duplicated(Gene.Ann$ENTREZID))){
    Gene.Ann <- Gene.Ann[-which(duplicated(Gene.Ann$ENTREZID)),]
  }
  doGO <- function(name,gnId, uni,db,pvalue){
    onto <- c('BP','MF','CC')
    for(i in onto){
      params <- new("GOHyperGParams",
                    geneIds = gnId,
                    universeGeneIds = uni,
                    annotation = db,
                    ontology = i,
                    pvalueCutoff = pvalue,
                    conditional = FALSE,
                    testDirection = "over")
      hgOver <- hyperGTest(params)
      htmlReport(hgOver, file = paste(name, i, "hgOver.html", sep = "."), summary.args = list("htmlLinks" = TRUE))
      
      testDirection(params) <- 'under'
      hgUnder <- hyperGTest(params)
      
      htmlReport(hgUnder, file = paste(name, i, "hgUnder.html", sep = "."), summary.args = list("htmlLinks" = TRUE))
    }
    params <- new("KEGGHyperGParams", geneIds = gnId,
                 annotation = db, universeGeneIds = uni,
                 pvalueCutoff = pvalue, testDirection = "over")
    hgOver <- hyperGTest(params)
    htmlReport(hgOver, file = paste(name, "KEGG.hgOver.html", sep = "."), summary.args = list("htmlLinks" = TRUE))
    
    testDirection(params) <- 'under'
    hgUnder <- hyperGTest(params)
    
    htmlReport(hgUnder, file = paste(name, "KEGG.hgUnder.html", sep = "."), summary.args = list("htmlLinks" = TRUE))
  
  }
  
  
  for(i in Universe.Sets){
  
    doGO(paste('Simplicity.GSEA.set',i, sep = "."),
         Gene.Ann$ENTREZID[which(Gene.Ann$Set == i)], universe,
         SP$DB,pvalue)
  }
  
}
