###############################
###############################
# NSilico Life Sciences Ltda  #
# author: Cintia Palu         #
#                             #
# GSEA made from interface    #
# input                       #
# 19/11/2018                  #
# version: 0.2                #
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

doGO <- function(name, gnId, uni, db, pvalue, 
  onto = c('BP','MF','CC')){
  
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
    #save(hgOver,file = paste(name, i, "hgOver.RData", sep = "."))
    # htmlReport(hgOver, file = paste(name, i, "hgOver.html", sep = "."), summary.args = list("htmlLinks" = TRUE))
    return(summary(hgOver,"htmlLinks" = TRUE))
    # testDirection(params) <- 'under'
    # hgUnder <- hyperGTest(params)
    # 
    # htmlReport(hgUnder, file = paste(name, i, "hgUnder.html", sep = "."), summary.args = list("htmlLinks" = TRUE))
  }
  # params <- new("KEGGHyperGParams", geneIds = gnId,
  #   annotation = db, universeGeneIds = uni,
  #   pvalueCutoff = pvalue, testDirection = "over")
  # hgOver <- hyperGTest(params)
  # htmlReport(hgOver, file = paste(name, "KEGG.hgOver.html", sep = "."), summary.args = list("htmlLinks" = TRUE))
  # 
  # testDirection(params) <- 'under'
  # hgUnder <- hyperGTest(params)
  # 
  # htmlReport(hgUnder, file = paste(name, "KEGG.hgUnder.html", sep = "."), summary.args = list("htmlLinks" = TRUE))
  
}

simplicity_gostats <- function(Gene.Ann, SP, Universe.Sets, Universe.Source, Universe.Gene.Ann = NULL, pvalue = 0.025){
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
  result <- list()
  if(is.null(Universe.Sets)){
    for (onto in c('BP','MF','CC')){#!!!
      result = c(result, list(as.data.frame(doGO(name = 'Simplicity.GSEA.singleSet',
                            gnId = Gene.Ann$ENTREZID, onto = onto,
                            uni =  universe, db = SP$DB, pvalue = pvalue))))
      names(result)[length(result)] <- paste0(tolower(onto),1)
      }
  }else{#
    for(i in Universe.Sets){
      for(onto in c('BP','MF','CC')){
        result <- c(result, list(as.data.frame(doGO(name = paste('Simplicity.GSEA.set',i, sep = "."),
                                gnId = Gene.Ann$ENTREZID[which(Gene.Ann$Set == i)],
                                uni =  universe, db = SP$DB, pvalue = pvalue, onto = onto))))
        names(result)[length(result)] <- paste0(tolower(onto),i)
      }
    }
  }
  
  return(result)
}
