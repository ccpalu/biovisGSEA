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
    
    return(summary(hgOver,"htmlLinks" = TRUE))
  }
  
}

simplicity_gostats <- function(Gene.Ann, SP, Universe.Sets, Universe.Source, Universe.Gene.Ann = NULL, pvalue = 0.025,
                              updateProgress = NULL, output = NULL){
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
    for (onto in c('BP','MF','CC')){
      if (is.function(updateProgress)) {
        text = switch (onto,
          'BP' = 'Step 1 of 3: GO Biological Processes',
          'MF' = 'Step 2 of 3: GO Molecular Functions',
          'CC' = 'Step 3 of 3: GO Cellular Components'
        )
        updateProgress(detail = text)
      }#if (is.function(updateProgress)) 
      
      result = c(result, list(as.data.frame(doGO(name = 'Simplicity.GSEA.singleSet',
                            gnId = Gene.Ann$ENTREZID, onto = onto,
                            uni =  universe, db = SP$DB, pvalue = pvalue))))
      names(result)[length(result)] <- paste0(tolower(onto),1)
      }
  }else{#
    steps <- length(Universe.Sets)*3
    countSteps <- 1
    for(i in Universe.Sets){
      for(onto in c('BP','MF','CC')){
        if (is.function(updateProgress)) {
          text = switch (onto,
            'BP' = paste('Step', countSteps, 'of', steps,': GO Biological Processes for gene set', i),
            'MF' = paste('Step', countSteps, 'of', steps,': GO Molecular Functions for gene set', i),
            'CC' = paste('Step', countSteps, 'of', steps,': GO Cellular Components for gene set', i)
          )
          updateProgress(value = countSteps / steps, detail = text)
          countSteps <- countSteps + 1
        }#if (is.function(updateProgress)) 
        
        result <- c(result, list(as.data.frame(doGO(name = paste('Simplicity.GSEA.set',i, sep = "."),
                                gnId = Gene.Ann$ENTREZID[which(Gene.Ann$Set == i)],
                                uni =  universe, db = SP$DB, pvalue = pvalue, onto = onto))))
        names(result)[length(result)] <- paste0(tolower(onto),i)
      }
    }
  }#if/else(is.null(Universe.Sets))

  return(result)
}
