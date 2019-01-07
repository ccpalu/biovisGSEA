###############################
###############################
# NSilico Life Sciences Ltda  #
# author: Cintia Palu         #
#                             #
# GSEA made from interface    #
# input                       #
# 07/01/2019                  #
# version: 0.9                #
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

doGO <- function(name, gnId, uni, db, pvalue, onto){
  
  params <- new("GOHyperGParams",
    geneIds = gnId,
    universeGeneIds = uni,
    annotation = db,
    ontology = onto,
    pvalueCutoff = pvalue,
    conditional = FALSE,
    testDirection = "over")
  
  hgOver = NULL;hgOver <- hyperGTest(params)
  
  return(hgOver)
}#doGo

get.parent <- function(g,family){
  # Function developed to retrieve the parent's of nodes and leaves that are significantly enriched
  p <- NULL
  e <- edges(g)
  
  for(i in 1:length(e)){
    
    if(any(is.element(e[[i]],family))){
      # If the node i is has any family as edge, we include it as
      # family membe. if it wasn't already, it will be added to the
      # end of the string and have its own parents checked from now own
      family <- unique(c(family,names(e)[i]))
    }
  }
  return(family)
}# get.parent

gogetgene <- function(l,gene){
  return(any(gene%in%l))
}

simplicity_gostats <- function(Gene.Ann, SP, Universe.Sets, Universe.Source, Universe.Gene.Ann = NULL, pvalue = 0.025,
                              updateProgress = NULL, output = NULL, data1){
  require(GOstats)
  require(Category)
  require(jsonlite)
  require(colorRamps)
  require(Rgraphviz)
  
  ################################
  # ############################ #
  # #                          # #
  # #       GO  Analysis       # #
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
      
      result = c(result, list(doGO(name = 'Simplicity.GSEA.singleSet',
                            gnId = Gene.Ann$ENTREZID, onto = onto,
                            uni =  universe, db = SP$DB, pvalue = pvalue)))
      names(result)[length(result)] <- paste0(tolower(onto),1)
      }
  }else{
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
        
        result <- c(result, list(doGO(name = paste('Simplicity.GSEA.set',i, sep = "."),
                                gnId = Gene.Ann$ENTREZID[which(Gene.Ann$Set == i)],
                                uni =  universe, db = SP$DB, pvalue = pvalue, onto = onto)))
        names(result)[length(result)] <- paste0(tolower(onto),i)
      }
    }
  }#if/else(is.null(Universe.Sets))
  
  
  #####################
  #       Graph       #
  #####################
  
  temp.graph <- list(nodes = NULL,links = NULL)

  for(I in 1:length(result)){
    ok <- summary(result[[I]], pvalue = 0.1)
    if(dim(ok)[1]>0){
      # selecting the graph element
      g <- result[[I]]@goDag
      #exploring it
      degree(g)
      nodes(g)
      edges(g)

      plot(subGraph(ok[,1],g))

      # selecting the GOId with significative p-value
      target <- ok[,1]

      # selecting the parent nodes for those targets
      g.size <- length(target)

      family <- get.parent(g,target)
      while(length(family)>length(target)){
        target <- family
        family <- get.parent(g,target)
      }
      target <- ok[,1]

      # Creating table with the data for each family 'member'
      allgo <- summary(result[[I]],pvalue = 1.1)
      dim(allgo)
      #73 7
      length(family)
      # 11
      family.table <- allgo[which(allgo[,1]%in%family),]
      dim(family.table)
      #[1] 11 7
      g1 <- subGraph(family,g)

      G <- result[[I]]

      go.genes <- NULL
      #The first version linked the Genes to all family GOs.
      for(i in which(names(g1@nodeData@data)%in%family)){
        # I will link only to the significative GOs
        #for(i in which(names(g1@nodeData@data)%in%target)){
        A <- as.array(g1@nodeData@data[[i]]$geneIds)
        A <- (A[which(g1@nodeData@data[[i]]$geneIds%in%G@geneIds)])
        A <- Gene.Ann$SYMBOL[which(Gene.Ann$ENTREZID%in%A)]
        go.genes <- c(go.genes, list(a = A))
        names(go.genes)[length(go.genes)] <- names(g1@nodeData@data)[i]
      }
      go.genes

      geneIds <- Gene.Ann$SYMBOL[which(Gene.Ann$ENTREZID%in%G@geneIds)]
      gene.edge <- NULL
      for (i in geneIds){
        ggg <- sapply(go.genes,gogetgene, gene = i)
        gene.edge <- c(gene.edge, list(a = list(edges = names(ggg)[which(ggg)])))
        names(gene.edge)[length(gene.edge)] = i
      }

      test <- addNode(node = geneIds, object = g1)

      for(i in geneIds){
        if(length(gene.edge[[i]]$edges)>0){
          test <- addEdge(from = i, to = gene.edge[[i]]$edges, test)
        }
      }
      test.gos <- grep('GO:', test@nodes)

      ################################
      # ############################ #
      # #                          # #
      # #      Creating JSON       # #
      # #                          # #
      # ############################ #
      ################################


      sub.gene <- Gene.Ann[which(Gene.Ann$SYMBOL%in%geneIds),-which(colnames(Gene.Ann) == 'Set')]
      if('logFC'%in%colnames(data1)){
        sub.gene$logFC <- data1[sub.gene[,1],'logFC']
        M <- round(max(abs(sub.gene$logFC)),digits = 1)
        COL <- matlab.like2(M*20+1)
        names(COL) <- ((-M*10):(M*10))/10
        sub.gene$color <- COL[as.character(round(sub.gene$logFC,digits = 1))]
          
      }else{
        sub.gene$logFC <- NA
        sub.gene$color <- 'grey'
        }
      
      if('FDR'%in%colnames(data1)){
        sub.gene$pValue <- data1[sub.gene[,1],'FDR']
      }else{sub.gene$pValue <- NA}


      colnames(sub.gene)[which(colnames(sub.gene) == 'SYMBOL')] <- 'name'
      colnames(sub.gene)[which(colnames(sub.gene) == 'GENENAME')] <- 'longname'
      sub.gene$type <- 'gene'

      temp.node <- as.data.frame(test@nodes[test.gos])
      names(temp.node) <- 'name'
      temp.node$color <- "#FFFFFF"
      temp.node$pValue <- NA
      for(i in family){
        temp.node$pValue[which(temp.node$name == i)] <- family.table$Pvalue[which(family.table[,1] == i)]
      }
      temp.node$type <- colnames(family.table)[1]
      rownames(family.table) <- family.table[,1]
      temp.node$longname <- family.table[as.character(temp.node$name),'Term']
      temp.node$pValue <- family.table[as.character(temp.node$name),'Pvalue']
      nodes <- merge(sub.gene, temp.node, all = TRUE, sort = FALSE)

      rownames(nodes) <- nodes$name

      # edges
      edges <- test@edgeL

      source <- NULL
   
      target <- NULL
      for(i in 1:length(edges)){
        if(length(edges[[i]]$edges)){
          for(j in edges[[i]]$edges){
            source[length(source)+1] <- names(edges)[i]
            target[length(target)+1] <- test@nodes[j]
          }
        }
      }
      temp.edges <- data.frame(source,target)
      temp.edges$value <- 1

      temp.graph$nodes <- rbind(temp.graph$nodes, nodes[setdiff(nodes$name,temp.graph$nodes$name),])
      temp.graph$links <- rbind(temp.graph$links, temp.edges)
    }
  }


  temp.graph$nodes <- temp.graph$nodes[order(temp.graph$nodes$name),]
  temp.graph$nodes <- temp.graph$nodes[order(temp.graph$nodes$type),]
  temp.graph$nodes$id <- 0:(length(temp.graph$nodes$name)-1)
  levels(temp.graph$links$source) <- temp.graph$nodes[levels(temp.graph$links$source),'id']
  levels(temp.graph$links$target) <- temp.graph$nodes[levels(temp.graph$links$target),'id']

  row.names(temp.graph$nodes) <- NULL

    #####################
  
  return(list(gsea=result, graph = temp.graph))#toJSON(temp.graph)))
  
}#simplicity_gostats
