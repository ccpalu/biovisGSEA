###############################
###############################
# NSilico Life Sciences Ltda  #
# author: Cintia Palu         #
#                             #
# GSEA interactive Go graphs  #
#                             #
# 26/02/2018                  #
# version: 0.8                #
###############################
###############################

# Input: genes list -> might be an R variable resulting from a 
#        differential expression analysis (from limma as example)
#        or the user will be able to input a simple gene list
# Gene IDs: the script recognises Entrez IDs, ...
# The script used Hypergeometric test to identify GO ontologies
# tht migth match the gene set role


### Seting the folder
local(setwd("C:\\Users\\cinti\\OneDrive\\Documents\\BioVis\\GSEA2018"))
folder=format(Sys.time(), "%Y%m%d_%H%M%S")
dir.create(folder)
local(setwd(paste0('.\\',folder)))
#source('GSEA\\gsea.interface.0.8.2.R')
######################################
######## Required Libraries###########
######################################

## Installing ##

# source("http://bioconductor.org/biocLite.R")
# 
# biocLite("GOstats")
# biocLite("Category")
# biocLite("RCyjs")
library(jsonlite)
library(colorRamps)
library(GOstats)
library(Category)
library(Rgraphviz)
#library(RCyjs)
##############################
##########  Input  ###########
##############################
# breast.hd=read.table(file = "../DEA.cancerBrCa.versus.cancerHD.txt",
#                      header = TRUE, sep = '\t', row.names = 1)
# 
# rownames(breast.hd)=gsub('[.][0-9]*','',rownames(breast.hd))
# write.table(x = breast.hd, file = "../DEA.cancerBrCa.versus.cancerHD.txt",col.names = NA,
#             sep = '\t', row.names=TRUE)
# breast.hd$set=NA
# 
# breast.hd$set[which(breast.hd$FDR<=0.05)]='sig'
# summary(as.factor(breast.hd$set))
# # sig NA's 
# # 6591 8832 
# write.table(x = breast.hd, file = "../cancerBrCa.versus.cancerHD.txt",col.names = NA,
#             sep = '\t', row.names=TRUE)
breast.hd=read.table(file = "../cancerBrCa.versus.cancerHD.txt",
                     header = TRUE, sep = '\t', row.names = 1)


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

universe= switch(Universe.Source,
                 'current' = Gene.Ann$ENTREZID,
                 'new.source' = Universe.Gene.Ann$ENTREZID,
                 #NCBI full collection is default -> it takes a long time
                 mappedkeys(get(gsub('.db','ACCNUM',SP$DB)) )
                 
)
# Removing any duplicates
if(any(duplicated(universe))){
  universe = universe[-which(duplicated(universe))]
}

#####################
#      GO GSEA      #
#####################
# Removing any duplicates
if(any(duplicated(Gene.Ann$ENTREZID))){
  Gene.Ann = Gene.Ann[-which(duplicated(Gene.Ann$ENTREZID)),]
}
doGO = function(name,gnId, uni,db,pvalue){
  onto=c('BP','MF','CC')
  r=NULL
  for(i in onto){
    params <- new("GOHyperGParams",
                  geneIds=gnId,
                  universeGeneIds=uni,
                  annotation=db,
                  ontology=i,
                  pvalueCutoff=pvalue,
                  conditional=FALSE,
                  testDirection="over")
    hgOver=NULL;hgOver <- hyperGTest(params)
    if(!is.null(hgOver)){
      htmlReport(hgOver, file=paste(name, i, "hgOver.html", sep="."), summary.args=list("htmlLinks"=TRUE))
    }
    r=c(r,Over = hgOver)
    names(r)[length(r)]=paste(i,'over',sep='_')
  }
  return(r)
}

if(is.null(Universe.Sets)){
  res=doGO(name='Simplicity.GSEA',
           gnId=Gene.Ann$ENTREZID, uni=universe,
           db=SP$DB,pvalue=1)
  save(res,file=paste('gsea_objects.RData',sep=''))
}else{
  
  for(i in Universe.Sets){
    res=doGO(name=paste('Simplicity.GSEA.set',i, sep="."),
             gnId=Gene.Ann$ENTREZID[which(Gene.Ann$Set==i)], uni=universe,
             db=SP$DB,pvalue=1)
    #GERAR UM FDR
    
    save(res,file=paste('set_',i,'_gsea_objects.RData',sep=''))
  }#for i in Universe.Sets
  
}#ifelse is.null

#################################################################################
get.parent=function(g,family){
  p=NULL
  e=edges(g)
  #  j=1
  
  for(i in 1:length(e)){
    #   print(paste(j,names(e)[i],':',is.element(e[[i]],target)))
    #    j=j+1
    if(any(is.element(e[[i]],family))){
      # If the node i is has any family as edge, we include it as
      # family membe. if it wasn't already, it will be added to the
      # end of the string and have its own parents checked from now own
      family=unique(c(family,names(e)[i]))
    }
  }
  return(family)
}
gogetgene = function(l,gene){
  return(any(gene%in%l))
}

temp.graph=list(nodes=NULL,links=NULL)

for(I in 1:length(res)){
  ok=summary(res[[I]], pvalue=0.1)
  if(dim(ok)[1]>0){
    # selecting the graph element
    g=res[[I]]@goDag
    #exploring it
    degree(g)
    nodes(g)
    edges(g)
  
    plot(subGraph(ok[,1],g))
  
    # selecting the GOId with significative p-value
    target=ok[,1]
  
    # selecting the parent nodes for those targets
    g.size=length(target)
  
    family=get.parent(g,target)
    while(length(family)>length(target)){
      target=family
      family=get.parent(g,target)
    }
    target=ok[,1]
    
    # Creating table with the data for each family 'member'
    allgo=summary(res[[I]],pvalue=1.1)
    dim(allgo)
    #73 7
    length(family)
    # 11
    family.table= allgo[which(allgo[,1]%in%family),]
    dim(family.table)
    #[1] 11  7
    g1=subGraph(family,g)
  
    G=res[[I]]
  
    go.genes=NULL
    #The first version linked the Genes to all family GOs.
    for(i in which(names(g1@nodeData@data)%in%family)){
    # I will link only to the significative GOs
    #for(i in which(names(g1@nodeData@data)%in%target)){
      A=as.array(g1@nodeData@data[[i]]$geneIds)
      A=(A[which(g1@nodeData@data[[i]]$geneIds%in%G@geneIds)])
      A=Gene.Ann$SYMBOL[which(Gene.Ann$ENTREZID%in%A)]
      go.genes=c(go.genes, list(a=A))
      names(go.genes)[length(go.genes)]=names(g1@nodeData@data)[i]
    }
    go.genes
  
    geneIds=Gene.Ann$SYMBOL[which(Gene.Ann$ENTREZID%in%G@geneIds)]
    gene.edge=NULL
    for (i in geneIds){
      ggg=sapply(go.genes,gogetgene, gene=i)
      gene.edge=c(gene.edge, list(a=list(edges=names(ggg)[which(ggg)])))
      names(gene.edge)[length(gene.edge)]=i
    }
    
    test=addNode(node = geneIds, object = g1)
  
    for(i in geneIds){
      if(length(gene.edge[[i]]$edges)>0){
        test= addEdge(from = i, to = gene.edge[[i]]$edges, test)
      }
    }
    test.gos=grep('GO:', test@nodes)
  
     # color=rep('white',length(nodes(test)))
     # color[which(nodes(test)%in%family)]='green'
     # color[which(nodes(test)%in%target)]='blue'
     # color[which(nodes(test)%in%geneIds)]='yellow'
     # nA=makeNodeAttrs(test, label=gsub('GO:','',nodes(test)),
     #                  fixedsize=F,  width=20,
     #                  height=10, font=20, fontsize=30, fillcolor=color)
     # plot(test,nodeAttrs=nA)
     # 
  
    ################################
    # ############################ #
    # #                          # #
    # #      Creating  JSON      # #
    # #                          # #
    # ############################ #
    ################################
    
    
    sub.gene=Gene.Ann[which(Gene.Ann$SYMBOL%in%geneIds),-which(colnames(Gene.Ann)=='Set')]
    sub.gene$logFC=data1[sub.gene[,1],'logFC']
    sub.gene$pValue=data1[sub.gene[,1],'FDR']
    
    M=round(max(abs(sub.gene$logFC)),digits = 1)
    COL=matlab.like2(M*20+1)#diverge_hcl(183,c=100,l=c(50,90),power=1)#primary.colors((M*20+1), steps = 2, no.white = FALSE)
    names(COL)=((-M*10):(M*10))/10
    
    colnames(sub.gene)[which(colnames(sub.gene)=='SYMBOL')]='name'
    colnames(sub.gene)[which(colnames(sub.gene)=='GENENAME')]='longname'
    sub.gene$color=COL[as.character(round(sub.gene$logFC,digits = 1))]
    sub.gene$type='gene'
    
    temp.node=as.data.frame(test@nodes[test.gos])
    names(temp.node)='name'
    temp.node$color="#FFFFFF"
    temp.node$pValue=NA
    for(i in family){
      temp.node$pValue[which(temp.node$name==i)]=family.table$Pvalue[which(family.table[,1]==i)]
    }
    temp.node$type=colnames(family.table)[1]
    rownames(family.table)=family.table[,1]
    temp.node$longname=family.table[as.character(temp.node$name),'Term']
    temp.node$pValue=family.table[as.character(temp.node$name),'Pvalue']
    nodes=merge(sub.gene, temp.node, all=TRUE, sort=FALSE)
    
    rownames(nodes)=nodes$name
    
    # edges
    edges=test@edgeL
    
    #temp.edges=data.frame(col.names=c('source','target'))
    source=NULL
    #TARGET=target
    target=NULL
    for(i in 1:length(edges)){
      if(length(edges[[i]]$edges)){
        for(j in edges[[i]]$edges){
          source[length(source)+1]=names(edges)[i]
          target[length(target)+1]=test@nodes[j]
        }  
      }
    }
    temp.edges=data.frame(source,target)
    temp.edges$value=1
    
    temp.graph$nodes=rbind(temp.graph$nodes, nodes[setdiff(nodes$name,temp.graph$nodes$name),])
    temp.graph$links=rbind(temp.graph$links, temp.edges)
    #temp.graph=list(nodes, temp.edges)
  }
}


temp.graph$nodes=temp.graph$nodes[order(temp.graph$nodes$name),]
temp.graph$nodes=temp.graph$nodes[order(temp.graph$nodes$type),]
temp.graph$nodes$id=0:(length(temp.graph$nodes$name)-1)
levels(temp.graph$links$source)=temp.graph$nodes[levels(temp.graph$links$source),'id']
levels(temp.graph$links$target)=temp.graph$nodes[levels(temp.graph$links$target),'id']

row.names(temp.graph$nodes)=NULL

temp.json=toJSON(temp.graph)

write(temp.json,'mygraph.json')



