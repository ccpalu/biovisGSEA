###############################
###############################
# NSilico Life Sciences Ltda  #
# author: Cintia Palu         #
#                             #
# Consulting wich genome wide #
# annotation are available on #
# BioConductor                #
#                             #
# 25/04/2017                  #
# version: 1.1                #
###############################
###############################

# Bioconductor is a dynamic site that grows its database.
# It will be good to consult it and update the genome-wide libraries on Simplicity-tools
# This script aims to automatise this process.
# This libraries are used on gene annotation and on GSEA analysis
# More info available at https://bioconductor.org/help/workflows/annotation/annotation/

### Seting the folder
#local(setwd("C:\\Users\\cinti\\OneDrive\\Documents\\BioVis"))
#setwd ("/media/Dados/PosDoc/Cintia/analise")

######################################
######## Required Libraries###########
######################################

library("AnnotationHub")
ah <- AnnotationHub()
db.date=snapshotDate(ah)
gw = ah[which(ah$preparerclass=="OrgDbFromPkgsImportPreparer")]
#mcols(gw)

genomes=data.frame(cbind(gw$species,gw$title, gw$description, gw$rdatadateadded))
# genomes=data.frame(cbind(gw$species,gw$title, gw$description, gw$dataprovider,
#                          gw$taxonomyid, gw$genome, gw$maintainer, gw$rdatadateadded,
#                          gw$preparerclass, gw$tags, gw$rdataclass, gw$sourceurl, gw$sourcetype))
colnames(genomes)=c("Species", "DB", "Description","Date")


genomes$Species=as.character(genomes$Species)
genomes$DB=as.character(genomes$DB)
genomes$Description=as.character(genomes$Description)
genomes$Date=as.Date.character(genomes$Date)
genomes$DB=gsub(".sqlite", "", genomes$DB)

# Checking for duplicates
while(any(duplicated(genomes$Species))){
  dp=which(duplicated(genomes$Species))
  if(length(dp>1)){
    dps=which(genomes$Species==genomes$Species[dp[1]])
  }else{
    dps=which(genomes$Species==genomes$Species[dp])
  }
  for(i in dps){
    genomes$Species[i]=paste(genomes$Species[i]," (",genomes$DB[i],")", sep="")
  }
}

# Sorting genomes alphabeticaly
genomes=genomes[order(genomes$Species),]
head(genomes)  
# Species                    DB
# 1       Anopheles gambiae   org.Ag.eg.db.sqlite
# 2    Arabidopsis thaliana org.At.tair.db.sqlite
# 3              Bos taurus   org.Bt.eg.db.sqlite
# 15 Caenorhabditis elegans   org.Ce.eg.db.sqlite
# 4        Canis familiaris   org.Cf.eg.db.sqlite
# 18            Danio rerio   org.Dr.eg.db.sqlite
# Description       Date
# 1       NCBI gene ID based annotations about Anopheles gambiae 2015-08-26
# 2    NCBI gene ID based annotations about Arabidopsis thaliana 2015-08-26
# 3              NCBI gene ID based annotations about Bos taurus 2015-08-26
# 15 NCBI gene ID based annotations about Caenorhabditis elegans 2015-08-26
# 4        NCBI gene ID based annotations about Canis familiaris 2015-08-26
# 18            NCBI gene ID based annotations about Danio rerio 2015-08-26

###########################
# Downloading all libraries
# Once this is incorporated to Simplicity I recomend to always download the newest libraries
issue=array()
j=1

for(i in genomes$DB){
  
  tryCatch({
    library(i,character.only=TRUE)
     issue[j]=i
     j=j+1
  }, error= function (e){
    tryCatch({  
     source("http://bioconductor.org/biocLite.R")
      biocLite(i, suppressUpdates =  TRUE)
      library(i,character.only=TRUE)
    
      issue[j]=i
      j=j+1
    }, error=function (e){
      
    })
  }
 )
}

# Removing samples that do not have annotation
genomes= genomes[which(genomes$DB %in% issue),]

## I want to display the references genome on the top of the list
model.org=read.table("model_organisms.txt", sep="\t")
mo=which(genomes$Species %in% model.org$V1)

# THe following organisms do not have annotations
setdiff(model.org$V1, genomes$Species)
# [1] "Anolis carolinensis"       "Chlamydomonas reinhardtii"
# [3] "Dictyostelium discoideum"  "Emiliania huxleyi"        
# [5] "Escherichia coli"          "Fundulus heteroclitus"    
# [7] "Nothobranchius furzeri"    "Oryzias latipes"          
# [9] "Physcomitrella patens"     "Schizosaccharomyces pombe"
# [11] "Tetrahymena thermophila"   "Xenopus tropicalis"       
   

genomes[2:(length(genomes$Species)+1),]=genomes[c(mo,setdiff(1:length(genomes$DB),mo)),]
genomes[1,]=c("-----------","",paste("There are ", length(genomes$Species)-1,
                                     " genomes available, as provided by AnnotationHub", 
                                      sep=""), db.date)

# for(i in 1:length(genomes$Species)){
#   genomes[i,4]=paste("<em>",genomes$Species[i],"</em>", sep="")
# }
# colnames(genomes)[4]="Italic"
head(genomes)
# Species             DB
# 1001             -----------               
#   1002    Arabidopsis thaliana org.At.tair.db
# 1003  Caenorhabditis elegans   org.Ce.eg.db
# 1015             Danio rerio   org.Dr.eg.db
# 1004 Drosophila melanogaster   org.Dm.eg.db
# 1018            Homo sapiens   org.Hs.eg.db
# Description       Date
# 1001 There are 19 genomes available, as provided by AnnotationHub 2016-10-11
# 1002    NCBI gene ID based annotations about Arabidopsis thaliana 2015-08-26
# 1003  NCBI gene ID based annotations about Caenorhabditis elegans 2015-08-26
# 1015             NCBI gene ID based annotations about Danio rerio 2015-08-26
# 1004 NCBI gene ID based annotations about Drosophila melanogaster 2015-08-26
# 1018            NCBI gene ID based annotations about Homo sapiens 2015-08-26


rm(list=ls()[-which(ls()=='genomes')])
write.table(genomes, sep="\t", file="reference.genomes.txt")
save.image("reference.genomes.RData")
