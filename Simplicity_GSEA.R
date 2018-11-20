###############################
###############################
# NSilico Life Sciences Ltda  #
# author: Cintia Palu         #
#                             #
# Shiny document for inputting#
# data to do GSEA analysis    #
#                             #
# 19/11/2018                  #
# version: 0.10               #
#         under development   #
###############################
###############################


options(repos = BiocInstaller::biocinstallRepos())

library(shiny)
library(AnnotationDbi)

#I have to call the libraries so all of them are installed in the ShinyApp
#library(org.Mmu.eg.db);library(org.Cf.eg); library(org.Ag.eg); library(org.Xl.eg);
### Error on 1st Nov 2018 ###
# org.Mmu.eg contains GO mappings based on older data because the Blast2GO data resource was removed from the public domain just before the
# most recent update was produced. We are working on an alternative means to get this kind of data before the next release.
library(org.Pt.eg.db); library(org.Pf.plasmo.db); library(org.Rn.eg.db); library(org.Ss.eg.db)
library(org.At.tair.db); library(org.Ce.eg.db); library(org.Dr.eg.db); library(org.Dm.eg.db); library(org.Hs.eg.db);
library(org.Mm.eg.db); library(org.Sc.sgd.db); library(org.Xl.eg.db); library(org.Ag.eg.db); library(org.Bt.eg.db);
library(org.Cf.eg.db); library(org.EcK12.eg.db); library(org.EcSakai.eg.db); library(org.Gg.eg.db); 

source('gsea_gostats.R')

load("reference.genomes.RData")

ui <- fluidPage(
  titlePanel("GSEA - Gene Set Enrichment Analysis"),
  tags$a(href="http://www.cit.ie/",
         tags$img(src="www/cit.jpg", height='40'),
         'CIT - Cork Institute of Technology'),
  tags$a(href="http://research.ie/",
         tags$img(src="www/irc_long.jpg", height='50'),
         'IRC - Irish Research Council'),
  tags$a(href="http://nsilico.com/",
         tags$img(src="www/nsilico.png", height='50'),
         'NSilico - Simplifying Scientific Proecesses'),
  tags$a(href="https://www.sfi.ie/",
         tags$img(src="www/sfi.jpg", height='50'),
         'SFI - Science Foundation Ireland'),
  tags$a(href="http://www.ucc.ie",
         tags$img(src="./www/ucc.png", height='50'),
         'UCC - University College Cork'),
  tags$h3("Uploading Files"),
  
  ################## 
  ##### STEP 1 #####  
  ##################
  
  sidebarLayout(
    sidebarPanel(
    
    tags$h4("1 - Annotation"),
      selectInput(inputId = "species", label = "Select Organism:",
                  choices = genomes$Species)
      
    , width=5),#sidebarPanel
    mainPanel(
      tags$hr(),
      textOutput("db.info"),
      tags$hr()
      , width=7)
  ),#sidebarLayout
  
  ################## 
  ##### STEP 2 #####  
  ##################
  
  sidebarLayout(
    sidebarPanel(
      tags$h4("2 - Gene set"),
      fileInput('geneSet.file', 'Choose CSV or tab-delimited File',
                accept=c('text/csv', 
                         'text/comma-separated-values,text/plain', 
                         '.csv')),
   
   #########################################
      tags$hr(),
      tags$h5(tags$em(textOutput("table.info"))),
      tags$hr(),
   
      checkboxInput('header', 'Header', TRUE),
      radioButtons('sep', 'Separator',
                   c(Comma=',',
                     Semicolon=';',
                     Tab='\t'),
                   ',',inline = TRUE),
      radioButtons('quote', 'Quote',
                   c(None='',
                     'Double Quote'='"',
                     'Single Quote'="'"),
                   '"',inline = TRUE),
   #########################################
    tags$hr(),
    uiOutput('colnames.gene'),
    radioButtons('gene.id', 'Gene ID type:',
                c('Symbol'='symbol',
                  'ENTREZID'="entrez",
                  'ENSEMBL'="ensembl",
                  'PMID'="pubmed",
                  'RefSeq'="refseq",
                  'UniGene'="unigene"),
                inline = TRUE),
    tags$h5(tags$em(textOutput('test.gene.id'))),
     radioButtons('collection', 'How many gene sets there are in your file?',
                  c('One'="set",'Multiple'='universe'),
                  'set',inline = TRUE),
     uiOutput('colnames.set'),
     uiOutput('total.set')
    
   , width=5),#sidebarPanel
    mainPanel(
      fluidRow(
        tags$hr(),
        tags$h5(tags$em("Before uploading the gene list, save", tags$strong("downregulated"),
                      "genes in a file", tags$strong("apart"), "from the ", tags$strong("upregulated"),
                      "ones. They should be", tags$strong("analysed separely.")))
      ),
      fluidRow(
        tags$hr(),
        dataTableOutput('dataset')
        )
      , width=7)#mainPanel
   
  ),#sidebarLayout

  ##################
  ##### STEP 3 #####
  ##################
  
  sidebarLayout(
    sidebarPanel(
      tags$h4("3 - Gene 'universe'"),
      #uiOutput('step3'),
      uiOutput('universe.source'),
      uiOutput('geneUniverse.file'),
      #########################################
      tags$hr(),

      uiOutput('header.universe'),
      uiOutput('sep.universe'),
      uiOutput('quote.universe'),
      #########################################
      tags$hr(),
      uiOutput('colnames.gene.universe'),
      uiOutput('gene.id.universe')
      , width=5),#sidebarPanel
    mainPanel(
      
      fluidRow(
        tags$hr(),
        tags$p("The gene set you submited will be compared to a larger gene collection,
               which can be:"),
        tags$div(tags$ol(
          tags$li(tags$span("all genes present on the annotation files,")),
          tags$li(tags$span("all genes present on the uploaded file in case of multiple datasets, or")),
          tags$li(tags$span("a new gene list provided by you.")))),
        tags$p()
        ),
      fluidRow(
        tags$hr(),
        dataTableOutput('dataset.universe')
        )
      , width=7)

  ),#sidebarLayout
  
  ##################
  ##### STEP 4 #####
  ##################
  wellPanel(
    tags$h4("4 - Input validation"),
    fluidRow(
      column(2,tags$strong("Step 1:")),#column
      column(10,textOutput('validate.step1')),
      tags$p()
    ),#fluidRow
    fluidRow(
      column(2,tags$strong("Step 2:")),#column
      column(10,textOutput('validate.step2')),
      tags$p()
    ),#fluidRow
    fluidRow(
      column(2,tags$strong("Step 3:")),#column
      column(10,textOutput('validate.step3')),
      tags$p()
    ),#fluidRow
    fluidRow(
      column(4),
      column(8,uiOutput("enable.gsea"))#column
    )#fluidRow
  ),#wellPanel
  
  ##################
  ##### STEP 5 #####
  ##################
  #https://shiny.rstudio.com/articles/generating-reports.html
  #http://shiny.rstudio.com/gallery/download-knitr-reports.html
  wellPanel(
    tags$h4("5 - Results"),
    uiOutput('resultsInTabs')
    
  )#wellPanel
      
)

server <- function(input, output, session) {
  #####################
  ##### Functions #####
  #####################
  
  
  ###############
  ##### Var #####
  ###############
  
  # Species selected on Step 1
  sp <- reactive({
    which(genomes$Species == input$species)
  })
  
  # Table uploaded on Step 2
  data1 <- reactive({
    #validate(need(input$file1, "Please upload a file"))
    inFile <- input$geneSet.file
    
    if (is.null(inFile))
      return(NULL)
    
    read.csv(inFile$datapath, header=input$header, sep=input$sep,
             quote=input$quote)
  })#data1
  
  
  # Gene set list
  universe.sets=reactive({
    if(length(input$header.set) && input$header.set!=""){
      return(unique(data1()[,input$header.set]))
    }else return(NULL)
  })
  
  # Gene list 
  gene.list <- reactive({
    if(input$header.gene !=""){
      return(as.character(unique(data1()[,input$header.gene])))
    }else return(NULL)
  })
  # Table Uploaded on Step3
  data2 <- reactive({
    
    inFile <- input$geneUniverse.file
    
    if (is.null(inFile) || (input$universe.source)!='new.source'){
      return(NULL)
    }else{
      read.csv(inFile$datapath, header=input$header.universe, sep=input$sep.universe,
               quote=input$quote.universe)
    }
  })#data2
  
  # Gene Universe list 
  gene.list.universe <- reactive({
    if(input$header.gene.universe !=""){
      return(as.character(unique(data2()[,input$header.gene.universe])))
    }else return(NULL)
  })
  
  # Collecting gene annotation
  gene.ann <- reactive({
    
    switch(input$gene.id,
    "symbol" = select(get(genomes$DB[sp()]), keys=gene.list(),
                                 columns= c("ENTREZID", "SYMBOL", "GENENAME"),
                                keytype="SYMBOL"),
    "entrez" = select(get(genomes$DB[sp()]), keys=gene.list(),
                                  columns= c("ENTREZID", "SYMBOL", "GENENAME"),
                                  keytype="ENTREZID"),
    "ensembl" = select(get(genomes$DB[sp()]), keys=as.character(gene.list()),
                                   columns= c("ENTREZID", "SYMBOL", "GENENAME","ENSEMBL"),
                                   keytype="ENSEMBL"),
    "pubmed" = select(get(genomes$DB[sp()]), keys=as.character(gene.list()),
                                  columns= c("ENTREZID", "SYMBOL", "GENENAME","PMID"),
                                  keytype="PMID"),
    "refseq" = select(get(genomes$DB[sp()]), keys=as.character(gene.list()),
                                  columns= c("ENTREZID", "SYMBOL", "GENENAME","REFSEQ"),
                                  keytype="REFSEQ"),
    "unigene" = select(get(genomes$DB[sp()]), keys=as.character(gene.list()),
                                   columns= c("ENTREZID", "SYMBOL", "GENENAME, UNIGENE"),
                                   keytype="UNIGENE"),
    ""
    )
  })#gene.ann()  
  gene.ann.universe <- reactive({
     
    switch(input$gene.id.universe,
           "symbol" = select(get(genomes$DB[sp()]), keys=gene.list.universe(),
                             columns= c("ENTREZID", "SYMBOL", "GENENAME"),
                             keytype="SYMBOL"),
           "entrez" = select(get(genomes$DB[sp()]), keys=gene.list.universe(),
                             columns= c("ENTREZID", "SYMBOL", "GENENAME"),
                             keytype="ENTREZID"),
           "ensembl" = select(get(genomes$DB[sp()]), keys=gene.list.universe(),
                              columns= c("ENTREZID", "SYMBOL", "GENENAME","ENSEMBL"),
                              keytype="ENSEMBL"),
           "pubmed" = select(get(genomes$DB[sp()]), keys=gene.list.universe(),
                             columns= c("ENTREZID", "SYMBOL", "GENENAME","PMID"),
                             keytype="PMID"),
           "refseq" = select(get(genomes$DB[sp()]), keys=gene.list.universe(),
                             columns= c("ENTREZID", "SYMBOL", "GENENAME","REFSEQ"),
                             keytype="REFSEQ"),
           "unigene" = select(get(genomes$DB[sp()]), keys=gene.list.universe(),
                              columns= c("ENTREZID", "SYMBOL", "GENENAME, UNIGENE"),
                              keytype="UNIGENE"),
          ""
    )
  })#gene.ann.universe() 
  
 ################## 
 ##### STEP 1 #####  
 ##################
  
  output$db.info<-renderText({
    if(sp()!=1){
      
      tryCatch({
        library(genomes$DB[sp()], character.only = TRUE)
        },error = function(e){
            e
        })
    }
    paste(as.character(genomes$Description[sp()]), ' in ', 
          as.character(genomes$Date[sp()]),".", sep="")
    
  })#output$db.info
  
  
  ################## 
  ##### STEP 2 #####  
  ##################
  
  output$dataset <- renderDataTable(data1(),options = list(
    pageLength = 10))# output$dataset
   
  output$colnames.gene <- renderUI({
    selectInput(inputId = "header.gene", label = "Gene IDs column:",
                choices = c("",colnames(data1())))
  })
  output$test.gene.id <-renderText({
    if(sp()==1){
      paste("Please select the organism on Step 1") 
    }else{
      
      if(length(data1())>1){#Needs the gene file
        if(input$header.gene==""){
          paste("Please indicate the column with the gene IDs and verify if the gene ID type is correct")
        }else{
          tryCatch({
               na=length(which(is.na(gene.ann()[,'ENTREZID'])))
               paste("We found annotation for",  dim(gene.ann())[1] - na, 'genes. There are', na, 'unrecognised ones.')
          },error=function(e){
            paste("Simplicity could not recognise any gene annotation.
                  Try choosing a different species, gene column or gene ID type.")
          })
        }
      }else{paste('Upload a file on Step 2')}
    }
  })
  output$colnames.set <- renderUI({
      if(input$collection=='universe'){
        opt = colnames(data1())
        selectInput(inputId = "header.set", label = "Set IDs column:",
                    choices = c("",opt[-which(opt==input$header.gene)]))
      }else return(NULL)
  })
  output$total.set <- renderUI({
    if(length(universe.sets())){
      paste("We have identified", length(universe.sets()), "gene sets.")
    }else return(NULL)
  })
 
  
  #################
  #### STEP 3 #####
  #################
    
  output$universe.source <- renderUI({
    if(length(data1())){
      if(input$collection=='universe'){
        radioButtons(inputId = 'universe.source', label = 'Source:',
                     choices = c('NCBI Annotation'='ncbi','Step 2 table'="current",
                     'New source'='new.source'), selected = "current",inline = TRUE)
      }else {
        radioButtons('universe.source', 'Gene collection source',
                     choices= c('NCBI Annotation'='ncbi',
                       'New source'='new.source'),selected = "ncbi",inline = TRUE)
      }
    }else {paste("This step will be enabled after a file is uploaded on Step 2.")}
  })
  
  output$table.info <- renderText({if(length(data1())){paste("Please, verify if we were able to recognise the 
                                                             columns correctly. If there any issue, try to 
                                                             change the options below.")}
  })
  
  output$geneUniverse.file <- renderUI({
    if(length(data1()) && length(input$universe.source)){ 
      if(input$universe.source=='new.source'){
        fileInput('geneUniverse.file', 'Choose CSV or tab-delimited File',
                  accept=c('text/csv', 'text/comma-separated-values,text/plain',
                         '.csv'))
        }else {
          return(NULL)
        }
      }else {
        return(NULL)}
  })
  
  output$header.universe <- renderUI({
    if(length(data1()) && length(input$universe.source)){
      if(input$universe.source=='new.source'){
        checkboxInput('header.universe', 'Header', input$header)
      }else return(NULL)
    }else {return(NULL)}
  })
  output$sep.universe <- renderUI({
    if(length(data1()) && length(input$universe.source)){
      if(input$universe.source=='new.source'){
        radioButtons('sep.universe', 'Separator',
                 c(Comma=',',
                   Semicolon=';',
                   Tab='\t'),
                 input$sep,inline = TRUE)
      }else {return(NULL)}
    }else {return(NULL)}
  })
  output$quote.universe <- renderUI({
    if(length(data1()) && length(input$universe.source)){
      if(input$universe.source=='new.source'){
        radioButtons('quote.universe', 'Quote',
                 c(None='',
                   'Double Quote'='"',
                   'Single Quote'="'"),
                 input$quote,inline = TRUE)
      }else {return(NULL)}
    }else {return(NULL)}
  })

  output$gene.id.universe <- renderUI({
    if(length(data1()) && length(input$universe.source)){
      if(input$universe.source=='new.source'){
        radioButtons('gene.id.universe', 'Gene ID type:',
                   c('Symbol'='symbol',
                      'ENTREZID'="entrez",
                      'ENSEMBL'="ensembl",
                      'PMID'="pubmed",
                      'RefSeq'="refseq",
                      'UniGene'="unigene"),
                   input$gene.id,inline = TRUE)
      }else {return(NULL)}
    }else {return(NULL)}
  })

 
  output$dataset.universe <- renderDataTable(data2(),
                          options = list(pageLength = 10))# output$dataset.universe
  
  ########################

  output$colnames.gene.universe <- renderUI({
    if(length(data2())){
      if(input$universe.source=='new.source'){
        selectInput(inputId = "header.gene.universe", label = "Gene IDs column:",
                    choices = c("",colnames(data2())))
        }else {return(NULL)}
    }else {return(NULL)}
  })
  
  output$gene.id.universe <- renderUI({
    if(length(data2()) && input$universe.source=='new.source'){
      radioButtons('gene.id.universe', 'Gene ID type:',
               c('Symbol'='symbol',
                  'ENTREZID'="entrez",
                  'ENSEMBL'="ensembl",
                  'PMID'="pubmed",
                  'RefSeq'="refseq",
                  'UniGene'="unigene"),
               'symbol',inline = TRUE)
    }
  })
  
  ##################
  ##### STEP 4 #####
  ##################
  
  output$validate.step1 <- reactive({
    validate(
      need(sp()!=1, "Please, select a species")
    )#validate
    paste("OK")
  })#output$validate.step1
  
  output$validate.step2 <- reactive({paste(step2())})
  
  step2 <- reactive({
    
   if(length(data1())>0){
      
    tryCatch({
      na=(which(is.na(gene.ann()[,'ENTREZID'])))
      if(length(na)>0){
        if(length(gene.ann()[-na,'ENTREZID'])<2){
           
           return('There is not enough annotated genes')
        }else{
          tryCatch({
            if(input$collection=='set'){
               
               return('OK')
            }else{
              if(length(universe.sets())>1){
                 
                 return("OK")
              }else{
                 
                 return("Are you sure you have more than one gene set in this file? 
                      Please, check the column or type of gene ID you selected.")
              }
            }
          }, error=function(e){
             
             return("Are you sure you have more than one gene set in this file? 
                      Please, check the column or type of gene ID you selected.")
          })
        }
      }else{
         
         return("We need information regarding the gene annotation in the inputed file.")
      }
    }, error=function(e){
       
       return("We need information regarding the gene annotation in the inputed file.")
    }
    )
   }else{
     return("Please, upload gene set(s)")
   }
  })#output$validate.step2

  step2.status <- reactive({
    return(step2()=='OK')
  })

  step3.status <- reactive({

    return(length(grep("We found annotation for",step3()))>0 || step3()=='OK')
  })

  output$validate.step3= reactive({paste(step3())})
  step3 <- reactive({
    tryCatch({
      
      if (length(gene.ann())){
        
        if(input$universe.source=="new.source"){
          tryCatch({
            if(length(data2())!=0){
               if (is.null(input$header.gene.universe) ||
                   input$header.gene.universe == ""){
                  return("Select the gene column") 
               }else{
                 
                 tryCatch({
                   uni= gene.ann.universe()[,'ENTREZID']
                   na= which(is.na(uni))
                   set= gene.ann()[,'ENTREZID']
                   na.set=which(is.na(set))
                   tryCatch({
                     lists.intersect =intersect(set,uni[-na])
                      
                     if (length(lists.intersect)==length(set[-na.set])){
                        return(paste("OK. We found annotation for",  length(uni)-length(na), 'genes and there are', 
                             length(na), 'not recognised. Nonetheless, all the', length(set[-na.set]), 'genes identified on step 2 are present.
                             We are ready to start your analysis.'))}
                     else if (length(lists.intersect)){
                        return(paste("OK. We found annotation for",  length(uni)-length(na), 'genes and there are', 
                             length(na), "not recognised.
                             The gene 'universe' contains only", length(lists.intersect), "out of the", length(set[-na.set]),
                             'genes identified on step 2. This means that the following gnes will not be used in the GSEA:',
                             setdiff(set, lists.intersect),". You may consider choosing a different source."))
                     }else{
                        
                        return("Simplicity could not find annotation for the genes, please change the parameters or choose another source.")
                     }
                     }, error= function(e){
                        
                        return("Review the information on Step 3")
                     })
                   
                 },error=function(e){
                    
                    return("Please verify the required parameters")
                 })
               }
            }else{
               
               return('Upload a file or choose other source')
            }
          },error = function(e){
             
             return("Please, check this step's parameters")
          })
        }else{
           
            return('OK')
        }
      }else{
         
         return('Please, complete step 2') 
      }
    }, error = function(e){
       
       return('Please, complete step 2')
    })
    
  })#output$validate.step3
  
  output$enable.gsea <- renderUI({
  
    tryCatch({
      if(step2.status() && step3.status()){
        actionButton("start.gsea", "Start GSEA")
      tryCatch({
        if(gene.ann()==""){
          paste("Please, complete all steps to enable GSEA analysis.")
        }else{
          actionButton("start.gsea", "Start GSEA")
        }}, error = function(e){
          paste("Please, complete all steps to enable GSEA analysis.")
        }
      )
    }else{
      paste("Complete all steps to enable GSEA analysis.")
    }}, error = function(e){
      paste("Complete all steps to enable GSEA analysis.")
      }
    )
    
    ###### Include the analysis option MP, CC, BP, p-valor
  })
  
  ##################
  ##### STEP 5 #####
  ##################
  
  observeEvent(input$start.gsea, {
    SP = genomes[sp(),]
    Universe.Sets = universe.sets()
    
    if(input$collection=='universe'){
      Data1=data1()
      Data1=Data1[,c(input$header.gene, input$header.set)]
      
      switch(input$gene.id,
                     "symbol" = {colnames(Data1)=c('SYMBOL', 'Set')},
                       
                     "entrez" = {colnames(Data1)=c('ENTREZID', 'Set')},
                  
                     "ensembl" = {colnames(Data1)=c('ENSEMBL', 'Set')},
                                        
                     "pubmed" = {colnames(Data1)=c('PMID', 'Set')},
                                       
                     "refseq" = {colnames(Data1)=c('REFSEQ', 'Set')},
                                       
                     "unigene" = {colnames(Data1)=c('UNIGENE', 'Set')}                        
              )
      
      Gene.Ann = merge(gene.ann(), Data1,by=colnames(Data1)[1])
      rm(Data1)
    }else{
      Gene.Ann = gene.ann()
      Gene.Ann[,length(colnames(Gene.Ann))+1]='NA'
      colnames(Gene.Ann)[length(colnames(Gene.Ann))]='Set'
    }
    
    Gene.Ann=Gene.Ann[-which(is.na(Gene.Ann$ENTREZID)),]
    
    Universe.Source = input$universe.source
    
    if(input$universe.source=='new.source'){
      Universe.Gene.Ann = gene.ann.universe()
      Universe.Gene.Ann = Universe.Gene.Ann[-which(is.na(Universe.Gene.Ann$ENTREZID)),]
    }else{
      Universe.Gene.Ann=NULL
    }
    
    result <- simplicity_gostats(Gene.Ann = Gene.Ann, SP = SP, Universe.Sets = Universe.Sets,
      Universe.Source = Universe.Source, Universe.Gene.Ann =  Universe.Gene.Ann, pvalue = 0.05,
      output)
    
    output$resultsInTabs <- renderUI({
      #resultTabs = lapply(paste0(c('bp','cc','mf'), rep(1:15,each=3)), tabPanel)
      if(is.null(universe.sets())){
        resultTabs = list(tabPanel("Biological Process", dataTableOutput('bp1')),
          tabPanel("Cellular Component", dataTableOutput('cc1')),
          tabPanel("Molecular Function", dataTableOutput('mf1')))
      }else{
        resultTabs <- list()
        for(i in universe.sets()){
          resultTabs <- c(resultTabs,list(tabPanel(paste("BP",i), dataTableOutput(paste0('bp',i))),
          tabPanel(paste("CC", i), dataTableOutput(paste0('cc',i))),
          tabPanel(paste("MF", i), dataTableOutput(paste0('mf',i)))))
        }
      }#if/else length(resultTabs)==0
      
      do.call(tabsetPanel, args = resultTabs, quote = TRUE)
    })#output$resultsInTabs
  
    for(i in names(result)){
      output[[i]] <- renderDataTable(result[[i]], escape = FALSE)
    }
    # 
    # save(list = ls(all.names = TRUE),
    #       file = sprintf("%s_%s.RData","GSEA", as.integer(Sys.time())))
    # 
  })#observeEvent(input$start.gsea
  # Add loading, try to include the results 'we we go', add download button
 
}#server

shinyApp(ui = ui, server = server)
