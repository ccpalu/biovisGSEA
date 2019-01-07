###############################
###############################
# NSilico Life Sciences Ltda  #
# author: Cintia Palu         #
#                             #
# Shiny document for inputting#
# data to do GSEA    #
#                             #
# 04/01/2019                  #
# version: 0.2               #
#         under development   #
###############################
###############################


options(repos = BiocInstaller::biocinstallRepos())

library(shiny)
library(AnnotationDbi)
library(jsonlite)
library(colorRamps)
library(GOstats)
library(Category)
library(Rgraphviz)
#library(DT)

#I have to call the libraries so all of them are installed in the ShinyApp
#library(org.Cf.eg); library(org.Ag.eg); library(org.Xl.eg);
### Error on 1st Nov 2018 ###
# org.Mmu.eg contains GO mappings based on older data because the Blast2GO data resource was removed from the public domain just before the
# most recent update was produced. We are working on an alternative means to get this kind of data before the next release.
library(org.Pt.eg.db); library(org.Pf.plasmo.db); library(org.Rn.eg.db); library(org.Ss.eg.db)
library(org.At.tair.db); library(org.Ce.eg.db); library(org.Dr.eg.db); library(org.Dm.eg.db); library(org.Hs.eg.db);
library(org.Mm.eg.db); library(org.Sc.sgd.db); library(org.Xl.eg.db); library(org.Ag.eg.db); library(org.Bt.eg.db);
library(org.Cf.eg.db); library(org.EcK12.eg.db); library(org.EcSakai.eg.db); library(org.Gg.eg.db); library(org.Mmu.eg.db);

source('gsea_gostats.R')

load("reference.genomes.RData")
foo = 10;

ui <- fluidPage(
  titlePanel("GSEA - Gene Set Enrichment Analysis"),
  tags$a(href="http://www.cit.ie/",
    #img(src="www/cit.jpg", height='40'),
    '|CIT - Cork Institute of Technology|'),
  tags$a(href="http://research.ie/",
    #img(src="www/irc_long.jpg", height='50'),
    '|IRC - Irish Research Council|'),
  tags$a(href="http://nsilico.com/",
    #tags$img(src="www/nsilico.png", height='50'),
    '|NSilico - Simplifying Scientific Proecesses|'),
  tags$a(href="https://www.sfi.ie/",
    #tags$img(src="www/sfi.jpg", height='50'),
    '|SFI - Science Foundation Ireland|'),
  tags$a(href="http://www.ucc.ie",
    #tags$img(src="./www/ucc.png", height='50'),
    '|UCC - University College Cork|'),
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
      column(2, tags$strong("Step 1:")),#column
      column(10, textOutput('validate.step1')),
      tags$p()
    ),#fluidRow
    fluidRow(
      column(2, tags$strong("Step 2:")),#column
      column(10, textOutput('validate.step2')),
      tags$p()
    ),#fluidRow
    fluidRow(
      column(2, tags$strong("Step 3:")),#column
      column(10, textOutput('validate.step3')),
      tags$p()
    ),#fluidRow
    fluidRow(
      column(4),#  actionButton("generategraph", "Generate Graph")),
      column(8, uiOutput("enable.gsea"))#column
    )#fluidRow
  ),#wellPanel
  
  ##################
  ##### STEP 5 #####
  ##################
  #https://shiny.rstudio.com/articles/generating-reports.html
  #http://shiny.rstudio.com/gallery/download-knitr-reports.html
  wellPanel(
    tags$h4("5 - Results in Tables"),
    uiOutput('resultsInTabs')
  ),#wellPanel
  
  ##################
  ##### STEP 6 #####
  ##################
  sidebarLayout(
    sidebarPanel(
      tags$h4("6 - Interactive Results"),
      tags$hr(),
      fluidRow(
        column(3,
          tags$p(),
          tags$div(align="center", 
            tags$h5(tags$b("Graph Control")),
            tags$hr(),
            HTML('
              <button id="graphunlink"style="color:grey;width:120px" disabled = true>Link/Unlink</button><br><br>
              <button id="graphstatic" style="color:grey;width:120px" disabled = true>Freeze/Unfreeze</button><br><br>
              <h6>p-value threshold</h6>
              <input type="number" id="pValue" value="0.1" step="0.01" max="0.1" min="0.0000000000" style="width:120px;">
              <button id="graphfilterp" style="color:grey;width:120px" disabled = true>Filter p-value</button><br><br>
              <input type="checkbox" id="filterGeneNodes" disabled = true"> Remove disconected genes<br><br>
              ')),
          tags$hr()
            ),#column
        column(1),
        column(8,
          tags$div(align = "center",
            tags$h5(tags$b("Legend")),
            
            tags$hr(),
            HTML('
              <svg id="legend-svg"></svg><br><br>
              <svg id="legendGOID-svg"></svg>
              '),#HTML
            tags$hr()
            )#div
          )#column  
          ),#fluidRow
      width = 5
        ),#sidebarPanel
    #tags$hr(),
    #tags$head(tags$script(src = "d3.v3.js"))
    mainPanel(
      tags$div(#includeHTML('final.html')
        HTML('
          <head>
          <script src="https://d3js.org/d3.v3.min.js"></script>
          <!--<script type="text/javascript" src="d3.v3.js"></script> -->
          <style type="text/css">
          .link {
          stroke: #ccc;
          }
          .link.selected{
          stroke: red;
          }
          
          .nodetext {
          pointer-events: none;
          font: 10px sans-serif;
          }
          
          /* Tooltip container */
          div.tooltip {   
          position: absolute;				 
          text-align: center;				 
          width: auto;										
          height: auto;								   
          padding: 2px;						   
          font: 12px sans-serif;		  
          background: lightsteelblue; 
          border: 0px;				
          border-radius: 8px;				 
          pointer-events: none;
          margin-left:30px;
          }
          
          /* Legend */
          .axis path,
          .axis line {
          fill: none;
          stroke: #000;
          shape-rendering: crispEdges;
          }
          
          circle {
          stroke-width: 1;
          stroke: white;
          }
          
          #legend-svg {
          vertical-align: bottom;
          }
          
          #controls {
          padding-top: 10px;
          }
          </style>
          </head>
          
          <body>
          
          <script type="text/javascript">
          //Loading json file
          // Error: Cross origin requests are only supported for protocol schemes: http, data, chrome, chrome-extension, https
          /*d3.json("graph.json", function(error, json) {
          if (error) throw error;
          
          root = json;
          update();
          });*/
          var maxlogFC = 5;
          var pValue = 0.1;
          var dim = 32;
          
          /* Gene Legend 
          adpted from https://bl.ocks.org/starcalibre/6cccfa843ed254aa0a0d*/
          var legendFullHeight = 60;
          var legendFullWidth = 320;
          // Margins
          var legendMargin = { top: 5, bottom: 25, left: 105, right: 5 };
          var legendWidth = legendFullWidth - legendMargin.left - legendMargin.right;
          var legendHeight = legendFullHeight - legendMargin.top - legendMargin.bottom;
          
          
          // color range used for  logFC
          var scale = ["#0000BF", "#0000CF", "#0000DF", "#0000EF", "#0000FF", "#0020FF", "#0040FF",		  								"#0060FF", "#0080FF", "#009FFF","#00BFFF", "#00DFFF", "#00FFFF", "#10FFEF",
          "#20FFDF", "#30FFCF", "#40FFBF", "#50FFAF", "#60FF9F", "#70FF8F","#80FF80",
          "#8FFF70", "#9FFF60", "#AFFF50", "#BFFF40", "#CFFF30", "#DFFF20", "#EFFF10",
          "#FFFF00", "#FFDF00","#FFBF00", "#FF9F00", "#FF8000", "#FF6000", "#FF4000",
          "#FF2000", "#FF0000", "#EF0000", "#DF0000", "#CF0000", "#BF0000"];
          
          var legendSvg = d3.select("#legend-svg")
          .attr("width", legendFullWidth)
          .attr("height", legendFullHeight)
          .append("g")
          .attr("transform", "translate(" + legendMargin.left + "," +
          legendMargin.top + ")");
          
          //Coloring Genes and adding rectangle//#aqui
          legendSvg.append("rect")
          .attr("x", -100)
          .attr("y", 0)
          .attr("width", dim*2.5)
          .attr("height", dim)
          .style("fill", "#b3b3b3");
          legendSvg.append("svg:text")
          .attr("class", "nodetext")
          .attr("dx", -97)
          .attr("dy", dim/2 - 3)
          .style("font-size", 12)
          .text("Gene Symbol");
          legendSvg.append("svg:text")
          .attr("class", "nodetext")
          .attr("dx", -97)
          .attr("dy", dim-5)
          .style("font-size", 10)
          .text("DGE p-value");
          
          var gradient = legendSvg.append("defs")
          .append("linearGradient")
          .attr("id", "gradient")
          .attr("x1", "0%") // bottom
          .attr("y1", "0%")
          .attr("x2", "100%") // to top
          .attr("y2", "0%")
          .attr("spreadMethod", "pad");
          
          var pct = linspace(0, 100, scale.length).map(function(d) {
          return Math.round(d) + "%";
          });
          
          var colourPct = d3.zip(pct, scale);
          colourPct.forEach(function(d) {
          gradient.append("stop")
          .attr("offset", d[0])
          .attr("stop-color", d[1])
          .attr("stop-opacity", 1);
          });
          
          legendSvg.append("rect")
          .attr("x1", legendMargin.left)
          .attr("y1", legendMargin.top)
          .attr("width", legendWidth)
          .attr("height", legendHeight)
          .style("fill", "url(#gradient)");
          
          // create a scale and axis for the legend
          var legendScale = d3.scale.linear()
          .domain([maxlogFC, -maxlogFC])
          .range([legendWidth, 0]);
          
          var legendAxis = d3.svg.axis()
          .scale(legendScale)
          .orient("bottom")
          .tickValues(d3.range(-maxlogFC, (maxlogFC + 1)))
          .tickFormat(d3.format("d"));
          
          legendSvg.append("g")
          .attr("class", "legend axis")
          .attr("transform", "translate(0," + legendHeight + ")")
          .call(legendAxis);
          
          legendSvg.append("svg:text")
          .attr("class", "nodetext")
          .attr("x", legendWidth/2 - 85)
          .attr("dy", legendHeight/2 + 5)
          .style("font-size", 12)
          .text("Gene color based on the logFC");
          /* end of Genes Legend */
          
          /* GO legend */
          var rad = 20;
          var gap = 25;
          var mar = 180;
          
          var GOtext = [
          { "x_axis": 23, "y_axis": 30, "text": "Gene Ontology Domain:", "id": "L.god"},
          { "x_axis": 28, "y_axis": gap*2 + rad -5, "text": "Cellular Component", "id": "L.cc"},
          { "x_axis": 60, "y_axis": gap*2 + rad + 12, "text": "(GOCCID)", "id": "L.gocc"},
          { "x_axis": 33, "y_axis": rad*2 + gap*3 - 5, "text": "Biological Process", "id": "L.bp"},
          { "x_axis": 60, "y_axis": rad*2 + gap*3  + 12, "text": "(GOBPID)", "id": "L.gobp"},
          { "x_axis": 28, "y_axis": rad*3 + gap*4 -5, "text": "Molecullar Function", "id": "L.mf"},
          { "x_axis": 60, "y_axis": rad*3 + gap*4 + 12, "text": "(GOMFID)", "id": "L.gomf"},
          { "x_axis": mar + (gap+rad), "y_axis": 20 , "text": "p-value", "id": "L.pv"},
          { "x_axis": mar + gap - 4, "y_axis": gap*1.5 , "text": "> " + pValue, "id": "L.Bpv"},
          { "x_axis": mar + gap*3 - 4, "y_axis": gap*1.5 , "text": "<= " + pValue, "id": "L.Spv"}
          ];
          
          var jsonCircles = [
          { "x_axis": gap + rad/2, "y_axis": gap*2 + rad, "radius": rad, "color" : "#ffe44d" },
          { "x_axis": (gap + rad)*2, "y_axis": gap*2 + rad, "radius": rad, "color" : "#ccae00" },
          { "x_axis": gap + rad/2, "y_axis": gap*3 + rad*2, "radius": rad, "color" : "#ceb6d1" },
          { "x_axis": (gap + rad)*2, "y_axis": gap*3 + rad*2, "radius": rad, "color" : "#a377a9" },
          { "x_axis": gap + rad/2, "y_axis": gap*4 + rad*3, "radius": rad, "color" : "#45ff9a" },
          { "x_axis": (gap + rad)*2, "y_axis": gap*4 + rad*3, "radius": rad, "color" : "#00c45a" }
          ];
          
          var legendSvg = d3.select("#legendGOID-svg")
          .attr("width", mar + (rad + gap) * 3)
          .attr("height", (gap + rad) * 4);
          
          var text = legendSvg.selectAll("text")
          .data(GOtext)
          .enter()
          .append("text");
          
          var textAttributes = text
          .attr("x", function (d) { return d.x_axis; })
          .attr("y", function (d) { return d.y_axis; })
          .text(function (d) { return d.text; })
          .attr("id", function (d) { return d.id; });
          
          var circles = legendSvg.selectAll("circle")
          .data(jsonCircles)
          .enter()
          .append("circle");
          
          var circleAttributes = circles
          .attr("cx", function (d) { return mar + d.x_axis; })
          .attr("cy", function (d) { return d.y_axis; })
          .attr("r", function (d) { return d.radius; })
          .style("fill", function(d) { return d.color; });
          
          legendSvg.selectAll("circle")
          .data(jsonCircles)
          .enter()
          .append("circle");
          /* end of GO legend */
          
          
          function linspace(start, end, n) {
          /* Support the creation of a colour scale 
          From https://bl.ocks.org/starcalibre/6cccfa843ed254aa0a0d */
          var out = [];
          var delta = (end - start) / (n - 1);
          
          var i = 0;
          while(i < (n - 1)) {
          out.push(start + (i * delta));
          i++;
          }
          
          out.push(end);
          return out;
          }
          
          var graph;  
          Shiny.addCustomMessageHandler("SendObjectToClientDynamicCallbackHandler",  
          function(variables) {  
          graph = JSON.parse(variables);
          function freeze(){
          graph.nodes.forEach(function(el) {
          el.fixed = !el.fixed;
          })
          };
          
          graph.links.forEach(function(el) {
          delete el.value;
          el.source = parseInt(el.source);
          el.target = parseInt(el.target);
          });
          
          var str = JSON.stringify(graph);
          var w = 1500//10 * Object.keys(graph["links"]).length;
          var h = 1200
          
          var linksToRemove = [];
          var linksAux = [];
          var svg = d3.select("body").append("svg:svg")
          .attr("width", w)
          .attr("height", h);
          
          svg.append("defs").append("marker")
          .attr({
          "id": "arrowhead",
          "viewBox": "-0 -5 10 10",
          "refX": 19,
          "refY": 0,
          "orient": "auto",
          "markerWidth": 13,
          "markerHeight": 13,
          "xoverflow": "visible"
          })
          .append("svg:path")
          .attr("d", "M 0,-5 L 10 ,0 L 0,5")
          .attr("fill", "#999")
          .style("stroke", "none");
          
          //Positioning genes on top and opening all
          var countGene = 0,
          countGeneLine = 0,
          countPValueAbove = 0,
          altAbove = 0;
          var level = 0,countNodesOnLevel=[];
          for (var i = 0; i < graph.nodes.length; i++) countNodesOnLevel[i] = 0;
          
          var link, node, node_drag, linkNodes = false;
          var div = d3.select("body").append("div")   
          .attr("class", "tooltip")
          .style("opacity", 0);
          init();
          
          function init() {
          
          //Dísplay nodes on the positions
          graph.nodes.forEach(function(el) {
          /*			if (el.type === "gene") {
          if(countGene === 12){
          countGene = 0;
          countGeneLine ++;
          }
          
          el.x = 50 + (90 * countGene);
          el.y = 100 + (40 * countGeneLine);
          el.fixed = true;
          el.opened = true;
          var aux = [];
          displayTreeRecursively(el);
          countGene++;
          } else { */
          if (el.pValue > pValue) {
          countPValueAbove++;
          altAbove++;
          if (altAbove === 4) altAbove = 0;
          el.x = 50 + 30 * countPValueAbove;
          el.y = 200 + (40 * countGeneLine) + 50 * altAbove;
          el.fixed = true;
          el.opened = false;
          }
          //	}// if (el.type === "gene")
          });//graph.nodes.forEach(function(el)
          
          var force = self.force = d3.layout.force()
          .nodes(graph.nodes)
          .links(graph.links)
          .gravity(1)
          .distance(100)
          .charge(-10000)
          .size([w, h])
          .start();
          
          link = svg.selectAll("line.link")
          .data(graph.links)
          .enter().append("svg:line")
          .attr("class", "link")
          .attr("marker-end", "url(#arrowhead)")
          .attr("x1", function(d) {
          return d.source.x;
          })
          .attr("y1", function(d) {
          return d.source.y;
          })
          .attr("x2", function(d) {
          return d.target.x;
          })
          .attr("y2", function(d) {
          return d.target.y;
          });
          //color for lines
          /*.style("stroke",function(d){
          n = graph.nodes.filter(function(el) {
          return el.id === d.source.id;
          })[0];
          if (n.type.indexOf("gene") >= 0) return "#ffd900";
          else {
          if (d.pValue <= pValue) {
          return "#e96a37";
          } else
          return "#50a161";
          } 
          })*/
          

          
          node_drag = d3.behavior.drag()
          .on("dragstart", dragstart)
          .on("drag", dragmove)
          .on("dragend", dragend);
          
          //Creating nodes with the data loaded from json. Ading click actions and tooltip
          node = svg.selectAll("g.node")
          .data(graph.nodes)
          .enter().append("svg:g")
          .attr("class", "node")
          .on("click", click)
          .call(node_drag)
          .on("mouseover", function(d) {
          d3.selectAll(".link").attr("class", "link");
          d3.selectAll(".link")
          .filter(function(l) {
          return (l.source === d) || (l.target === d);
          })
          .attr("class", "link selected");
          
          div.transition()				
          .duration(200)		  
          .style("opacity", .9);
          
          var str = "<b>Name:</b> " + d.name +"<br/>"+
          d.longname + "<br/>";
          
          if(d.type == "gene"){
          str += "<b>DGE p-value:</b> " + d.pValue + "<br/><b>LogFC: </b>"+d.logFC; 
          }else{
          str += "<b>GSEA p-value:</b> " + d.pValue + "<br/><b>Type: </b>"+d.type; 
          }
          div .html(str)  
          .style("left", (d3.event.pageX) + "px")		 
          .style("top", (d3.event.pageY - 28) + "px");
          })	//end of .on("mouseover", function(d)							  
          .on("mouseout", function(d) {		
          d3.selectAll(".link").attr("class", "link");				
          div.transition()				
          .duration(500)		  
          .style("opacity", 0);   
          });
          
          
          //Add a circle to the node and define its color
          node.filter(function(d){
          return d.type !== "gene";})
          .append("circle")
          .attr("r", 14)
          .style("cursor", "pointer")
          .style("fill", function(d, i) {		 
          if(d.type==="GOCCID"){
          if(d.pValue>pValue) return "#ffe44d";
          else return "#ccae00"
          }else if(d.type==="GOBPID"){
          if(d.pValue>pValue) return "#ceb6d1";
          else return "#a377a9";
          }else if(d.type==="GOMFID"){
          if(d.pValue>pValue) return "#45ff9a";
          else return "#00c45a";
          }		
          });
          
          var dim = 32
          
          //Coloring Genes and adding rectangle
          node.filter(function(d){
          return d.type === "gene";})
          .append("rect")
          .attr("x", -dim)
          .attr("y", -dim)
          .attr("width", dim*2.5)
          .attr("height", dim)
          .style("cursor", "pointer")
          .style("fill", function(d, i) {
          return d.color;
          });
          
          
          //Signal that indicate if the node it"s opened or not
          node.filter(function(d){
          return d.type !== "gene";})
          .append("svg:text")
          .attr("class", "openCircle")
          .attr("dx", -4)
          .attr("dy", 5)
          .style("font-weight", "bold")
          .style("font-size", "16")
          .style("cursor", "pointer")
          .text(function(d) {
          if (d.opened) {
          return "-";
          } else {
          return "+";
          }
          });
          
          //Adding longName and name to the nodes
          node.append("svg:text")
          .attr("class", "nodetext")
          .attr("dx", function(d) {
          if(d.type === "gene"){
          return -dim +5 
          }else{
          return 15 
          }								
          })
          .attr("dy", function(d) {
          if(d.type === "gene"){
          return -(dim/2)-2 
          }else{
          return 5 
          }
          })
          .style("font-size", function(d) {
          if(d.type === "gene"){
          return 12
          }else{
          return 12
          }
          })
          .text(function(d) {
          if(d.type === "gene"){
          return d.name
          }else{
          return d.longname.substring(0,15) + "...";										   
          }
          });
          
          //Adding pValue into the nodes
          node.append("svg:text")//#AQUI DONE
          .attr("class", "nodetext")
          .attr("dx", function(d) {
          if(d.type === "gene"){
          return -dim +5 
          }else{
          return 15 
          }								})
          .attr("dy", function(d) {
          if(d.type === "gene"){
          return -5 
          }else{
          return 18 
          }								})
          .style("font-size", "10")
          .text(function(d) {
          return d.pValue
          });
          
          force.on("tick", tick);
          
          //Opening nodes, fixing nodes in the screen and positioning
          graph.nodes.forEach(function(el) {
          el.fixed = false;
          el.opened = false;
          
          /* if (el.type === "gene") {
          countGene++;
          el.x = 150 * countGene;
          el.y = 200 + (40 * countGeneLine);
          el.fixed = true;
          el.opened = true;
          } else {
          if (el.pValue <= pValue) {
          if (!el.x || !el.y) displayTreeRecursively(el);
          
          el.fixed = true;
          el.opened = true;
          } else {
          countPValueAbove++;
          altAbove++;
          
          el.fixed = true;
          el.opened = false;
          }
          }*/
          });//end of	graph.nodes.forEach(function(el)
          
          updateNodesLinks();
          
          }//end of function init()
          
          function displayTreeRecursively(node) {//AQUI review
          var linksFromThatNode = [];
          //Search for the links that are source from the node clicked
          graph.links.forEach(function(link) {
          if (node.id === link.source) {
          var n = graph.nodes.filter(function(n) {
          return n.id === link.target;
          })[0];
          if (n && n.pValue <= pValue) {
          link.target = n;
          linksFromThatNode.push(link);
          }
          }
          });//end of graph.links.forEach(function(link) 
          
          //Sort the array utilizing it"s pValue
          linksFromThatNode.sort(function(a, b) {
          return (-1) * (a.target.pValue - b.target.pValue);
          });
          
          //Find the node that the link target to
          linksFromThatNode.forEach(function(link) {
          var n = graph.nodes.filter(function(n) {
          return n.id === link.target.id;
          })[0];
          level++;
          if (!n.x || !n.y) {
          countNodesOnLevel[level]++;
          n.x = 300+ 100*countNodesOnLevel[level];
          n.y = 400+ 50*level;
          }
          
          
          displayTreeRecursively(n);
          level--;
          });//linksFromThatNode.forEach(function(link)
          }//function displayTreeRecursively(node)
          
          function dragstart(d, i) {
          force.stop() // stops the force auto positioning before you start dragging
          }
          
          function dragmove(d, i) {
          d.px += d3.event.dx;
          d.py += d3.event.dy;
          d.x += d3.event.dx;
          d.y += d3.event.dy;
          tick(); // this is the key to make it work together with updating both px,py,x,y on d !
          }
          
          function dragend(d, i) {
          d.fixed = true; // of course set the node to fixed so the force doesn"t include the node in its auto positioning stuff
          tick();
          force.resume();
          }
          
          function tick() {
          link.attr("x1", function(d) {
          return d.source.x;
          })
          .attr("y1", function(d) {
          return d.source.y;
          })
          .attr("x2", function(d) {
          return d.target.x;
          })
          .attr("y2", function(d) {
          return d.target.y;
          });
          
          node.attr("transform", function(d) {
          return "translate(" + d.x + "," + d.y + ")";
          });
          }
          
          //Called when there is a click on node
          function click(d) {
          if (d3.event.defaultPrevented) return; // ignore drag
          updateNodesLinks(d);
          }
          
          function update(nodes, links) {//#AQUI
          
          // Update links.
          link = link.data(links);
          //link = link.data(links, function(d) { return d.source.id; });
          
          link.exit().remove();
          
          link.enter().insert("line", ".node")
          .attr("class", "link")
          .attr("marker-end", "url(#arrowhead)");
          
          // Update nodes.
          // IMPORTANT - Force node to mantain its index id
          node = node.data(nodes, function(d) {
          return d.id;
          });
          
          node.exit().remove();
          
          var nodeEnter = node.enter().append("g")
          .attr("class", "node")
          .on("click", click)
          .call(node_drag)
          .on("mouseover", function(d) {
          d3.selectAll(".link").attr("class", "link")
          .filter(function(l) {
          return (l.source === d) || (l.target === d);
          })
          .attr("class", "link selected");
          
          div.transition()				
          .duration(200)		  
          .style("opacity", .9);
          var str = "<b>Name:</b> " + d.name +"<br/>" + d.longname + "<br/>" ;
          if(d.type == "gene"){
          str += "<b>DGE p-value:</b> " + d.pValue + "<br/><b>LogFC: </b>"+d.logFC; 
          }else{
          str += "<b>GSEA p-value:</b> " + d.pValue + "<br/><b>Type: </b>"+d.type; 
          }
          div .html(str)  
          .style("left", (d3.event.pageX) + "px")		 
          .style("top", (d3.event.pageY - 28) + "px");
          })//end of .on("mouseover", function(d) 								  
          .on("mouseout", function(d) {  
          d3.selectAll(".link").attr("class", "link");				
          div.transition()				
          .duration(500)		  
          .style("opacity", 0);   
          });//end of .on("mouseout", function(d)
          
          //Add a circle to the node and define its color
          nodeEnter.filter(function(d){
          return d.type !== "gene";})
          .append("circle")
          .attr("r", 14)
          .style("cursor", "pointer")
          .style("fill", function(d, i) {		 
          if(d.type==="GOCCID"){
          if(d.pValue>pValue) return "#ffe44d";
          else return "#ccae00"
          }else if(d.type==="GOBPID"){
          if(d.pValue>pValue) return "#ceb6d1";
          else return "#a377a9";
          }else if(d.type==="GOMFID"){
          if(d.pValue>pValue) return "#45ff9a";
          else return "#00c45a";
          }		
          });
          
          var dim = 32;
          
          //Coloring Genes and adding rectangle//#aqui
          nodeEnter.filter(function(d){
          return d.type === "gene";})
          .append("rect")
          .attr("x", -dim)
          .attr("y", -dim)
          .attr("width", dim*2.5)
          .attr("height", dim)
          .style("cursor", "pointer")
          .style("fill", function(d, i) {
          return d.color;
          });
          
          
          //Signal that indicate if the node it"s opened or not
          nodeEnter.filter(function(d){
          return d.type !== "gene";})
          .append("svg:text")
          .attr("class", "openCircle")
          .attr("dx", -4)
          .attr("dy", 5)
          .style("font-size", "16")
          .style("font-weight", "bold")
          .style("cursor", "pointer")
          .text(function(d) {
          if (d.opened) {
          return "";
          } else {
          return "+";
          }
          });
          nodeEnter.append("svg:text")
          .attr("class", "nodetext")
          .attr("dx", function(d) {
          if(d.type === "gene"){
          return -dim +5 
          }else{
          return 15 
          }								
          })
          .attr("dy", function(d) {
          if(d.type === "gene"){
          return -(dim/2)-2 
          }else{
          return 5 
          }								})
          .style("font-size", function(d) {
          if(d.type === "gene"){
          return 12
          }else{
          return 12
          }
          })
          .text(function(d) {
          if(d.type === "gene"){
          return d.name
          }else{
          return d.longname.substring(0,15) + "...";										   
          }
          });
          nodeEnter.append("svg:text")
          .attr("class", "nodetext")
          .attr("dx", function(d) {
          if(d.type === "gene"){
          return -dim +5 
          }else{
          return 15 
          }								})
          .attr("dy", function(d) {
          if(d.type === "gene"){
          return -5 
          }else{
          return 18 
          }								})
          .style("font-size", "10")
          .text(function(d) {
          return d.pValue
          });
          
          
          //Remove all the signs to open node, to later add them only on the ones realy opened
          node.selectAll(".openCircle").remove();
          
          //Signal that indicate if the node it"s opened or not
          node.filter(function(d){
          return d.type !== "gene";})
          .append("svg:text")
          .attr("class", "openCircle")
          .attr("dx", -4)
          .attr("dy", 5)
          .style("cursor", "pointer")
          .style("font-weight", "bold")
          .style("font-size", "16")
          .text(function(d) {
          if (d.opened) {
          return "";
          } else {
          return "+";
          }
          });
          //Color for line
          /*link.style("stroke",function(d){
          n = nodes.filter(function(el) {
          return el.id === d.source.id;
          })[0];
          if (n.type.indexOf("gene") >= 0) return "#ffd900";
          else {
          if (d.pValue <= pValue) {
          return "#2daec6";
          } else
          return "#50a161";
          }
          })*/
          tick();
          }//end of function update(nodes, links)
          
          function updateNodesLinks(d) {
          var i, j;
          var node = [];
          var link = [];
          var linksForThatNode = [];
          
          //Defining if the node are suposed to open/close (toggle)
          if (d) {
          graph.nodes.forEach(function(el) {
          if (d.name === el.name) {
          if (el.opened) {
          el.opened = 0;
          closeAndOpenNodesRecursively(d, false);
          } else {
          el.opened = 1;
          //closeAndOpenNodesRecursively(d,true);
          }
          }
          });
          }
          
          //Put links to a auxiliary array
          linksAux = [];
          graph.links.forEach(function(l) {
          linksAux.push(l);
          });
          linksToRemove = [];
          //Iterate through nodes to close one by one
          graph.nodes.forEach(function(n) {// #AQUI
          if (!n.opened) {
          //Close node and return links to be removed from the link list
          var linkFinal = removeNodesLinks(n);
          linksToRemove = linksToRemove.concat(linkFinal);
          //Remove links from list
          linksToRemove.forEach(function(l) {
          var index = linksAux.indexOf(l)
          if (index >= 0) {
          linksAux.splice(index, 1);
          }
          });
          }
          });
          
          //Final array with the links maintained
          link = link.concat(linksAux)
          
          
          //Managing the nodes that are suposed to be Added according to the links
          graph.nodes.forEach(function(n) {
          
          //If it isn"t gene, verify if it"s a node without links - Don"t close gene nodes
          if (document.getElementById("filterGeneNodes").checked || n.type !== "gene") {
          if (n.pValue > pValue || n.type === "gene") {
          var addIt = false;
          link.forEach(function(l) {
          
          if (n.id === l.source.id || n.id === l.target.id) {
          addIt = true;
          }
          });
          if (addIt) {
          node.push(n);
          }
          } else {
          node.push(n);
          }
          } else {
          node.push(n);
          }
          
          });
          
          //Add links between genes and closed nodes that are visible. 
          var linksNodesToGenes = [];
          var linksNodesToGenesFound = []
          node.forEach(function(n) {
          linksNodesToGenes = graph.links.filter(function(l) {
          //if(n.name === l.target.name && l.source.type === "gene" &&
          //(l.target.type === "GOBPID" || l.target.type === "GOCCID" || l.target.type === "GOMFID") ){ #AQUI REMOVED FOR ALLOWING LINKS TO ANY GO DISPLAYED//n.children === null
          if((n.name === l.target.name) && 
          (l.target.type === "GOBPID" || l.target.type === "GOCCID" || l.target.type === "GOMFID") ){
          return true;
          }
          return false;
          });
          linksNodesToGenes.forEach(function(link1) {
          var found = false;
          link.forEach(function(link2) {
          if(link1.source.name == link2.source.name && link1.target.name == link2.target.name){
          found = true;
          }
          });
          if(!found){
          linksNodesToGenesFound.push(link1);
          }
          });
          });
          
          //Add just links that have a visible nodes 
          for(var x = linksNodesToGenesFound.length-1; x >= 0;x--){
          for(t = 0; t < node.length;t++){
          if(node[t].name === linksNodesToGenesFound[x].source.name){
          link.push(linksNodesToGenesFound[x]);
          break;
          }
          }
          }
          
          
          //Update nodes and links in the d3 lib
          update(node, link);
          d3.selectAll(".link").attr("class", "link");	
          }//end of function updateNodesLinks(d)
          
          function closeAndOpenNodesRecursively(node, open) {//#AQUI
          var linksFromThatNode = [];
          //Search for the links that target to the node clicked
          linksFromThatNode = graph.links.filter(function(link) {
          if(node.name === link.target.name){
          return true;
          }
          return false;
          });
          //Sort the array utilizing it"s pValue
          linksFromThatNode.sort(function(a, b) {
          return (a.target.pValue - b.target.pValue);
          });
          //Find the node that the link is source from
          linksFromThatNode.forEach(function(link) {
          var n = graph.nodes.filter(function(n) {
          return n.id === link.source.id;
          })[0];
          n.opened = open;
          closeAndOpenNodesRecursively(n, open);
          });
          }//end of function closeAndOpenNodesRecursively(node, open)
          
          function removeNodesLinks(node) {
          var linksTOThatNode = [];
          //Search for the links that are target from the node clicked.
          //If the connection is between genes and nodes above pValue, the connection can be removed
          linksTOThatNode = linksAux.filter(function(link) {
          return (node.name === link.target.name && link.source.type !== "gene") || ((node.name === link.target.name) && node.pValue > pValue);
          });
          //Sort the array utilizing it"s pValue
          linksTOThatNode.sort(function(a, b) {
          return (a.source.pValue - b.source.pValue);
          });
          //Find the node that the link target to
          linksTOThatNode.forEach(function(linkToN) {
          var n = graph.nodes.filter(function(n) {
          return n.id === linkToN.source.id;
          })[0];
          //Verify if there is connection source from this node
          if (!verifyConnection(node, n)) {
          var removeIt = false;
          //Verify if that node has links to be removed
          linksAux.forEach(function(l) {
          //Verify if there"s a source from that node
          if (n.id === l.source.id) {
          removeIt = true;
          }
          });
          if (removeIt) {
          //Navigate to the target node recursively to search for more links to be removed
          linksToRemove = linksToRemove.concat(removeNodesLinks(n));
          }
          }
          });
          
          return linksTOThatNode;
          }//function removeNodesLinks(node)
          
          function verifyConnection(nodeSource, nodeVerified) {
          var connected = false;
          var arrayAux = [];
          arrayAux = arrayAux.concat(graph.links);
          //Remove links from the original list to get an updated list
          linksToRemove.forEach(function(l) {
          var index = arrayAux.indexOf(l)
          if (index >= 0) {
          arrayAux.splice(index, 1);
          }
          });
          arrayAux.forEach(function(l) {
          //Verify if there"s a target to that node differente from the source #AQUI!!
          if ((nodeVerified.id === l.target.id) || (nodeVerified.id === l.source.id)) {
          
          connected = true;
          }
          });
          return connected;
          }
          
          function linkUnlink() {
          linkNodes = !linkNodes;
          if (linkNodes) {
          graph.nodes.forEach(function(el) {
          el.opened = true;
          });
          } else {
          graph.nodes.forEach(function(el) {
          el.opened = false;
          });
          
          }
          updateNodesLinks();
          }
          
          function setPValue(){
          pValue= parseFloat(document.getElementById("pValue").value);
          document.getElementById("L.Bpv").innerHTML="> " + pValue;
          document.getElementById("L.Spv").innerHTML="<= " + pValue;
          level = 0,countNodesOnLevel=[];
          for (var i = 0; i < graph.nodes.length; i++) countNodesOnLevel[i] = 0;
          node = svg.selectAll("g.node").remove();
          link = svg.selectAll("line.link").remove();
          
          init();
          }
          
          function activateButton(b){
              b.disabled = false
              b.style.color="black"
          }
          
          activateButton(document.getElementById("graphunlink"));
          document.getElementById("graphunlink").onclick = linkUnlink;

          activateButton(document.getElementById("graphstatic"));
          //document.getElementById("graphstatic").addEventListener("click", freeze);
          document.getElementById("graphstatic").onclick = freeze;
          

          activateButton(document.getElementById("graphfilterp"));
          document.getElementById("graphfilterp").onclick = setPValue;
          
          activateButton(document.getElementById("filterGeneNodes"));
          document.getElementById("filterGeneNodes").onchange = updateNodesLinks;          
          });  
          </script>
          </body>
          
          
          ')#HTML
          )
          )
          )#wellPanel      
          )

server <- function(input, output, session) {
  ###############
  ##### Var #####
  ###############
  
  # Species selected on Step 1
  sp <- reactive({
    which(genomes$Species == input$species)
  })
  
  # Table uploaded on Step 2 #AQUI
  data1 <- reactive({
   
    inFile <- input$geneSet.file
    
    if (is.null(inFile))
      return(NULL)
    
    as.data.frame(read.csv(inFile$datapath, header=input$header, sep=input$sep,
      quote=input$quote))
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
      "symbol" = select(get(genomes$DB[sp()]), keys=as.character(gene.list.universe()),
        columns= c("ENTREZID", "SYMBOL", "GENENAME"),
        keytype="SYMBOL"),
      "entrez" = select(get(genomes$DB[sp()]), keys=as.character(gene.list.universe()),
        columns= c("ENTREZID", "SYMBOL", "GENENAME"),
        keytype="ENTREZID"),
      "ensembl" = select(get(genomes$DB[sp()]), keys=as.character(gene.list.universe()),
        columns= c("ENTREZID", "SYMBOL", "GENENAME","ENSEMBL"),
        keytype="ENSEMBL"),
      "pubmed" = select(get(genomes$DB[sp()]), keys=as.character(gene.list.universe()),
        columns= c("ENTREZID", "SYMBOL", "GENENAME","PMID"),
        keytype="PMID"),
      "refseq" = select(get(genomes$DB[sp()]), keys=as.character(gene.list.universe()),
        columns= c("ENTREZID", "SYMBOL", "GENENAME","REFSEQ"),
        keytype="REFSEQ"),
      "unigene" = select(get(genomes$DB[sp()]), keys=as.character(gene.list.universe()),
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
  
  output$colnames.gene <- renderUI({ #AQUI
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
          input$gene.id.universe,inline = TRUE)
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
    if (sp()==1){
      paste("Please, select a species")
    }else{
      paste("OK")
    }
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
                        The genes' 'universe' contains only", length(lists.intersect), "out of the", length(set[-na.set]),
                        'genes identified on step 2. This means that the following genes will not be used in the GSEA:',
                        setdiff(set, lists.intersect),". You may consider choosing a different source."))
                    }else{
                      return("None of the genes listed on Step 2 were found in the 'Gene Universe' provided on Step 3. Please, change the Step 3 parameters or choose another source.")
                      # return("Simplicity could not find annotation for the genes, please change the parameters or choose another source.")
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
        if(gene.ann() == ""){
          paste("Please, complete all steps to enable GSEA.")
        }else{
          actionButton("start.gsea", "Start GSEA")
        }}, error = function(e){
          paste("Please, complete all steps to enable GSEA.")
        }
      )
    }else{
      paste("Complete all steps to enable GSEA.")
    }}, error = function(e){
      paste("Complete all steps to enable GSEA.")
      }
    )
    
  })
  
  ##################
  ##### STEP 5 #####
  ##################
  
  observeEvent(input$start.gsea, {
    output$enable.gsea <- renderUI({paste("Your analysis has started!")})
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
    
    ####### Progress Bar #########
    # Extracted from https://shiny.rstudio.com/articles/progress.html
    
    # Create a Progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Computing hypergeometric tests", value = 0)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    
    # Create a callback function to update progress.
    # Each time this is called:
    # - If `value` is NULL, it will move the progress bar 1/3 of the remaining
    #   distance. If non-NULL, it will set the progress to that value.
    # - It also accepts optional detail text.
    updateProgress <- function(value = NULL, detail = NULL) {
      if (is.null(value)) {
        value <- progress$getValue()
        value <- value + (progress$getMax() - value) / 3
      }
      progress$set(value = value, detail = detail)
    }
    
    ##############################
    
    results <- simplicity_gostats(Gene.Ann = Gene.Ann, SP = SP, Universe.Sets = Universe.Sets,
      Universe.Source = Universe.Source, Universe.Gene.Ann =  Universe.Gene.Ann, pvalue = 0.05, updateProgress =  updateProgress, data1 = data1())
    if(is.null(results)){
      output$enable.gsea <- renderUI({paste("The analysis didn't find any enriched GO.")})
    }else{
      
    # Plotting for Graph Step 6
    observe({
      session$sendCustomMessage(type="SendObjectToClientDynamicCallbackHandler", as.character(toJSON(results$graph)))#     
    })#observe
    
    
    #Tabs
    result<-results$gsea
    
    output$resultsInTabs <- renderUI({
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
    
    lapply(names(result), function(f) {
      output[[f]] <- renderDataTable(as.data.frame(summary(result[[f]],"htmlLinks" = TRUE)), escape = FALSE)
    })
    output$enable.gsea <- renderUI({paste("Your analysis is finished. Scroll down to see the results.")})
    }#if(is.null(results))
  })#obserEvent
  
  ##################
  ##### STEP 6 #####
  ##################
  
  }#server

shinyApp(ui = ui, server = server)
