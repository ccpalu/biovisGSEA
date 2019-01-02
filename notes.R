install.packages("BiocManager")
BiocManager::valid()
BiocManager::install(c(
"ade4", "annotate", "AnnotationDbi", "AnnotationForge",
"AnnotationHub", "Biobase", "BiocGenerics", "BiocInstaller",
"Category", "cli", "dendextend", "digest", "dplyr", "edgeR",
"evaluate", "fansi", "genefilter", "ggplot2", "GO.db",
"GOstats", "graph", "GSEABase", "interactiveDisplayBase",
"IRanges", "later", "limma", "mime", "NMF", "org.Ag.eg.db",
"org.At.tair.db", "org.Bt.eg.db", "org.Ce.eg.db",
"org.Cf.eg.db", "org.Dm.eg.db", "org.Dr.eg.db",
"org.EcK12.eg.db", "org.EcSakai.eg.db", "org.Gg.eg.db",
"org.Hs.eg.db", "org.Mm.eg.db", "org.Mmu.eg.db",
"org.Pf.plasmo.db", "org.Pt.eg.db", "org.Rn.eg.db",
"org.Sc.sgd.db", "org.Ss.eg.db", "org.Xl.eg.db", "pkgconfig",
"R6", "RBGL", "Rcpp", "Rgraphviz", "rlang", "rstudioapi",
"S4Vectors", "scales", "stringi", "tidyr", "tidyselect",
"tinytex", "xfun", "XML", "xtable"
), update = TRUE, ask = FALSE)
library(rsconnect)
deployApp()
install.packages("digest")
install.packages("Rcpp")
install.packages(c("Rgraphviz", "XML", "genefilter", "later", "mime", "rlang"))
install.packages("RBGL")
install.packages("annotation")
install.packages("annotate")
install.packages("GSEABase")
install.packages("AnnotationForge")
install.packages("interactiveDisplayBase")


install.packages("IRanges")
install.packages("Biobase")
install.packages("AnnotationDbi")
install.packages("AnnotationForge")


###################################################################################################
# OK. We found annotation for 2173 genes and there are 188 not recognised. The gene 'universe' contains only 2 out of the 9 genes identified on step 2. This means that the following gnes will not be used in the GSEA: 14609 . You may consider choosing a different source.,OK. We found annotation for 2173 genes and there are 188 not recognised. The gene 'universe' contains only 2 out of the 9 genes identified on step 2. This means that the following gnes will not be used in the GSEA: 18619 . You may consider choosing a different source.,OK. We found annotation for 2173 genes and there are 188 not recognised. The gene 'universe' contains only 2 out of the 9 genes identified on step 2. This means that the following gnes will not be used in the GSEA: 20519 . You may consider choosing a different source.,OK. We found annotation for 2173 genes and there are 188 not recognised. The gene 'universe' contains only 2 out of the 9 genes identified on step 2. This means that the following gnes will not be used in the GSEA: 22152 . You may consider choosing a different source.,OK. We found annotation for 2173 genes and there are 188 not recognised. The gene 'universe' contains only 2 out of the 9 genes identified on step 2. This means that the following gnes will not be used in the GSEA: NA . You may consider choosing a different source.,OK. We found annotation for 2173 genes and there are 188 not recognised. The gene 'universe' contains only 2 out of the 9 genes identified on step 2. This means that the following gnes will not be used in the GSEA: 54215 . You may consider choosing a different source.,OK. We found annotation for 2173 genes and there are 188 not recognised. The gene 'universe' contains only 2 out of the 9 genes identified on step 2. This means that the following gnes will not be used in the GSEA: 668421 . You may consider choosing a different source.,OK. We found annotation for 2173 genes and there are 188 not recognised. The gene 'universe' contains only 2 out of the 9 genes identified on step 2. This means that the following gnes will not be used in the GSEA: 72719 . You may consider choosing a different source.

##########################################################################################
NEXT:
https://stackoverflow.com/questions/17502661/shiny-open-new-browser-tab-from-within-shiny-app

##########################################################################################
OTHER

THe colapsed version of the tree is ignoring the relationships between the GO nodes.
