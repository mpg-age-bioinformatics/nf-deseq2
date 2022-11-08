#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process cellplot {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    path inFile
    val filetype
    val category
    val nterms

  when:
    ( ! file("${params.project_folder}/deseq2_output/annotated/${inFile}".replace(".tsv",".${category}.cellplot.pdf")).exists() ) 
  
  script:
  """
    #!/usr/bin/Rscript

    # USAGE: `Rscript david_to_cellplot.R /beegfs/group_bit/data/projects/departments/Thomas_Langer/TL_Kai_RNAseq_5lines/adiff_output/Ka5R/wt_gfp.vs.cko_nlrp.DAVID.tsv tsv KEGG_PATHWAY 10`
    # works with csv, tsv, txt and xlsx

    # need to have CellPlot installed: devtools::install_github("dieterich-lab/CellPlot", build_vignettes = TRUE)
    # and openxlsx: install.library('readxl')

    rm(list = ls())

    args<-commandArgs(TRUE)


    #if (!require(devtools)) install.packages('devtools')
    library(devtools)

    #if (!require(readxl)) install.packages('readxl')
    library(readxl)

    #if (!require(CellPlot)) devtools::install_github("dieterich-lab/CellPlot")
    library(CellPlot)
    
    setwd("/workdir/deseq2_output/annotated/")

    # read in data, either excel or tsv
    # reformat input data
    inFile <- toString("$inFile")
    filetype <- toString("$filetype")
    category <- toString("$category")
    nterms <- as.numeric("$nterms")

    filetype_map <- c("xlsx" = 'xlsx',  'tsv' = '\t', 'csv' = ',', 'txt'=" ")

    if(filetype == 'xlsx'){
        D <- read_excel(inFile, sheet = category)
        D <- as.data.frame(D)
      } else {
        D <- read.csv(inFile, header = TRUE, sep = filetype_map[filetype], as.is = TRUE)
    }

    D<-D[D["categoryName"] == category, ]
    if ( nrow(D) == 0 ) {
      quit(save="no")
    }

    D\$ease <- as.numeric(as.character(D\$ease))
    D\$foldEnrichment <- as.numeric(as.character(D\$foldEnrichment))
    D\$listHits <- as.numeric(as.character(D\$listHits))
    D <- D[order(D\$ease),]

    # subset to number of rows to plot..if specified number is larger than number of rows.
    if (nterms >= nrow(D)) nterms = nrow(D)
    D <- D[1:nterms,]

    # log2FoldChange as list
    D\$log2fc <-lapply( gsub('inf', 'Inf', D\$log2fc), function(x) as.numeric(as.character(unlist(strsplit(toString(x), ", ")))))

    # cellplot    
    x <- D

    pdf(gsub(filetype, paste(category, '.cellplot.pdf', sep = ''), inFile))
    cell.plot(x = setNames(-log10(D\$ease), D\$termName), 
                  cells = D\$log2fc, 
                  main ="GO enrichment",
                  xlab ="-log10(P.Value)", 
                  x.mar = c(max(unlist(lapply(D\$termName, function(x) nchar(x))))/100 + 0.1,0),
                  key.n = 7, 
                  y.mar = c(0.1, 0.1), 
                  cex = 1.6, 
                  cell.outer = 3, 
                  bar.scale = .7, 
                  space = .2)

    dev.off()

    # symplot
    pdf(gsub(filetype, paste(category, '.symplot.pdf', sep = ''), inFile))
    sym.plot(x = setNames(-log10(D\$ease), D\$termName), 
                cells = D\$log2fc, 
                x.annotated = D\$listHits, 
                main = "GO enrichment",
                key.lab = "-log10(P.Value)",
                x.mar = c(max(unlist(lapply(D\$termName, function(x) nchar(x))))/100 + 0.1, 0),
                y.mar = c(0.2,0.1),
                key.n = 7, 
                cex = 1.6, 
                axis.cex = .8, 
                group.cex = .7) 

    dev.off()
   
  """
}


workflow david_bp {
  DAVID = channel.fromPath( "${params.project_folder}/deseq2_output/annotated/*.DAVID.tsv" )
  cellplot(DAVID,"tsv", "GOTERM_BP_FAT", "15") 
} 

workflow david_kegg {
  DAVID = channel.fromPath( "${params.project_folder}/deseq2_output/annotated/*.DAVID.tsv" )
  cellplot(DAVID,"tsv", "KEGG_PATHWAY", "15") 
}

workflow topgo_bp {
  topGO = channel.fromPath( "${params.project_folder}/deseq2_output/annotated/*.topGO.tsv" )
  cellplot(topGO,"tsv", "GOTERM_BP_FAT", "15") 
}

workflow topgo_kegg {
  topGO = channel.fromPath( "${params.project_folder}/deseq2_output/annotated/*.topGO.tsv" )
  cellplot(topGO,"tsv", "KEGG_PATHWAY", "15")
}

workflow {
  david_bp()
  david_kegg()
  topgo_bp()
  topgo_kegg()
}

