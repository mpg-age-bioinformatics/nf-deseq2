#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process topgo {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    path f

  when:
    ( ! file("${params.project_folder}/deseq2_output/annotated/${f.baseName}.tsv".replace(".results.tsv",".topGO.tsv")).exists() ) 
  
  script:
  """
    #!/usr/bin/Rscript
    library(topGO)
    library(biomaRt)
    library(plyr)
    library(openxlsx)
    ## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    fin="/workdir/deseq2_output/annotated/${f.baseName}.tsv"
    Din = read.delim(fin, as.is = TRUE)
    # head(Din)
    # make geneList input for topGO
    sigGenes = subset(Din, padj <= 0.05)[, 'ensembl_gene_id']
    geneList = factor(as.integer(Din[,'ensembl_gene_id'] %in% sigGenes))
    names(geneList) = Din[,'ensembl_gene_id']
    results = list(genes = subset(Din, padj <= 0.05))
    ## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # This part if you are working in R version 4 or higher
    # biomartCacheClear()
    # ensembl <- useEnsembl(biomart = "genes", dataset = "${params.biomart_dataset}", host="${params.biomart_host}")
    host=stringr::str_split("${params.biomart_host}", "/biomart")[[1]][[1]]
    ensembl = useMart("ensembl",dataset="${params.biomart_dataset}", host=host)
    goterms = getBM(attributes = c("ensembl_gene_id", 
                                            "external_gene_name",
                                            "go_id", "name_1006"),  
                              filters = 'ensembl_gene_id',
                              values = Din[, 'ensembl_gene_id'],
                              mart = ensembl)
    ## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    geneID2GO = list()
    for(i in 1:nrow(goterms)){
      if(goterms[i, 'ensembl_gene_id'] %in% names(geneID2GO)){
        geneID2GO[[goterms[i, 'ensembl_gene_id']]] = c(geneID2GO[[goterms[i, 'ensembl_gene_id']]], 
                                                                goterms[i, 'go_id'])
      } else {
        geneID2GO[[goterms[i, 'ensembl_gene_id']]] = c(goterms[i, 'go_id'])
      }
    }
    GOTERM_categories = list('BP' = 'GOTERM_BP_FAT',
                            'MF' = 'GOTERM_MF_FAT', 
                            'CC' = 'GOTERM_CC_FAT')
    for(go_cat in names(GOTERM_categories)){
      ## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
      sampleGOdata <- new("topGOdata",
                          description = "Simple session",
                          ontology = go_cat,
                          allGenes = geneList,
                          nodeSize = 10,
                          annot = annFUN.gene2GO, 
                          gene2GO = geneID2GO)
      
      
      ## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
      resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
      
      
      ## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
      allRes = GenTable(sampleGOdata, classicFisher = resultFisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = sum(score(resultFisher) <= 0.1))
      names(allRes) <- c("GO.ID", "Term", "Annotated", "Significant", "Expected", "classicFisher")
      ## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
      gg = genesInTerm(sampleGOdata)
      ggsig = lapply(gg, function(x) x[x %in% sigGenes])
      geneIds = lapply(ggsig, function(x) paste(x, collapse = ', '))
      gn = lapply(ggsig,  function(x) paste(mapvalues(x, from = Din[, 'ensembl_gene_id'], to = Din[, 'gene_name'], warn_missing = FALSE), collapse = ', '))
      gfc = lapply(ggsig, function(x) paste(mapvalues(x, from = Din[, 'ensembl_gene_id'], to = Din[, 'log2FoldChange'], warn_missing = FALSE), collapse = ', '))
      
      D_gogenes = data.frame(goid = names(gg), geneIds = unlist(geneIds), gn = unlist(gn), gfc = unlist(gfc))
      
      
      ## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
      D_out = merge(allRes, D_gogenes, by.x = 'GO.ID', by.y = 'goid', all.x = TRUE, sort = FALSE)
      D_out[, 'categoryName'] = GOTERM_categories[[go_cat]]
      D_out[, 'termName'] = paste(D_out[,'GO.ID'], D_out[, 'Term'], sep = '~')
      D_out[,'percent'] = NA
      D_out[, 'listTotals'] = length(sigGenes)
      D_out[, 'listTotals_used'] = length(sigGenes(sampleGOdata))
      D_out[, 'popTotals'] = length(allGenes(sampleGOdata))
      D_out[, 'popTotals_used'] = numGenes(sampleGOdata)
      D_out[, 'foldEnrichment'] = (D_out[,'Significant']/D_out[,'listTotals_used']) / (D_out[, 'Annotated']/D_out[, 'popTotals_used'])
      D_out[, 'bonferroni'] = p.adjust(D_out[,'classicFisher'], method = "bonferroni")
      D_out[, 'benjamini'] = p.adjust(D_out[,'classicFisher'], method = "BY")
      D_out[, 'afdr'] = p.adjust(D_out[,'classicFisher'], method = "fdr")
      
      
      ## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
      D_out = D_out[, c('categoryName', 'termName', 'Significant', 'percent', 'classicFisher', 'geneIds',
                        'listTotals', 'listTotals_used', 'Annotated', 'popTotals', 'popTotals_used', 
                        'Expected',  'foldEnrichment', 'bonferroni', 'benjamini', 'afdr', "gn", "gfc")]
      
      names(D_out) <- c('categoryName', 'termName', 'listHits', 'percent', 'ease', 'geneIds',
                        'listTotals', 'listTotals_used', 'popHits', 'popTotals', 'popTotals_used',
                        'Expected', 'foldEnrichment', 'bonferroni', 'benjamini', 'afdr', 'genes name', 'log2fc')
      
      results[[GOTERM_categories[[go_cat]]]] = D_out
      
    }
    # save list as excel workbook
    write.xlsx(results, gsub("results.tsv","topGO.xlsx", fin), row.names = FALSE)
    # rbind and save as tsv
    results.tab = do.call(rbind, results[2:4])
    write.table(results.tab, gsub("results.tsv","topGO.tsv", fin), row.names = FALSE, sep = '\\t', quote = FALSE)
  """
}

workflow {
  data = channel.fromPath( "${params.project_folder}/deseq2_output/annotated/*.results.tsv" )
  topgo(data) 
} 