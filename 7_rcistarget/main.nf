#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process rcistarget {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    path f

  when:
    ( ! file("${params.project_folder}/deseq2_output/annotated/${f}".replace(".results.tsv",".RcisTarget.xlsx")).exists() ) 
  
  script:
  """
    #!/usr/bin/Rscript
    library(RcisTarget)
    library(openxlsx)
    setwd("/workdir/deseq2_output/annotated/")
    fout=stringr::str_split("${f}", ".results.tsv")[[1]][[1]]
    fout=paste0(fout,".RcisTarget.xlsx")
    # Load gene sets to analyze. e.g.:
    deg = read.delim("${f}", as.is = TRUE)
    sigGene = subset(deg, padj <= 0.05)[, 'gene_name']
    sigGeneLists <- list(geneListName=sigGene)
    motif_files = list(chip = c("/workdir/${params.rcis_db}/chip_collection.feather", 
                                "workdir/${params.rcis_db}/chip_annotation.tsv"),
                      tf   = c("/workdir/${params.rcis_db}/tf_collection.feather", 
                                "/workdir/${params.rcis_db}/tf_annotation.tsv"))
    results = list()
    for(db in names(motif_files)){
        if(file.exists(motif_files[[db]][1])){
            ## 0. load motif data
            motifRankings <- importRankings(motif_files[[db]][1])
            motif_anno = read.delim(motif_files[[db]][2], as.is = TRUE)
      
            ## 1. Calculate AUC
            motifs_AUC <- calcAUC(sigGeneLists, motifRankings)
            
            ## 2. Select significant motifs, add TF annotation & format as table
            motifEnrichmentTable <- addMotifAnnotation(motifs_AUC, nesThreshold=3)
            
            ## 3. Identify significant genes for each motif
            motifEnrichmentTable_wGenes <- addSignificantGenes(motifEnrichmentTable, 
                                                              geneSets=sigGeneLists,
                                                              rankings=motifRankings)
            ## 4. Add motif name
            motifEnrichmentTable_wGenes = merge(motif_anno, motifEnrichmentTable_wGenes, by.x = 'X.motif_id', by.y = 'motif', all.y = TRUE)
      
            ## 5. Store data
            results[[db]] = motifEnrichmentTable_wGenes[order(motifEnrichmentTable_wGenes[, 'NES'], decreasing = TRUE), ]
      
            # make flattened table
            toFlat = list()
            for(i in 1:nrow(motifEnrichmentTable_wGenes)){
                genes_column = unlist(strsplit(motifEnrichmentTable_wGenes[i,'enrichedGenes'], ';'))
                toFlat[[i]] = data.frame(motifID = motifEnrichmentTable_wGenes[i,'X.motif_id'],
                                        motifName = motifEnrichmentTable_wGenes[i,'gene_name'],
                                        targetGene = genes_column,
                                        NES = motifEnrichmentTable_wGenes[i,'NES'])
            }
            results[[paste0('flat_', db)]] = do.call(rbind, toFlat)
        }
    }
    write.xlsx(results, fout, row.names = FALSE)
  """
}

workflow {
  if ( 'rcis_db' in params.keySet() ) {
    data = channel.fromPath( "${params.project_folder}/deseq2_output/annotated/*.results.tsv" )
    rcistarget(data) 
} 