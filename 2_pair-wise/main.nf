#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process tx2gene {
  stageInMode 'symlink'
  stageOutMode 'move'

  output:
    val "/workdir/deseq2_output/deseq2.part1.Rdata", emit: rdata

  // when:
  //   ( ! file("${params.project_folder}/deseq2_output/deseq2.part1.Rdata").exists() ) 

  script:
  """
    #!/usr/bin/Rscript

    library(tidyverse)
    library(tximportData)
    library(tximport)
    library(DESeq2)
    library(rhdf5)
    library(readr)
    library(apeglm)
    tx2gene <- read_csv("/workdir/deseq2_output/tx2gene.csv")
    tx2gene <- tx2gene[rowSums(is.na(tx2gene)) == 0,]
    tx2gene <- droplevels(tx2gene)
    save.image("/workdir/deseq2_output/deseq2.part1.Rdata")
    sessionInfo()
  """
}

process deseq2 {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val input_file
    val rdata

  output:
    val input_file

  when:
    ( ! file("${params.project_folder}/deseq2_output/${input_file}".replace(".input.txt",".results.txt")).exists() ) 

  script:
  """
    #!/usr/bin/Rscript
    load("/workdir/deseq2_output/deseq2.part1.Rdata")
    library(tximportData)
    library(tximport)
    library(DESeq2)
    library(rhdf5)
    library(readr)
    library(apeglm)

    print("${input_file.baseName}.tsv")
    
    ref=stringr::str_split("${input_file.baseName}.tsv", "_vs_")[[1]][[2]]
    ref=stringr::str_split(ref, ".input.tsv")[[1]][[1]]

    coef=stringr::str_split("${input_file.baseName}.tsv", ".input.tsv")[[1]][[1]]

    out=paste("/workdir/deseq2_output/",coef,".results.tsv", sep="")

    # filein
    sampleTable<-read.delim2("/workdir/deseq2_output/${input_file.baseName}.tsv",sep = "\t", row.names = 1)
    samples<-row.names(sampleTable)

    # dir
    dir<-"/workdir/kallisto_output/"
    files<-file.path(dir, samples, "abundance.tsv")
    print(files)
    names(files)<-samples
    txi <- tximport(files, type = "kallisto", txOut = TRUE)
    txi <- summarizeToGene(txi, tx2gene = tx2gene)

    # model
    dds <- DESeqDataSetFromTximport(txi, sampleTable, ~group )

    # if circRNA file is present, add circRNA counts to dds, else proceed
    # here, get count table,
    # add circRNA counts
    circRNA_folder="${params.circRNA}"
    if(circRNA_folder != "None"){
      # gene count table
      gene_count = as.data.frame(counts(dds, normalized=FALSE))
      # circRNA table
      circRNA = read.delim('/workdir/${params.circRNA}/CircRNACount', check.names = FALSE, as.is = TRUE )
      coord = read.delim('/workdir/${params.circRNA}/CircCoordinates', check.names = FALSE, as.is = TRUE)
      # reformat circRNA header
      names(circRNA) <- gsub('.Chimeric.out.junction', '', names(circRNA))
      row.names(circRNA) <- paste0('circ_', coord[, 'Chr'], ':', coord[,'Start'], '-', coord[,'End'], '|', coord[,'Strand'], '|', coord[,'Gene'])
      circRNA = circRNA[,names(gene_count)]
      # add circRNAs to gene count table
      gene_count = rbind(gene_count, circRNA)
      # create DESeqDataSetFromMatrix
      dds <- DESeqDataSetFromMatrix(gene_count, sampleTable, ~group)
    }

    # model
    dds\$group <- relevel(dds\$group, ref = ref)
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]
    dds <- DESeq(dds)
    res.counts<-counts(dds, normalized=TRUE)

    # coef
    resLFC <- lfcShrink(dds,  coef=coef, type="apeglm")
    counts.dge<-merge(res.counts,resLFC, by=0, all=FALSE)
    row.names(counts.dge)<-counts.dge\$Row.names
    counts.dge<-counts.dge[, !(colnames(counts.dge) %in% c("Row.names"))]
    counts.dge <- counts.dge[order(counts.dge\$pvalue),]

    # file out
    write.table(counts.dge, out, sep="\\t")

    sessionInfo()
  """
}

process mastertable {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val input_files

  when:
    ( ! file("${params.project_folder}/deseq2_output/all_results_stats.xlsx").exists() ) 

  script:
  """
    #!/usr/bin/Rscript
    load("/workdir/deseq2_output/deseq2.part1.Rdata")
    library(tidyverse)
    library(tximportData)
    library(tximport)
    library(DESeq2)
    library(rhdf5)
    library(readr)
    library(apeglm)

    sampleTable<-read.delim2("/workdir/deseq2_output/samples_MasterTable.txt",sep = "\\t", row.names = 1)
    samples<-row.names(sampleTable)
    # dir
    dir <- "/workdir/kallisto_output"
    files<-file.path(dir, samples, "abundance.tsv")
    names(files)<-samples
    txi <- tximport(files, type = "kallisto", txOut = TRUE)
    txi <- summarizeToGene(txi, tx2gene = tx2gene)
    # model
    dds <- DESeqDataSetFromTximport(txi, sampleTable, ~ group )
    # if circRNA file is present, add circRNA counts to dds, else proceed
    # here, get count table, 
    # add circRNA counts
    circRNA_folder="${params.circRNA}"
    if(circRNA_folder != "None"){
      # gene count table
      gene_count = as.data.frame(counts(dds, normalized=FALSE))
      # circRNA table
      circRNA = read.delim('/workdir/${params.circRNA}/CircRNACount', check.names = FALSE, as.is = TRUE )
      coord = read.delim('/workdir/${params.circRNA}/CircCoordinates', check.names = FALSE, as.is = TRUE)
      
      # reformat circRNA header
      names(circRNA) <- gsub('.Chimeric.out.junction', '', names(circRNA))
      row.names(circRNA) <- paste0('circ_', coord[, 'Chr'], ':', coord[,'Start'], '-', coord[,'End'], '|', coord[,'Strand'], '|', coord[,'Gene'])
      
      circRNA = circRNA[,names(gene_count)]
      
      # add circRNAs to gene count table
      gene_count = rbind(gene_count, circRNA)
      # create DESeqDataSetFromMatrix
      dds <- DESeqDataSetFromMatrix(gene_count, sampleTable, ~ group)
    } 
    dds <- estimateSizeFactors(dds)
    res.counts <- counts(dds, normalized=TRUE)
    res_counts <- as.data.frame(res.counts)
    openxlsx::write.xlsx(res_counts, "/workdir/deseq2_output/all_res_counts.xlsx", row.names = TRUE, col.names = TRUE)
    result_tables <- list.files('/workdir/deseq2_output/', pattern = '.results.tsv')
    for(f in result_tables){
      tmp <- read.delim(paste0('/workdir/deseq2_output/', f))
      tmp <- tmp[,c('log2FoldChange', 'pvalue', 'padj')]
      names(tmp) <- paste(names(tmp), gsub('.results.tsv', '', gsub('Group_', '', f)), sep = '.')
      res_counts <- merge(res_counts, tmp, by.x = 'row.names', by.y = 'row.names')
      #print(nrow(res_counts))
      names(res_counts)[names(res_counts) == "Row.names"] <- "ensembl_gene_id"
      print(length(res_counts[,'ensembl_gene_id']))
      row.names(res_counts) <- res_counts[,'ensembl_gene_id']
      res_counts <- res_counts[, !duplicated(colnames(res_counts))]
    }
    # calculate fpkm values
    fpkm.deseq = as.data.frame(fpkm(dds))
    names(fpkm.deseq) = paste0('fpkm.', names(fpkm.deseq))
    res_counts = merge(res_counts, fpkm.deseq, by = 'row.names', all.x = TRUE)
    res_counts = res_counts[,-1]
    write.table(res_counts, "/workdir/deseq2_output/all_results_stats.tsv", sep = "\\t", quote = F, row.names = F)
    openxlsx::write.xlsx(res_counts, "/workdir/deseq2_output/all_results_stats.xlsx", row.names = FALSE, col.names = TRUE)

    sessionInfo()
  """
}

workflow {
  tx2gene()
  data = channel.fromPath( "${params.project_folder}/deseq2_output/*.input.tsv" )
  deseq2( data, tx2gene.out.collect() )
  mastertable( deseq2.out.collect() )
}