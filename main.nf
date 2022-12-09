#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process get_images {
  stageInMode 'symlink'
  stageOutMode 'move'

  script:
    """

    if [[ "${params.run_type}" == "r2d2" ]] || [[ "${params.run_type}" == "raven" ]] ; 

      then

        cd ${params.image_folder}

        if [[ ! -f rnaseq.python-3.8-1.sif ]] ;
          then
            singularity pull rnaseq.python-3.8-1.sif docker://index.docker.io/mpgagebioinformatics/rnaseq.python:3.8-1
        fi

        if [[ ! -f deseq2-1.38.0.sif ]] ;
          then
            singularity pull deseq2-1.38.0.sif docker://index.docker.io/mpgagebioinformatics/deseq2:1.38.0
        fi

        if [[ ! -f topgo-2.50.0.sif ]] ;
          then
            singularity pull topgo-2.50.0.sif docker://index.docker.io/mpgagebioinformatics/topgo:2.50.0
        fi

        if [[ ! -f cellplot-ea2dbc4.sif ]] ;
          then
            singularity pull cellplot-ea2dbc4.sif docker://index.docker.io/mpgagebioinformatics/cellplot:ea2dbc4
        fi

        if [[ ! -f rcistarget-1.17.0.sif ]] ;
          then
            singularity pull rcistarget-1.17.0.sif docker://index.docker.io/mpgagebioinformatics/rcistarget:1.17.0
        fi
    fi

    if [[ "${params.run_type}" == "local" ]] ; 

      then

        docker pull mpgagebioinformatics/rnaseq.python:3.8-1
        docker pull mpgagebioinformatics/deseq2:1.38.0
        docker pull mpgagebioinformatics/topgo:2.50.0
        docker pull mpgagebioinformatics/cellplot:ea2dbc4
        docker pull mpgagebioinformatics/rcistarget:1.17.0

    fi

    """

}

process tx2gene {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val gtf

  when:
    ( ! file("${params.project_folder}/deseq2_output/tx2gene.csv").exists() ) 
  
  script:
  """
    #!/usr/local/bin/python
    import AGEpy as age
    GTF=age.readGTF("/workdir/${gtf}")
    GTF["gene_id"]=age.retrieve_GTF_field(field="gene_id",gtf=GTF)
    GTF["transcript_id"]=age.retrieve_GTF_field(field="transcript_id",gtf=GTF)
    tx2gene=GTF[["transcript_id","gene_id"]].drop_duplicates().dropna()
    tx2gene.columns=["TXNAME","GENEID"]
    tx2gene[["TXNAME","GENEID"]].to_csv("/workdir/deseq2_output/tx2gene.csv", quoting=1, index=None)
  """
}

process annotations {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val gtf

  when:
    ( ! file("${params.project_folder}/deseq2_output/annotated/biotypes_go.txt").exists() ) 
  
  script:
  """
    #!/usr/local/bin/python
    import AGEpy as age
    import pandas as pd
    from biomart import BiomartServer
    attributes=["ensembl_gene_id","go_id","name_1006"]
    server = BiomartServer( "${params.biomart_host}" )
    organism=server.datasets["${params.biomart_dataset}"]
    response=organism.search({"attributes":attributes})
    response=response.content.decode().split("\\n")
    response=[s.split("\\t") for s in response ]
    bio_go=pd.DataFrame(response,columns=attributes)
    bio_go.to_csv("/workdir/deseq2_output/annotated/biotypes_go_raw.txt", index=None, sep="\t")
    bio_go.columns = ["ensembl_gene_id","GO_id","GO_term"]
    def CombineAnn(df):
        return pd.Series(dict(ensembl_gene_id = "; ".join([ str(s) for s in list(set(df["ensembl_gene_id"]))  if str(s) != "nan" ] ) ,\
                           GO_id = "; ".join([ str(s) for s in list(set(df["GO_id"])) if str(s) != "nan" ] ) ,\
                           GO_term = "; ".join([ str(s) for s in list(set(df["GO_term"])) if str(s) != "nan" ] ) ,\
                          ) )
    bio_go=bio_go.groupby(by="ensembl_gene_id", as_index=False).apply(CombineAnn)
    bio_go.reset_index(inplace=True, drop=True)
    
    GTF=age.readGTF("/workdir/${gtf}")
    GTF["gene_id"]=age.retrieve_GTF_field(field="gene_id",gtf=GTF)
    GTF["gene_biotype"]=age.retrieve_GTF_field(field="gene_biotype",gtf=GTF)
    GTF=GTF[["gene_id","gene_biotype"]].drop_duplicates()
    GTF.columns=["ensembl_gene_id","gene_biotype"]
    bio_go=pd.merge(GTF,bio_go,on=["ensembl_gene_id"],how="outer")
    bio_go.to_csv("/workdir/deseq2_output/annotated/biotypes_go.txt", sep= "\\t", index=None)
  """
}

process parse_submission {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val samplestable

  when:
    ( ! file("${params.project_folder}/deseq2_output/models.txt").exists() ) 
  
  script:
  """
    #!/usr/local/bin/python
    import pandas as pd
    from biomart import BiomartServer
    import itertools
    import AGEpy as age
    sfile="/workdir/${samplestable}"
    sdf=pd.read_excel(sfile)
    sam_df=sdf.copy()
    files_col=sam_df.columns.tolist()[0]
    cond_col=sam_df.columns.tolist()[1]
    sam_df.index=[s.split(".READ_1.fastq.gz")[0] for s in sam_df[files_col].tolist()]
    sam_df=sam_df[[cond_col]]
    sam_df.to_csv("/workdir/deseq2_output/samples_MasterTable.txt", sep="\t")
    fs=sdf.columns.tolist()[0]
    sdf[fs]=sdf[fs].apply(lambda x: x.split(".READ_1.fastq.gz")[0] )
    sdf.index=sdf[fs].tolist()
    sdf=sdf.drop([fs],axis=1)
    cols=sdf.columns.tolist()
    mods=[ [x,y] for x in cols for y in cols ]
    interactions=[]
    for c in mods:
        if c not in interactions:
            if [c[1], c[0]] not in interactions:
                if c[0] != c[1]:
                    interactions.append(c)
    single_models=cols
    textout=[]
    for m in single_models:
        variants=[ s for s in cols if s != m ] 
        tmp=sdf.copy()
        tmp["_group_"]=m
        for c in variants:
            tmp["_group_"]=tmp["_group_"]+"."+c+"_"+tmp[c]
        for g in list(set(tmp["_group_"].tolist())):
            outdf=tmp[tmp["_group_"]==g]
            model_data=list(set(outdf[m].tolist()))
            model_pairs=[ list(set([x,y])) for x in model_data for y in model_data ]
            model_pairs_=[ pair for pair in model_pairs if len(pair) > 1 ]
            model_pairs=[]
            for pair in model_pairs_:
                if pair not in model_pairs:
                    model_pairs.append(pair)
            for pair in model_pairs:
                outdf_=outdf[outdf[m].isin(pair)]
                outdf_=outdf_.drop(["_group_"],axis=1)
                coef=outdf_[m].tolist()
                ref=coef[0]
                target=[ t for t in coef if t != ref ][0]
                coef=m+"_"+str(target)+"_vs_"+str(ref)
                filename=coef+g.split(m)[-1]
                outdf_.to_csv("/workdir/deseq2_output/"+filename+".input.tsv", sep="\t")
                text=[ filename, m, g, str(ref), coef ]
                text="\\t".join(text)
                textout.append(text)
            
    with open("/workdir/deseq2_output/models.txt", "w") as mout:
        mout.write("\\n".join(textout) + "\\n")
  """
}

process tx2gene_proc {
  stageInMode 'symlink'
  stageOutMode 'move'

  output:
    val "/workdir/deseq2_output/deseq2.part1.Rdata", emit: rdata

  when:
    ( ! file("${params.project_folder}/deseq2_output/deseq2.part1.Rdata").exists() ) 

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
    ( ! file("${params.project_folder}/deseq2_output/${input_file}".replace(".input.tsv",".results.tsv")).exists() ) 

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


process annotator {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val gtf

  when:
    ( ! file("${params.project_folder}/deseq2_output/annotated/masterTable_annotated.xlsx").exists() ) 
  
  script:
  """
    #!/usr/local/bin/python
    import pandas as pd
    import os
    import AGEpy as age
    GTF=age.readGTF("/workdir/${gtf}")
    GTF["gene_id"]=age.retrieve_GTF_field(field="gene_id",gtf=GTF)
    GTF["gene_name"]=age.retrieve_GTF_field(field="gene_name",gtf=GTF)
    id_name=GTF[["gene_id","gene_name"]].drop_duplicates()
    id_name.reset_index(inplace=True, drop=True)
    id_name.columns=["ensembl_gene_id","gene_name"]
    #id_name=pd.read_table("/workdir/kallisto_index/cdna.norRNA.tsv")
    #id_name=id_name[["gene_id","gene_symbol","description"]]
    #id_name.columns=["ensembl_gene_id","gene_name","description"]
    #id_name=id_name.drop_duplicates()
    #id_name.reset_index(inplace=True,drop=True)
    if "${params.biomart_host}" != "None":
        bio_go=pd.read_csv("/workdir/deseq2_output/annotated/biotypes_go.txt", sep="\t")
    else:
        bio_go=pd.DataFrame(columns=["ensembl_gene_id"])
    deg_files=os.listdir("/workdir/deseq2_output/")
    deg_files=[ s for s in deg_files if "results.tsv" in s ]
    i=1
    s=[]
    dfs={}
    for f in deg_files:
        df=pd.read_table("/workdir/deseq2_output/"+f)
        df=pd.merge(id_name,df,left_on=["ensembl_gene_id"],right_index=True, how="right") # change to gene_id
        df=pd.merge(df,bio_go,on=["ensembl_gene_id"],how="left")
        df=df.sort_values(by=["padj"],ascending=True)
        df.to_csv("/workdir/deseq2_output/annotated/"+f, sep="\\t",index=None)
        df.to_excel("/workdir/deseq2_output/annotated/"+f.replace('.tsv', '.xlsx'), index=None)
        n=f.split(".results.tsv")[0]
        s.append([i,n])
        df=df[df["padj"]<0.05]
        df.reset_index(inplace=True, drop=True)
        dfs[i]=df
        i=i+1
    sdf=pd.DataFrame(s,columns=["sheet","comparison"])
    EXC=pd.ExcelWriter("/workdir/deseq2_output/annotated/significant.xlsx")
    sdf.to_excel(EXC,"summary",index=None)
    for k in list(dfs.keys()):
        dfs[k].to_excel(EXC, str(k),index=None)
    EXC.close()
    mt=pd.read_csv("/workdir/deseq2_output/all_results_stats.tsv", sep="\\t")
    mt_ann=pd.merge(id_name,mt,on=["ensembl_gene_id"], how="right")
    mt_ann=pd.merge(mt_ann,bio_go,on=["ensembl_gene_id"],how="left")
    mt_ann.to_csv("/workdir/deseq2_output/annotated/masterTable_annotated.tsv", sep="\\t",index=None)
    mt_ann.to_excel("/workdir/deseq2_output/annotated/masterTable_annotated.xlsx", index=None)
  """
}


process david_proc {
  stageInMode 'symlink'
  stageOutMode 'move'

  when:
    ( ! file("${params.project_folder}/deseq2_output/annotated/david.completed.txt").exists() )
  
  script:
  """
    #!/usr/local/bin/python
    import pandas as pd
    import AGEpy as age
    import os 
    import sys
    deseq2="/workdir/deseq2_output/annotated/"
    files=os.listdir(deseq2)
    files=[ s for s in files if ".results.tsv" in s ]
    for f in files:
        if os.path.isfile(deseq2+f.replace("results.tsv","DAVID.xlsx")):
            continue
        df=pd.read_csv(deseq2+f,sep="\t")
        df=df[df["padj"]<0.05]
        if df.shape[0] > 0:
            dics=df[["ensembl_gene_id","gene_name","log2FoldChange"]]
            dics["ensembl_gene_id"]=dics["ensembl_gene_id"].apply(lambda x: x.upper())
            dics.index=dics["ensembl_gene_id"].tolist()
            names_dic=dics[["gene_name"]].to_dict()["gene_name"]
            exp_dic=dics[["log2FoldChange"]].to_dict()["log2FoldChange"]
            
            genes=df["ensembl_gene_id"].tolist()
            
            DAVID=age.DAVIDenrich(database="${params.daviddatabase}",\
                        categories="GOTERM_BP_FAT,GOTERM_CC_FAT,GOTERM_MF_FAT,KEGG_PATHWAY,PFAM,PROSITE,GENETIC_ASSOCIATION_DB_DISEASE,OMIM_DISEASE",\
                      user="${params.DAVIDUSER}",\
                      ids=genes, verbose=True)
            
            if type(DAVID) == type(pd.DataFrame()):
                #for c in DAVID.columns.tolist():
                #    DAVID[c]=DAVID[c].apply(lambda x: x.decode())
                #print(DAVID.head(),DAVID["geneIds"].tolist(),names_dic )
                DAVID["genes name"]=DAVID["geneIds"].apply(lambda x: ", ".join([ str(names_dic[s.upper()]) for s in x.split(", ") ] ) )
                DAVID["log2fc"]=DAVID["geneIds"].apply(lambda x: ", ".join([ str(exp_dic[s.upper()]) for s in x.split(", ") ] ) )
                
                DAVID.to_csv(deseq2+f.replace("results.tsv","DAVID.tsv"), sep="\\t", index=None)
                EXC=pd.ExcelWriter(deseq2+f.replace("results.tsv","DAVID.xlsx"))
                df.to_excel(EXC,"genes",index=None)
                for cat in DAVID["categoryName"].tolist():
                    tmp=DAVID[DAVID["categoryName"]==cat]
                    tmp.to_excel(EXC,cat,index=None)
                EXC.close()
        else:
            print ("No significant differentially expressed genes")
    with open(deseq2+"david.completed.txt", "w") as fout:
      fout.write("completed")
  """
}


process topgo_proc {
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
    GOTERM_categories = list('BP' = 'GOTERM_BP',
                            'MF' = 'GOTERM_MF', 
                            'CC' = 'GOTERM_CC')
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
      
      names(D_out) <- c('categoryName', 'termName', 'listHits', 'percent', 'classicFisher', 'geneIds',
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

    if( "classicFisher" %in%  names(D) ){
      D\$ease <- as.numeric(as.character(D\$classicFisher))
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

process rcistarget_proc {
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


process qc_plots {
  stageInMode 'symlink'
  stageOutMode 'move'

  when:
    ( ! file("${params.project_folder}/qc_plots/pca_all_samples.pdf").exists() ) 
  
  script:
  """
  #!/bin/bash
  mkdir -p /workdir/qc_plots/
  QC_plots -o /workdir/qc_plots/ -de /workdir/deseq2_output/ -s /workdir/deseq2_output/samples_MasterTable.txt -t /workdir -sp ${params.spec}
  """
}


process get_ip {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    path cytoscape_host

  output:
    path cytoscape_host

  when:
    ( ! file("${params.project_folder}/deseq2_output/annotated/cytoscape.completed.txt").exists() )
  
  script:
  """
    #!/bin/bash
    echo "waiting for cytoscape to be available"
    while [[ ! -f /workdir/${cytoscape_host} ]] ; do 
      sleep 3\$((RANDOM % 9))
    done
    mv /workdir/${cytoscape_host} /workdir/${cytoscape_host}_inuse 
  """


}

process string {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    path cytoscape_host
  
  output:
    path cytoscape_host
  
  when:
      ( ! file("${params.project_folder}/deseq2_output/annotated/cytoscape.completed.txt").exists() ) 
  
  script:
  """
    #!/usr/local/bin/python
    import pandas as pd
    import numpy as np
    import AGEpy as age
    import sys
    import os
    from py2cytoscape import cyrest
    from py2cytoscape.cyrest.base import *
    import paramiko
    from time import sleep
    import matplotlib
    import matplotlib.pyplot as plt
    import tempfile
    ################# in values ################################
    with open("/workdir/${cytoscape_host}_inuse" , "r") as hostfile:
      host=hostfile.readlines()[0].split("\\n")[0]
    species="${params.species}"
    biomarthost="${params.biomart_host}"
    ###########################################################
    input_files=os.listdir("/workdir/deseq2_output/annotated/")
    input_files=[s for s in input_files if ".results.tsv" in s ]
    input_files=[ os.path.join("/workdir/deseq2_output/annotated/",s) for s in input_files if ".results.tsv" in s ]
    for fin in input_files:
        python_output="/".join(fin.split("/")[:-1])
        target=fin.replace("results.tsv","cytoscape")
        if os.path.isfile(target+".cys"):
           continue
        taxons={"caenorhabditis elegans":"6239","drosophila melanogaster":"7227",\
              "mus musculus":"10090","homo sapiens":"9606", "saccharomyces cerevisiae": "4932", "nothobranchius furzeri": "105023"}
        tags={"caenorhabditis elegans":"CEL","drosophila melanogaster":"DMEL",\
              "mus musculus":"MUS","homo sapiens":"HSA"}
        taxon_id=taxons[species]
        aging_genes = []
        ### ATTENTION ### if you are using yeast, you will need to uncomment the follwing lines 
        if species in tags.keys():
            organismtag=tags[species]
            
            if not os.path.isfile(python_output+"/homdf.txt"):
                print("Could not find ageing evidence table. Using biomart to create one.")
                sys.stdout.flush()
                homdf,HSA,MUS,CEL,DMEL=age.FilterGOstring(host=biomarthost)
                homdf.to_csv(python_output+"/homdf.txt", index=None,sep="\t")
            else:
                print("Found existing ageing evidence table.")
                sys.stdout.flush()
            homdf=pd.read_csv(python_output+"/homdf.txt", sep="\t")
            aging_genes=homdf[[organismtag+"_ensembl_gene_id","evidence"]].dropna()
            aging_genes=aging_genes[aging_genes[organismtag+"_ensembl_gene_id"]!="None"]
            aging_genes=aging_genes[organismtag+"_ensembl_gene_id"].tolist()
        ### till here
        dfin=pd.read_csv(fin, sep="\\t")
        cytoscape=cyrest.cyclient(host=host)
        cytoscape.version()
        cytoscape.session.new()
        # cytoscape.vizmap.apply(styles="default")
        # Annotate aging evindence
        def CheckEvidence(x,aging_genes=aging_genes):
            if x in aging_genes:
                res="aging_gene"
            else:
                res="no"
            return res
        ### also comment this line
        dfin["evidence"]=dfin["ensembl_gene_id"].apply(lambda x:CheckEvidence(x) )
        dfin["baseMean"]=dfin["baseMean"].apply(lambda x: np.log10(x))
        qdf=dfin[dfin["padj"]<0.05]
        if qdf.shape[0] == 0:
                sys.exit()
        qdf=qdf.sort_values(by=["padj"],ascending=True)
        query_genes=qdf["ensembl_gene_id"].tolist()[:1000]
        limit=int(len(query_genes)*.25)
        response=api("string", "protein query",\
                              {"query":",".join(query_genes),\
                              "cutoff":str(0.4),\
                              "species":species,\
                              "limit":str(limit),\
                              "taxonID":taxon_id},\
                              host=host, port="1234")
        cytoscape.layout.force_directed(defaultSpringCoefficient=".000004", defaultSpringLength="5")
        defaults_dic={"NODE_SHAPE":"ellipse",\
                      "NODE_SIZE":"60",\
                      "NODE_FILL_COLOR":"#AAAAAA",\
                      "EDGE_TRANSPARENCY":"120"}
        defaults_list=cytoscape.vizmap.simple_defaults(defaults_dic)
        NODE_LABEL=cytoscape.vizmap.mapVisualProperty(visualProperty="NODE_LABEL",\
                                                      mappingType="passthrough",\
                                                      mappingColumn="display name")
        cytoscape.vizmap.create_style(title="dataStyle",\
                                      defaults=defaults_list,\
                                      mappings=[NODE_LABEL])
        sleep(4)
        cytoscape.vizmap.apply(styles="dataStyle")
        uploadtable=dfin[dfin["padj"]<0.05][["ensembl_gene_id","baseMean","log2FoldChange","evidence"]].dropna()
        # uploadtable=dfin[dfin["padj"]<0.05][["ensembl_gene_id","baseMean","log2FoldChange"]].dropna() ### use this line if you are using yeast
        cytoscape.table.loadTableData(uploadtable,df_key="ensembl_gene_id",table_key_column="query term")
        sleep(10)
        cmap = matplotlib.cm.get_cmap("bwr")
        norm = matplotlib.colors.Normalize(vmin=-4, vmax=4)
        min_color=matplotlib.colors.rgb2hex(cmap(norm(-4)))
        center_color=matplotlib.colors.rgb2hex(cmap(norm(0)))
        max_color=matplotlib.colors.rgb2hex(cmap(norm(4)))  
        NODE_FILL_COLOR=cytoscape.vizmap.mapVisualProperty(visualProperty="NODE_FILL_COLOR",mappingType="continuous",\
                                                          mappingColumn="log2FoldChange",\
                                                        lower=[-4,min_color],\
                                                          center=[0.0,center_color],\
                                                          upper=[4,max_color])
        ### do not do this if you are using yeast ...
        # apply diamond shape and increase node size to nodes with aging evidence
        NODE_SHAPE=cytoscape.vizmap.mapVisualProperty(visualProperty="NODE_SHAPE",mappingType="discrete",mappingColumn="evidence",\
                                                      discrete=[ ["aging_gene","no"], ["DIAMOND", "ellipse"] ])
        NODE_SIZE=cytoscape.vizmap.mapVisualProperty(visualProperty="NODE_SIZE",mappingType="discrete",mappingColumn="evidence",\
                                                    discrete=[ ["aging_gene","no"], ["100.0","60.0"] ])
        ###
        cytoscape.vizmap.update_style("dataStyle",mappings=[NODE_SIZE,NODE_SHAPE,NODE_FILL_COLOR])
        # cytoscape.vizmap.update_style("dataStyle",mappings=[NODE_FILL_COLOR]) # if using yeast
        cytoscape.vizmap.apply(styles="dataStyle")
        network = "current"
        namespace='default'
        PARAMS=set_param(["columnList","namespace","network"],["SUID",namespace,network])
        network=api(namespace="network", command="get attribute",PARAMS=PARAMS, host=host,port='1234',version='v1')
        network=int(network[0]["SUID"])
        basemean = cytoscape.table.getTable(table="node",columns=["baseMean"], network = network)
        min_NormInt = min(basemean.dropna()["baseMean"].tolist())
        max_NormInt = max(basemean.dropna()["baseMean"].tolist())
        cent_NormInt = np.mean([min_NormInt,max_NormInt])
        cmap = matplotlib.cm.get_cmap("Reds")
        norm = matplotlib.colors.Normalize(vmin=min_NormInt, vmax=max_NormInt)
        min_color=matplotlib.colors.rgb2hex(cmap(norm(np.mean([min_NormInt,max_NormInt]))))
        center_color=matplotlib.colors.rgb2hex(cmap(norm(cent_NormInt)))
        max_color=matplotlib.colors.rgb2hex(cmap(norm(max_NormInt)))  
        NODE_BORDER_PAINT=cytoscape.vizmap.mapVisualProperty(visualProperty="NODE_BORDER_PAINT",\
                                                            mappingType="continuous",\
                                                            mappingColumn="baseMean",\
                                                            lower=[min_NormInt,min_color],\
                                                            center=[np.mean([min_NormInt,max_NormInt]),center_color],\
                                                            upper=[max_NormInt,max_color])
        cytoscape.vizmap.update_style("dataStyle",mappings=[NODE_BORDER_PAINT])
        NODE_BORDER_WIDTH=cytoscape.vizmap.mapVisualProperty(visualProperty="NODE_BORDER_WIDTH",\
                                                            mappingType="continuous",\
                                                            mappingColumn="baseMean",\
                                                            lower=[min_NormInt,2],\
                                                            center=[np.mean([min_NormInt,max_NormInt]),4],\
                                                            upper=[max_NormInt,8])
        cytoscape.vizmap.update_style("dataStyle",mappings=[NODE_BORDER_WIDTH])
        cytoscape.vizmap.apply(styles="dataStyle")
        cytoscape.network.rename(name="main String network")
        cytoscape.network.select(edgeList="all", extendEdges="true")
        cytoscape.network.create(source="current",nodeList="selected")
        cytoscape.network.rename(name="main String network (edges only)")
        cytoscape.network.set_current(network="main String network (edges only)")
        log2FoldChange = cytoscape.table.getTable(table="node",columns=["log2FoldChange"])
        if int(len(log2FoldChange)*.10) > 0:
            log2FoldChange["log2FoldChange"]=log2FoldChange["log2FoldChange"].apply(lambda x: abs(x))
            log2FoldChange=log2FoldChange.sort_values(by=["log2FoldChange"],ascending=False)
            top_nodes=log2FoldChange.index.tolist()[:int(len(log2FoldChange)*.10)]
            cytoscape.network.set_current(network="main String network (edges only)")
            cytoscape.network.select(nodeList="name:"+",".join(top_nodes))
            cytoscape.network.select(firstNeighbors="any",network="current")
            sleep(5)
            cytoscape.network.create(source="current",nodeList="selected")
            cytoscape.network.rename(name="top "+str(int(len(log2FoldChange)*.10))+" changed firstNeighbors")
        def MAKETMP():
            (fd, f) = tempfile.mkstemp()
            f="/tmp/"+f.split("/")[-1]
            return f
        cys=MAKETMP()
        cyjs=MAKETMP()
        main_png=MAKETMP()
        main_pdf=MAKETMP()
        edg_png=MAKETMP()
        edg_pdf=MAKETMP()
        neig_png=MAKETMP()
        neig_pdf=MAKETMP()
        cytoscape.session.save_as(session_file=cys)
        cytoscape.network.export(options="CYJS",OutputFile=cyjs)
        cytoscape.network.set_current(network="main String network")
        cytoscape.network.deselect(edgeList="all",nodeList="all")
        cytoscape.view.export(options="PNG",outputFile=main_png)
        cytoscape.view.export(options="PDF",outputFile=main_pdf)
        cytoscape.network.set_current(network="main String network (edges only)")
        cytoscape.network.deselect(edgeList="all",nodeList="all")
        cytoscape.view.export(options="PNG",outputFile=edg_png)
        cytoscape.view.export(options="PDF",outputFile=edg_pdf)
        if int(len(log2FoldChange)*.10) > 0:
            cytoscape.network.set_current(network="top "+str(int(len(log2FoldChange)*.10))+" changed firstNeighbors")
            cytoscape.network.deselect(edgeList="all",nodeList="all")
            sleep(5)
            cytoscape.view.export(options="PNG",outputFile=neig_png)
            cytoscape.view.export(options="PDF",outputFile=neig_pdf)
        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh.connect(host, username="bioinf")
        ftp_client=ssh.open_sftp()
        for f, extension, local in zip([cys,cyjs,main_png,main_pdf,edg_png,edg_pdf,neig_png,neig_pdf],\
                                        [".cys",".cyjs",".png",".pdf",".png",".pdf",".png",".pdf" ],\
                                        [target+".cys",target+".cyjs",target+".main.png",target+".main.pdf",\
                                        target+".main.edges.png",target+".main.edges.pdf",\
                                        target+".topFirstNeighbors.png",target+".topFirstNeighbors.pdf"]):
            try:
                ftp_client.get(f+extension,local)
                ssh_stdin, ssh_stdout, ssh_stderr = ssh.exec_command("rm "+f+extension )
            except:
                print("No "+local)
                sys.stdout.flush()
        print(f"Done with cytoscape for {fin}.")
    with open("/workdir/deseq2_output/annotated/cytoscape.completed.txt", "w") as f:
      f.write("cytoscape completed")

    os.rename("/workdir/${cytoscape_host}_inuse","/workdir/${cytoscape_host}")
  """
}

process release_ip {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    path cytoscape_host

  when:
    file("${params.project_folder}/deseq2_output/annotated/cytoscape.completed.txt").exists()
  
  script:
  """
    #!/bin/bash
    mv /workdir/${cytoscape_host}_inuse /workdir/${cytoscape_host} 
  """

}

process upload_paths {
  stageInMode 'symlink'
  stageOutMode 'move'

  script:
  """
    rm -rf upload.txt

    cd ${params.project_folder}/deseq2_output/annotated

    for f in \$(ls *.results.xlsx) ; do echo "deseq2 \$(readlink -f \${f})" >>  upload.txt_ ; done
    echo "deseq2 \$(readlink -f significant.xlsx)" >>  upload.txt_
    echo "deseq2 \$(readlink -f masterTable_annotated.xlsx)" >>  upload.txt_

    if [[ \$(ls  | grep cytoscape) ]] ; then
      for f in \$(ls *.cytoscape.* ) ; do echo "cytoscape \$(readlink -f \${f})" >>  upload.txt_ ; done
    fi

    if [[ \$(ls  | grep DAVID) ]] ; then
      for f in \$(ls *.DAVID.* ) ; do echo "david \$(readlink -f \${f})" >>  upload.txt_ ; done
    fi

    if [[ \$(ls  | grep RcisTarget) ]] ; then
      for f in \$(ls *.RcisTarget.* ) ; do echo "rcistarget \$(readlink -f \${f})" >>  upload.txt_ ; done
    fi

    if [[ \$(ls  | grep topGO) ]] ; then
      for f in \$(ls *.topGO.* ) ; do echo "togo \$(readlink -f \${f})" >>  upload.txt_ ; done
    fi

    uniq upload.txt_ upload.txt 
    rm upload.txt_

    cd ${params.project_folder}/qc_plots
    rm -rf upload.txt 
    for f in \$(ls *.* | grep -v upload.txt) ; do echo "qc_plots \$(readlink -f \${f})" >>  upload.txt_ ; done

    uniq upload.txt_ upload.txt 
    rm upload.txt_

  """
}

workflow images {
  main:
    get_images()
}

workflow upload {
  main:
    upload_paths()
}

workflow preprocess {

  folder=file("${params.project_folder}/deseq2_output/annotated")
  if( ! folder.isDirectory() ) {
    folder.mkdirs()
  }

  if ( 'gtf' in params.keySet() ) {
    gtf=${params.gtf}
  } else {
    gtf="kallisto_index/${params.organism}.${params.release}.no.rRNA.gtf"
  }

  tx2gene(gtf)

  if ( 'biomart_host' in params.keySet() ) { 
    annotations(gtf)
  }

  if ( 'samplestable' in params.keySet() ) {
    samplestable=${params.samplestable}
  } else {
    samplestable="sample_sheet.xlsx"
  }
  parse_submission(samplestable) 

}


workflow pairwise {
  tx2gene_proc()
  data = channel.fromPath( "${params.project_folder}/deseq2_output/*.input.tsv" )
  deseq2( data, tx2gene_proc.out.collect() )
  mastertable( deseq2.out.collect() )
}


workflow annotate {
  if ( 'gtf' in params.keySet() ) {
    gtf=${params.gtf}
  } else {
    gtf="kallisto_index/${params.organism}.${params.release}.no.rRNA.gtf"
  }

  annotator(gtf) 

} 

workflow david {
  // data = channel.fromPath( "${params.project_folder}/deseq2_output/annotated/*.results.tsv" )
  david_proc() 
}

workflow topgo {
  data = channel.fromPath( "${params.project_folder}/deseq2_output/annotated/*.results.tsv" )
  topgo_proc(data) 
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
  cellplot(topGO,"tsv", "GOTERM_BP", "15") 
}

workflow cellplots {
  david_bp()
  david_kegg()
  topgo_bp()
}

workflow rcistarget {
  if ( 'rcis_db' in params.keySet() ) {
    data = channel.fromPath( "${params.project_folder}/deseq2_output/annotated/*.results.tsv" )
    rcistarget_proc(data) 
  }
}

workflow qc {
  qc_plots()
}

workflow string_cytoscape {
  if ( 'cytoscape_host' in params.keySet() ) {
    if ( "${params.cytoscape_host}" != "None" ) {
        printf "${params.cytoscape_host} will be renamed to ${params.cytoscape_host}_inuse"
        get_ip("${params.cytoscape_host}")
        string(get_ip.out.collect())
        release_ip(string.out.collect()) 
    }
  }

} 
