#!/usr/bin/env nextflow
nextflow.enable.dsl=2

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

workflow {

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