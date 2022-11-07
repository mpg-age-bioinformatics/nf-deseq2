#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process annotate {
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

workflow {
  if ( 'gtf' in params.keySet() ) {
    gtf=${params.gtf}
  } else {
    gtf="kallisto_index/${params.organism}.${params.release}.no.rRNA.gtf"
  }

  annotate(gtf) 

} 