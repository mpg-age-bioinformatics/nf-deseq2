#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process david {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val f

  when:
    ( ! file("${params.project_folder}/deseq2_output/annotated/${f}".replace(".results.tsv",".DAVID.tsv")).exists() ) 
  
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
  """
}

workflow {
  david() 
} 