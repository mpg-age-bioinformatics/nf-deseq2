#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process tx2gene {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    path gtf

  when:
    ( ! file("/workdir/deseq2_output/tx2gene.csv").exists() ) 
  
  script:
    """
    #!/usr/local/bin/python
    GTF=age.readGTF("${gtf}")
    GTF["gene_id"]=age.retrieve_GTF_field(field="gene_id",gtf=GTF)
    GTF["transcript_id"]=age.retrieve_GTF_field(field="transcript_id",gtf=GTF)
    tx2gene=GTF[["transcript_id","gene_id"]].drop_duplicates().dropna()
    tx2gene.columns=["TXNAME","GENEID"]
    tx2gene[["TXNAME","GENEID"]].to_csv("/workdir/deseq2_output/tx2gene.csv", quoting=1, index=None)
  """
}

workflow {
  folder=file("${params.project_folder}/deseq2_output/")
  if( ! folder.isDirectory() ) {
    folder.mkdirs()
  }

  if ( 'gtf' in params.keySet() ) {
    gtf=${params.gtf}
  } else {
    gtf="/workdir/kallisto_index/${params.organism}.${params.release}.no.rRNA.gtf"
  }
  tx2gene(gtf)
  

} 