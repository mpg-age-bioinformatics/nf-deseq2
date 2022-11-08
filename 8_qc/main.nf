#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process qc {
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

workflow {
  qc()
} 