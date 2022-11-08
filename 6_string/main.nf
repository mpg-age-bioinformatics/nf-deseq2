#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process get_ip {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    path cytoscape_host

  output:
    path cytoscape_host
  
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
  """
}

process release_ip {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    path cytoscape_host
  
  script:
  """
    #!/bin/bash
    mv /workdir/${cytoscape_host}_inuse /workdir/${cytoscape_host} 
  """

}

workflow {
  if ( 'cytoscape_host' in params.keySet() ) {
    if ( "${params.cytoscape_host}" != "None" ) {
        printf "${params.cytoscape_host} will be renamed to ${params.cytoscape_host}_inuse"
        get_ip("${params.cytoscape_host}")
        string(get_ip.out.collect())
        release_ip(string.out.collect()) 
    }
  }

} 