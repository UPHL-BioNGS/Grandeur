process datasets_summary {
  tag           "${taxon}"
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/ncbi-datasets:14.3.0'
  maxForks      10
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-large'
  
  input:
  val(taxon)

  output:
  path "datasets/${taxon}_genomes.csv"                          , emit: genomes
  path "logs/${task.process}/${taxon}.${workflow.sessionId}.log", emit: log

  shell:
  '''
    mkdir -p datasets logs/!{task.process}
    log_file=logs/!{task.process}/!{taxon}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    datasets --version 2>> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    taxon=$(echo !{taxon} | tr "_" " ")

    datasets summary genome taxon "$taxon" --reference  --limit !{params.datasets_max_genomes} --as-json-lines | \
      dataformat tsv genome --fields accession,assminfo-refseq-category,assminfo-level,organism-name,assmstats-total-ungapped-len | \
      tr '\\t' ',' \
      > datasets/!{taxon}_genomes.csv
  '''
}

process datasets_download {
  tag           "Downloading Genomes"
  // because there's no way to specify threads
  label         "medcpus"
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/ncbi-datasets:14.3.0'
  maxForks      10
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-large'
  
  input:
  file(ids)

  output:
  path "datasets/fastani_refs.tar.gz"                                    , emit: tar
  path "genomes"                                                         , emit: genomes
  path "logs/${task.process}/datasets_download.${workflow.sessionId}.log", emit: log

  shell:
  '''
    mkdir -p datasets genomes logs/!{task.process}
    log_file=logs/!{task.process}/datasets_download.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    datasets --version 2>> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    grep -h -v Accession !{ids} | cut -f 1 -d , | sort | uniq > id_list.txt

    datasets download genome accession --inputfile id_list.txt --filename ncbi_dataset.zip

    unzip -o ncbi_dataset.zip
    
    fastas=$(ls ncbi_dataset/data/*/*.fna)
    
    for fasta in ${fastas[@]}
    do
      accession=$(echo $fasta | cut -f 3 -d / )
      organism=$(grep -h $accession !{ids} | cut -f 4 -d , | sed 's/ /_/g' )
      cat $fasta | sed 's/ /_/g' | sed 's/,//g' > genomes/${organism}_${accession}.fna
    done

    tar -czvf datasets/fastani_refs.tar.gz genomes/
  '''
}