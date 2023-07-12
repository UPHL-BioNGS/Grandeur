process datasets_summary {
  tag           "${taxon}"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/ncbi-datasets:15.2.0'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '10m'
  
  input:
  val(taxon)

  output:
  path "datasets/*_genomes.csv"                                 , emit: genomes
  path "logs/${task.process}/*.${workflow.sessionId}.log", emit: log

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

    taxon="$(echo !{taxon} | tr '_' ' ' | sed 's/[//g' | sed 's/]//g' )"
    echo "the taxon is now $taxon"

    datasets summary genome taxon "$taxon" --reference --limit !{params.datasets_max_genomes} --as-json-lines | \
      dataformat tsv genome --fields accession,assminfo-refseq-category,assminfo-level,organism-name,assmstats-total-ungapped-len | \
      grep -v Homo | \
      tr '\\t' ',' \
      > datasets/!{taxon}_genomes.csv

    datasets summary genome taxon "$taxon" --limit !{params.datasets_max_genomes} --as-json-lines | \
      dataformat tsv genome --fields accession,assminfo-refseq-category,assminfo-level,organism-name,assmstats-total-ungapped-len | \
      grep -v Homo | \
      grep -v "Assembly Accession" | \
      tr '\\t' ',' \
      >> datasets/!{taxon}_genomes.csv
  '''
}

// It is faster if datasets can download the entire list at a time, but there is a timeout for downloading that is about 20 minutes.
// The '||' is to allow each genome to be downloaded on its own, which is longer overall but each genome should be less than 20 minutes.
process datasets_download {
  tag           "Downloading Genomes"
  // because there's no way to specify threads
  label         "medcpus"
  publishDir = [ path: "${params.outdir}", mode: 'copy', pattern: "{logs/*/*log,datasets/fastani_refs.tar.gz}" ]
  container     'staphb/ncbi-datasets:15.2.0'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '1h'
  
  input:
  file(ids)
  file(genomes)

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

    cut -f 1 !{genomes} > all_runs.txt
    grep -h -v Accession !{ids} | cut -f 1 -d , | sort | uniq > this_run.txt

    cat all_runs.txt this_run.txt | sort | uniq > id_list.txt

    ( datasets download genome accession --inputfile id_list.txt --filename ncbi_dataset.zip unzip -o ncbi_dataset.zip ) || \
    ( while read line ; do echo "Downloading $line" ; datasets download genome accession $line --filename dataset.zip ; unzip -o dataset.zip ; done < id_list.txt )

    fastas=$(ls ncbi_dataset/data/*/*.fna )

    for fasta in ${fastas[@]}
    do
      echo "Copying $fasta to genomes"
      accession=$(echo $fasta | cut -f 4 -d / | cut -f 1,2 -d _ )
      organism=$(head -n 1 $fasta | awk '{print $2 "_" $3 }' | sed 's/,//g' )
      cat $fasta | sed 's/ /_/g' | sed 's/,//g' > genomes/${organism}_${accession}.fna
    done  

    rm -rf genomes/*:_*

    tar -czvf datasets/fastani_refs.tar.gz genomes/
  '''
}
