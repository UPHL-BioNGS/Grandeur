process summary {
  tag           "Creating summary files"
  publishDir    params.outdir, mode: 'copy'
  container     'quay.io/uphl/seaborn:0.12.2-2'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '10m'

  input:
  file(input)

  output:
  path "grandeur_summary.tsv"                                  , emit: summary_tsv
  path "grandeur_summary.txt"                                  , emit: summary_txt
  path "summary/grandeur_extended_summary.tsv"                 , emit: extended_tsv
  path "summary/grandeur_extended_summary.txt"                 , emit: extended_txt
  path "logs/${task.process}/summary.${workflow.sessionId}.log", emit: log

  shell:
  '''
    mkdir -p summary logs/!{task.process}
    log_file=logs/!{task.process}/summary.!{workflow.sessionId}.log

    date > $log_file
    echo "container : !{task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    python summary.py | tee -a $log_file
  '''
}

process names {
  tag           "${sample}"
  container     'quay.io/uphl/seaborn:0.12.2-2'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '10m'
  
  input:
  tuple val(sample), file(input)

  output:
  path "summary/${sample}_names.csv"                              , emit: collect
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log" , emit: log

  shell:
  '''
    mkdir -p summary logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    date > $log_file
    echo "container : !{task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    echo "sample,file,version" > summary/!{sample}_names.csv
    echo "!{sample},!{input},!{workflow.manifest.version}" >> summary/!{sample}_names.csv
  '''
}