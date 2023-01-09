process summary {
  tag           "Creating summary files"
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  publishDir    params.outdir, mode: 'copy'
  container     'quay.io/biocontainers/pandas:1.1.5'
  maxForks      10
  
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

    python summary.py

    mv *extended* summary/.
  '''
}

process names {
  tag           "${sample}"
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  publishDir    params.outdir, mode: 'copy'
  container     'quay.io/biocontainers/pandas:1.1.5'
  maxForks      10
  
  input:
  tuple val(sample), file(input), val(reads), val(phix)

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

    echo "sample,file,version,reads,phix_reads" > summary/!{sample}_names.csv
    echo "!{sample},!{input},!{workflow.manifest.version},!{reads},!{phix}" >> summary/!{sample}_names.csv
  '''
}