process prokka {
  tag           "${sample}"
  label         "maxcpus"
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/prokka:1.14.5'
  maxForks      10
  
  input:
  tuple val(sample), file(contigs), val(genus), val(species)

  output:
  path "prokka/${sample}/*"                                      , emit: prokka_files
  path "prokka/${sample}/${sample}.txt"                          , emit: for_multiqc
  path "gff/${sample}.gff"                                       , emit: gffs
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log", emit: log

  shell:
  '''
    mkdir -p prokka gff logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    prokka -v >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    prokka !{params.prokka_options} \
      --cpu !{task.cpus} \
      --outdir prokka/!{sample} \
      --prefix !{sample} \
      --genus !{genus} \
      --species !{species} \
      --force !{contigs} \
      | tee -a $log_file

    cp prokka/!{sample}/!{sample}.gff gff/!{sample}.gff
  '''
}
