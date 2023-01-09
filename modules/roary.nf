process roary {
  tag           "Core Genome Alignment"
  label         'maxcpus'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/roary:3.13.0'
  maxForks      10
  
  input:
  file(contigs)

  output:
  path "roary/*"                                                       , emit: roary_files
  path "roary/fixed_input_files/*"                                     , emit: roary_input_files
  path "roary/core_gene_alignment.aln"                                 , emit: core_gene_alignment
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log", emit: log_files

  shell:
  '''
    mkdir -p logs/!{task.process}
    log_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    roary -a >> $log_file
    echo "container : !{task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    echo "There are $(ls *gff | wc -l) files for alignment" >> $log_file

    roary !{params.roary_options} \
      -p !{task.cpus} \
      -f roary \
      -e -n \
      *.gff \
      | tee -a $log_file

    # do something with summary
    # check number of samples
    # check number of core genes
    exit 1
  '''
}
