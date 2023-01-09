process multiqc {
  tag           "multiqc"
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  publishDir    params.outdir, mode: 'copy'
  container     'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0'
  maxForks      10
  
  input:
  file(input)

  output:
  path "multiqc/multiqc_report.html", optional: true                   , emit: report
  path "multiqc/multiqc_data/*"     , optional: true                   , emit: data_folder
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log", emit: log_files

  shell:
  '''
    mkdir -p multiqc quast logs/!{task.process}
    log_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    multiqc --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    for quast_file in $(ls *_quast_report.tsv)
    do
      sample=$(echo $quast_file | sed 's/_quast_report.tsv//g' | head -n 1 )
      mkdir -p quast/$sample
      mv $quast_file quast/$sample/report.tsv
    done

    multiqc !{params.multiqc_options} \
      --outdir multiqc \
      --cl_config "prokka_fn_snames: True"  \
      . \
      | tee -a $log_file
  '''
}
