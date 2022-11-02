process mlst {
  tag "${sample}"

  when:
  params.contig_processes =~ /mlst/

  input:
  tuple val(sample), file(contig)

  output:
  path "mlst/${sample}_mlst.txt"                                       , emit: collect
  tuple val(sample), env(mlst)                                         , emit: mlst
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log", emit: log

  shell:
  '''
    mkdir -p mlst logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    mlst --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    mlst !{params.mlst_options} !{contig} > mlst/!{sample}_mlst.txt

    mlst=$(awk '{ print $2 ":" $3 }' mlst/!{sample}_mlst.txt)
  '''
}
