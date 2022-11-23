process fasterqdump {
  tag "${SRR}"
  label "maxcpus"
  errorStrategy { task.exitStatus == 21 ? 'ignore' : 'terminate' }

  input:
  val(SRR)

  output:
  path "reads/${sample}/*"                                   , emit: reads
  path "logs/${task.process}/${SRR}.${workflow.sessionId}.log", emit: log

  shell:
  '''
    mkdir -p reads logs/!{task.process}
    log_file=logs/!{task.process}/!{SRR}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    fasterq-dump --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    fasterq-dump \
        -A !{SRR} \
        --split-files \
        --threads !{task.cpus} \
        --outdir reads \
        | tee -a $log_file
  '''
}
