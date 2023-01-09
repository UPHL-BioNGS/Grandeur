process fasterqdump {
  tag           "${SRR}"
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  publishDir    params.outdir, mode: 'copy'
  container     'ncbi/sra-tools:3.0.1'
  maxForks      10
  
  input:
  val(SRR)

  output:
  tuple val(sample), file("reads/${SRR}_*.fastq.gz")            , emit: reads
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
      tee -a $log_file
  '''
}
