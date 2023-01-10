process download_sra {
  tag           "${SRR}"
  //errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/sratoolkit:2.9.2'
  maxForks      10
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-large'
  
  input:
  val(SRR)

  output:
  tuple val(sample), file("reads/${SRR}_*.fastq.gz")          , emit: fastq
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

