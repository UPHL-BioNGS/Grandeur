process download_sra {
  tag           "${SRR}"
  label         "process_single"
  publishDir    params.outdir, mode: 'copy'
  container     'quay.io/biocontainers/sra-tools:2.11.0--pl5321ha49a11a_3'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '2h'
  
  input:
  val(SRR)

  output:
  tuple val(SRR), file("reads/${SRR}_{1,2}.fastq.gz")         , emit: fastq
  path "logs/${task.process}/${SRR}.${workflow.sessionId}.log", emit: log
  path  "versions.yml"                          , emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
      def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p reads logs/${task.process}
    log_file=logs/${task.process}/${SRR}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : ${task.container}" >> $log_file
    fasterq-dump --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    fasterq-dump ${params.fasterqdump_options} \
      ${SRR} \
      --split-files \
      --threads ${task.cpus} \
      --outdir reads | \
      tee -a $log_file

    gzip reads/${SRR}_1.fastq
    gzip reads/${SRR}_2.fastq
  """
}

