process mykrobe {
  tag           "$meta.id"
  label         "process_medium"
  stageInMode   "copy"
  publishDir    path: params.outdir, mode: 'copy'
  container     'staphb/mykrobe:0.13.0'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '24h'

  when:
  flag =~ 'found'

  input:
  tuple val(meta), file(contigs), val(flag)

  output:
  path "mykrobe/*.csv"                                   , emit: collect
  path "mykrobe/*.json"
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml"                          , emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
      def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p tmp/ logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : ${task.container}" >> $log_file
    mykrobe --version >> $log_file
    mykrobe panels describe >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    mykrobe predict ${params.mykrobe_options} \
      --species tb \
      --prefix ${prefix} \
      --output mykrobe/${prefix} \
      --seq ${contigs} \
      --threads ${task.cpus} \
      --format json_and_csv \
      | tee -a $log_file
  """
}