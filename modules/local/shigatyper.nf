process shigatyper {
  tag           "${meta.id}"
  label         "process_medium"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/shigatyper:2.0.5'
  stageInMode   'copy'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-large'
  //#UPHLICA memory 26.GB
  //#UPHLICA cpus 7
  //#UPHLICA time '10m'
  
  input:
  tuple val(meta), file(input), val(flag), file(script)

  output:
  path "shigatyper/*_shigatyper.tsv",      optional: true, emit: files
  path "shigatyper/*_shigatyper-hits.tsv", optional: true, emit: collect
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  (task.ext.when == null || task.ext.when) && flag =~ 'found'

  shell:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p shigatyper logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : ${task.container}" >> $log_file
    shigatyper --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    shigatyper ${params.shigatyper_options} \
      --SE ${input} \
      --name ${prefix} \
      | tee -a $log_file

    python3 ${script} ${prefix}-hits.tsv shigatyper/${prefix}_shigatyper-hits.tsv shigatyper ${prefix}

    if [ -f "${prefix}.tsv" ] ; then cp ${prefix}.tsv shigatyper/${prefix}_shigatyper.tsv ; fi
  """
}
