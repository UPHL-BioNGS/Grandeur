process phytreeviz {
  tag           "${analysis}"
  label         "process_medium"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/phytreeviz:0.1.0'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA cpus 14
  //#UPHLICA memory 60.GB
  //#UPHLICA time '24h'
  
  input:
  tuple val(analysis), file(newick)

  output:
  path "phytreeviz/${analysis}_tree.png", emit: for_multiqc
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p phytreeviz logs/${task.process}
    log_file=logs/${task.process}/${analysis}.${task.process}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    phytreeviz -v >> $log_file
    echo "container : ${task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    phytreeviz ${params.phytreeviz_options} \
        -i ${newick} \
        -o phytreeviz/${analysis}_tree.png
  """
}
