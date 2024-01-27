process pbptyper {
  tag           "$meta.id"
  label         "process_medium"
  stageInMode   "copy"
  publishDir    path: params.outdir, mode: 'copy'
  container     'staphb/pbptyper:1.0.4'
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
  path "pbptyper/*.tsv"                                  , emit: collect
  path "pbptyper/*"                                     , emit: all
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml"                          , emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p pbptyper logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : ${task.container}" >> $log_file
    pbptyper --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file
    
    pbptyper ${params.pbptyper_options} \
      --assembly ${contigs} \
      --prefix ${prefix} \
      --outdir pbptyper \
      | tee -a $log_file 
  """
}