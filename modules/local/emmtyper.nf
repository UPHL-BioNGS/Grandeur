process emmtyper {
  tag           "$meta.id"
  label         "process_medium"
  stageInMode   "copy"
  publishDir    path: params.outdir, mode: 'copy'
  container     'staphb/emmtyper:0.2.0'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '24h'
  
  when:
  flag =~ 'found'

    when:
  task.ext.when == null || task.ext.when

  input:
  tuple val(meta), file(contigs), val(flag), file(script)

  output:
  path "emmtyper/${prefix}_emmtyper.txt"                         , emit: collect
  path "emmtyper/*"                                              , emit: everything
  path "logs/${task.process}/${prefix}.${workflow.sessionId}.log", emit: log
  path  "versions.yml"                          , emit: versions

  shell:
      def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p emmtyper logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    PATH=/opt/conda/envs/emmtyper/bin:$PATH
    date > $log_file
    echo "container : ${task.container}" >> $log_file
    emmtyper --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file
    
    emmtyper ${params.emmtyper_options} \
      --output-format 'verbose' \
      ${contigs} \
      | tee -a $log_file \
      > ${prefix}_emmtyper.txt

    python3 ${script} ${prefix}_emmtyper.txt emmtyper/${prefix}_emmtyper.txt emmtyper ${prefix}
  """
}
