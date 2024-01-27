process kaptive {
  tag           "$meta.id"
  label         "process_medium"
  stageInMode   "copy"
  publishDir    path: params.outdir, mode: 'copy'
  container     'staphb/kaptive:2.0.5'
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
  path "kaptive/${prefix}_kaptive.csv"                           , emit: collect
  path "logs/${task.process}/${prefix}.${workflow.sessionId}.log", emit: log
  path  "versions.yml"                          , emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
      def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
  """
    echo "1"
    mkdir -p kaptive logs/${task.process}
    echo "2"
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    echo "3"
    # time stamp + capturing tool versions
    echo "4"
    date > $log_file
    echo "5"
    echo "container : ${task.container}" >> $log_file
    echo "6"
    kaptive --version >> $log_file
    echo "7"
    echo "Nextflow command : " >> $log_file
    echo "8"
    cat .command.sh >> $log_file
    
    echo "9"
    for reference in /kaptive/reference_database/*
    do
      echo $reference    
        kaptive ${params.kaptive_options} \
            --k_refs $reference \
            --assembly ${contigs} \
            | tee -a $log_file 

    done
    exit 1
  """
}