process elgato {
  tag           "$meta.id"
  label         "process_medium"
  publishDir    path: params.outdir, mode: 'copy'
  container     'staphb/elgato:1.15.2'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '24h'

  when:
  flag =~ 'found'

  input:
  tuple val(sample), file(contigs), val(flag)

  output:
  path "legsta/${sample}_legsta.csv"                             , emit: collect
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log", emit: log
  path  "versions.yml"                          , emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
      def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
  '''
    mkdir -p elgato logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    elgato --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file
    
    el_gato.py !{params.elgato_options} \
        --assembly !{contigs} \
        --out elgato/!{sample} \
        --threads !{task.cpus} \
        | tee -a $log_file

    exit 1
  '''
}