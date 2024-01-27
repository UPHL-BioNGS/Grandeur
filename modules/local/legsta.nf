process legsta {
  tag           "$meta.id"
  label         "process_medium"
  stageInMode   "copy"
  publishDir    path: params.outdir, mode: 'copy'
  container     'staphb/legsta:0.5.1'
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
  tuple val(meta), file(contigs), val(flag)

  output:
  path "legsta/${prefix}_legsta.csv"                             , emit: collect
  path "logs/${task.process}/${prefix}.${workflow.sessionId}.log", emit: log
  path  "versions.yml"                          , emit: versions

  shell:
      def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p legsta logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : ${task.container}" >> $log_file
    legsta --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file
    
    cat ${contigs} | awk '{($1=$1); print $0}' > ${prefix}.fasta.tmp
    mv ${prefix}.fasta.tmp ${prefix}.fasta

    legsta ${params.legsta_options} \
      ${prefix}.fasta \
      --csv \
      | tee -a $log_file \
      | awk -v sample=${prefix} '{print sample "," $0 }' \
      | sed '0,/${prefix}/{s/${prefix}/sample/}' > legsta/${prefix}_legsta.csv
  """
}