process mlst {
  tag           "$meta.id"
  label         "process_medium"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/mlst:2.23.0-2024-01'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '10m'

  input:
  tuple val(meta), file(contig), file(script)

  output:
  path "mlst/${prefix}_mlst.tsv"                                 , emit: collect
  path "logs/${task.process}/${prefix}.${workflow.sessionId}.log", emit: log
  path  "versions.yml"                          , emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
      def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p mlst logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : ${task.container}" >> $log_file
    mlst --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    mlst ${params.mlst_options} \
      --threads ${task.cpus} \
      ${contig} \
      > ${prefix}_mlst.txt

    python3 ${script} ${prefix}_mlst.txt mlst/${prefix}_mlst.tsv mlst ${prefix}
  """
}
