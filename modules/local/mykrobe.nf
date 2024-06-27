process mykrobe {
  tag           "${meta.id}"
  label         "process_medium"
  stageInMode   "copy"
  publishDir    path: params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/mykrobe:0.13.0'
  time          '1h'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}

  input:
  tuple val(meta), file(contigs)

  output:
  path "mykrobe/*.csv", emit: collect
  path "mykrobe/*.json", emit: json
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  (task.ext.when == null || task.ext.when)

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p tmp/ logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    mykrobe predict ${args} \
      --species tb \
      --sample ${prefix} \
      --output mykrobe/${prefix} \
      --seq ${contigs} \
      --threads ${task.cpus} \
      --format json_and_csv \
      | tee -a \$log_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mykrobe: \$(echo \$(mykrobe --version 2>&1) | sed 's/^.*mykrobe v//' )
    END_VERSIONS
  """
}