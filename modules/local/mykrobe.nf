process MYKROBE {
  tag           "${meta.id}"
  label         "process_medium"
  container     'staphb/mykrobe:0.13.0'

  input:
  tuple val(meta), file(contigs)

  output:
  path "mykrobe/*.csv", optional: true, emit: collect
  path "mykrobe/*.json", optional: true, emit: json
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: '--species tb --format json_and_csv'
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p tmp/ logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    mykrobe predict ${args} \
      --sample ${prefix} \
      --output mykrobe/${prefix} \
      --seq ${contigs} \
      --threads ${task.cpus} \
      | tee -a \$log_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mykrobe: \$(echo \$(mykrobe --version 2>&1) | sed 's/^.*mykrobe v//' )
    END_VERSIONS
  """
}