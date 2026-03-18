process PLASMIDFINDER {
  tag           "${meta.id}"
  label         "process_low"
  container     'staphb/plasmidfinder:3.0.2'

  input:
  tuple val(meta), file(file)

  output:
  tuple val(meta), file("plasmidfinder/*/*"), emit: files
  path "plasmidfinder/*/*json", emit: collect, optional: true
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    export HOME=\$PWD
    git config --global --add safe.directory /database
    
    mkdir -p plasmidfinder/${prefix} logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    python -m plasmidfinder ${args} \
      -i ${file} \
      -o plasmidfinder/${prefix} \
      -j plasmidfinder/${prefix}/results_${prefix}_plasmidfinder.json \
      | tee -a \$log_file

    rm -rf plasmidfinder/${prefix}/tmp
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plasmidfinder: \$(echo \$(python -m plasmidfinder -v))
    END_VERSIONS
  """
}
