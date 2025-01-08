process PLASMIDFINDER {
  tag           "${meta.id}"
  label         "process_medium"
  container     'staphb/plasmidfinder:2.1.6_2024-03-07'


  input:
  tuple val(meta), file(file), file(script)

  output:
  path "plasmidfinder/*/*"                , emit: files
  path "plasmidfinder/*_plasmidfinder.tsv", emit: collect, optional: true
  path "logs/${task.process}/*.log"       , emit: log
  path "versions.yml"                     , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p plasmidfinder/${prefix} logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    plasmidfinder.py ${args} \
      -i ${file} \
      -o plasmidfinder/${prefix} \
      --extented_output \
      | tee -a \$log_file

    python3 ${script} plasmidfinder/${prefix}/results_tab.tsv plasmidfinder/${prefix}_plasmidfinder.tsv plasmidfinder ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plasmidfinder: "${task.container}"
    END_VERSIONS
  """
}
