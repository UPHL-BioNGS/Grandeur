process SEROTYPEFINDER {
  tag           "${meta.id}"
  label         "process_medium"
  container     'staphb/serotypefinder:2.0.2'
 

  input:
  tuple val(meta), file(file), file(script)

  output:
  path "serotypefinder/*/*"                 , emit: files
  path "serotypefinder/*_serotypefinder.tsv", emit: collect, optional: true
  path "logs/${task.process}/*.log"         , emit: log
  path "versions.yml"                       , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p serotypefinder/${prefix} logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    serotypefinder.py ${args} \
      -i ${file} \
      -o serotypefinder/${prefix} \
      -x \
      | tee -a \$log_file

    cp serotypefinder/${prefix}/results_tab.tsv serotypefinder/${prefix}_serotypefinder.tsv

    python3 ${script} serotypefinder/${prefix}/results_tab.tsv serotypefinder/${prefix}_serotypefinder.tsv serotypefinder ${prefix}

    rm -rf serotypefinder/${prefix}/tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      serotypefinder.py: ${task.container}
    END_VERSIONS
  """
}
