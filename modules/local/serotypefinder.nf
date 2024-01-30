process serotypefinder {
  tag           "${meta.id}"
  label         "process_medium"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/serotypefinder:2.0.1'
  maxForks      10
  time          '10m'
  //errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}  

  input:
  tuple val(meta), file(file), val(flag), file(script)

  output:
  path "serotypefinder/*/*"                 , emit: files
  path "serotypefinder/*_serotypefinder.tsv", emit: collect, optional: true
  path "logs/${task.process}/*.log"         , emit: log
  path "versions.yml"                       , emit: versions

  when:
  (task.ext.when == null || task.ext.when) && flag =~ 'found'

  shell:
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      serotypefinder.py: ${task.container}
    END_VERSIONS
  """
}
