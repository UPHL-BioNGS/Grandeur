process ELGATO {
  tag           "${meta.id}"
  label         "process_medium"
  container     'staphb/elgato:1.20.1'

  input:
  tuple val(meta), file(contigs)

  output:
  path "elgato/*/possible_mlsts.txt", emit: collect
  path "logs/${task.process}/*.log" , emit: log
  path "versions.yml"               , emit: versions
  tuple val(meta)                   , emit: meta

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p elgato logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    el_gato.py ${args} \
      --header \
      --assembly ${contigs} \
      --sample ${prefix} \
      --out elgato/${prefix} \
      --threads ${task.cpus} \
      | tee -a \$log_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      elgato: "${task.container}"
    END_VERSIONS
  """
}
