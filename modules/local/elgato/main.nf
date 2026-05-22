process ELGATO {
  tag           "${meta.id}"
  label         "process_medium"
  container     'staphb/elgato:1.22.0'

  input:
  tuple val(meta), file(contigs)

  output:
  path "elgato/*/*_possible_mlsts.txt", emit: collect, optional: true
  path "elgato/*/*", emit: results, optional: true
  path "logs/${task.process}/*.log" , emit: log
  path "versions.yml"               , emit: versions
  val meta                          , emit: meta

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

    if [ -f "elgato/${prefix}/possible_mlsts.txt" ]
    then
      cp elgato/${prefix}/possible_mlsts.txt elgato/${prefix}/${prefix}_possible_mlsts.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      elgato: "${task.container}"
    END_VERSIONS
  """
}
