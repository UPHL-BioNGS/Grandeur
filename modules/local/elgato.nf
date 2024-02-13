process elgato {
  tag           "${meta.id}"
  label         "process_medium"
  publishDir    path: params.outdir, mode: 'copy'
  container     'staphb/elgato:1.15.2'
  time          '10m'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}

  input:
  tuple val(meta), file(contigs), val(flag)

  output:
  path "elgato/*/possible_mlsts.txt", emit: collect
  path "logs/${task.process}/*.log" , emit: log
  path "versions.yml"               , emit: versions

  when:
  (task.ext.when == null || task.ext.when) && flag =~ 'found'

  shell:
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