process PHYTREEVIZ {
  tag           "${analysis}"
  label         "process_medium"
  container     'staphb/phytreeviz:0.2.0'

  
  input:
  tuple val(analysis), file(newick)

  output:
  path "phytreeviz/${analysis}_tree.png", emit: for_multiqc
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${analysis}"
  """
    mkdir -p phytreeviz logs/${task.process}
    log_file=logs/${task.process}/${analysis}.${task.process}.${workflow.sessionId}.log

    phytreeviz ${args} \
      -i ${newick} \
      -o phytreeviz/${prefix}_tree.png \
      | tee -a \$log_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      phytreeviz: \$(phytreeviz --version)
    END_VERSIONS
  """
}
