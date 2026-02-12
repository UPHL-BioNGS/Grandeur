process GOTREE {
  tag           "${analysis}"
  label         "process_medium"
  container     'staphb/gotree:latest'

  
  input:
  tuple val(analysis), file(newick)

  output:
  path "gotree/${analysis}_tree.png", emit: for_multiqc
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${analysis}"
  """
    mkdir -p gotree logs/${task.process}
    log_file=logs/${task.process}/${analysis}.${task.process}.${workflow.sessionId}.log

    gotree ${args} \
      -i ${newick} \
      -o gotree/${prefix}_tree.png \
      | tee -a \$log_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      gotree: \$(gotree --version)
    END_VERSIONS
  """
}