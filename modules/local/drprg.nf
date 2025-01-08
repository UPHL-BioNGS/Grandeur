process DRPRG {
  tag           "${meta.id}"
  label         "process_medium"
  container     'staphb/drprg:0.1.1'

  input:
  tuple val(meta), file(contigs)

  output:
  tuple val(meta), file("drprg/*/*.drprg.json"), emit: json
  path "drprg/*/*", emit: results
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p drprg logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log
    
    drprg predict ${args} \
      -x /drprg/mtb/mtb \
      -i ${contigs} \
      -o drprg/${prefix} \
      --sample ${prefix} \
      | tee =a \$log_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drprg: \$( drprg --version | awk '{print \$NF}')
    END_VERSIONS
  """
}