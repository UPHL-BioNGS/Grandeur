process PBPTYPER {
  tag           "${meta.id}"
  label         "process_medium"
  container     'staphb/pbptyper:2.0.0'


  input:
  tuple val(meta), file(contigs)

  output:
  path "pbptyper/${meta.id}.tsv"   , emit: collect
  path "pbptyper/*"                , emit: all
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml"              , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p pbptyper logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log
    
    pbptyper ${args} \
      --input ${contigs} \
      --prefix ${prefix} \
      --outdir pbptyper \
      | tee -a \$log_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbptyper: \$(echo \$(pbptyper --version 2>&1) | sed 's/^.*pbptyper, version //;' )
    END_VERSIONS
  """
}
