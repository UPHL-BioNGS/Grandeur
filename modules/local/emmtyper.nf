process EMMTYPER {
  tag           "${meta.id}"
  label         "process_medium"
  container     'staphb/emmtyper:0.2.0'

  input:
  tuple val(meta), file(contigs), file(script)

  output:
  path "emmtyper/*_emmtyper.txt"   , emit: collect
  path "emmtyper/*"                , emit: everything
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml"              , emit: versions
  tuple val(meta)                  , emit: meta

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p emmtyper logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    PATH=/opt/conda/envs/emmtyper/bin:\$PATH
    export PATH
    
    emmtyper ${args} \
      --output-format 'verbose' \
      ${contigs} \
      | tee -a \$log_file \
      > ${prefix}_emmtyper.txt

    python3 ${script} ${prefix}_emmtyper.txt emmtyper/${prefix}_emmtyper.txt emmtyper ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      emmtyper: \$( echo \$(emmtyper --version 2>&1) | sed 's/^.*emmtyper v//' )
    END_VERSIONS
  """
}
