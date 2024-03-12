process pbptyper {
  tag           "${meta.id}"
  label         "process_medium"
  stageInMode   "copy"
  publishDir    path: params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/pbptyper:1.0.4'
  time          '1h'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}

  input:
  tuple val(meta), file(contigs), val(flag)

  output:
  path "pbptyper/${meta.id}.tsv"   , emit: collect
  path "pbptyper/*"                , emit: all
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml"              , emit: versions

  when:
  (task.ext.when == null || task.ext.when) && flag =~ 'found'

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p pbptyper logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log
    
    pbptyper ${args} \
      --assembly ${contigs} \
      --prefix ${prefix} \
      --outdir pbptyper \
      | tee -a \$log_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbptyper: \$(echo \$(pbptyper --version 2>&1) | sed 's/^.*pbptyper, version //;' )
    END_VERSIONS
  """
}