process drprg {
  tag           "${meta.id}"
  label         "process_medium"
  stageInMode   "copy"
  publishDir    path: params.outdir, mode: 'copy'
  container     'staphb/drprg:0.1.1'
  time          '10m'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}

  input:
  tuple val(meta), file(contigs), val(flag)

  output:
  tuple val(meta), val("drprg"), file("drprg/*/*.drprg.json"), emit: json
  path "drprg/*/*", emit: results
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  (task.ext.when == null || task.ext.when) && flag =~ 'found'

  shell:
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
        drprg: \$( drprg --version )
    END_VERSIONS
  """
}