process heatcluster {
  tag           "HeatCluster"
  label         "process_single"
  publishDir    params.outdir, mode: 'copy'
  container     'quay.io/uphl/heatcluster:1.0.2-2023-12-19'
  time          '10m'
  //errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}

  input:
  file(matrix)

  output:
  path "heatcluster/heatcluster*"   , optional : true
  path "heatcluster/heatcluster.png", optional : true, emit: for_multiqc
  path "logs/${task.process}/*.log", emit: log_files
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: '-t png'
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p heatcluster logs/${task.process}
    log_file=logs/${task.process}/heatcluster.${workflow.sessionId}.log

    heatcluster.py ${args} \
        -i ${matrix} \
        -o heatcluster/heatcluster \
        | tee -a \$log_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        heatcluster: \$(echo \$(heatcluster.py -v 2>&1))
    END_VERSIONS
  """
}
