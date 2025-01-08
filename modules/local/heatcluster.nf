process HEATCLUSTER {
  tag           "HeatCluster"
  label         "process_single"
  container     'staphb/heatcluster:1.0.2c'

  input:
  file(matrix)

  output:
  path "heatcluster/*", optional : true
  path "heatcluster/heatcluster.png", optional : true, emit: for_multiqc
  path "logs/${task.process}/*.log", emit: log_files
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: '-t png'
  """
    mkdir -p heatcluster logs/${task.process}
    log_file=logs/${task.process}/heatcluster.${workflow.sessionId}.log

    heatcluster.py ${args} \
        -i ${matrix} \
        -o heatcluster/heatcluster \
        | tee -a \$log_file

    if [ -f "sorted_matrix.csv" ]; then cp sorted_matrix.csv heatcluster/. ; fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      heatcluster: \$(echo \$(heatcluster.py --version | grep -v DeprecationWarning | grep -i heatcluster | awk '{print \$NF}' ))
    END_VERSIONS
  """
}
