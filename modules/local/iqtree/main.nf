process IQTREE {
  tag           "${msa.baseName}"
  label         "process_high"
  container     'staphb/iqtree3:3.1.3'
  
  input:
  file(msa)

  output:
  path "iqtree/iqtree*", emit: tree
  path "iqtree/*nwk", optional: true, emit: newick
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: '-t RANDOM -m GTR+F+I -bb 1000 -alrt 1000'
  def prefix = task.ext.prefix ?: "iqtree_${msa.baseName}"

  """
    mkdir -p iqtree logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${task.process}.${workflow.sessionId}.log

    iqtree3 ${args} \
      -s ${msa} \
      -pre iqtree/${prefix} \
      -nt AUTO \
      -ntmax ${task.cpus} \
      | tee -a \$log_file

    if [ -f "iqtree/${prefix}.treefile" ]; then cp iqtree/${prefix}.treefile iqtree/${prefix}.treefile.nwk ; fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree: \$(echo \$( iqtree3 --version | head -n 1 | awk '{print \$3}'))
    END_VERSIONS
  """
}
