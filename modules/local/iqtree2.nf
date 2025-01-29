process IQTREE2 {
  tag           "Phylogenetic analysis"
  label         "process_high"
  container     'staphb/iqtree2:2.3.6'
  
  input:
  file(msa)

  output:
  path "iqtree2/iqtree*"                                               , emit: tree
  tuple val("iqtree"), file("iqtree2/iqtree.contree"), optional: true  , emit: newick
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log", emit: log
  path  "versions.yml"                                                 , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args     = task.ext.args ?: '-t RANDOM -m GTR+F+I -bb 1000 -alrt 1000'
  def outgroup = params.iqtree2_outgroup ? "-o ${params.iqtree2_outgroup}" : "" 
  """
    mkdir -p iqtree2 logs/${task.process}
    log_file=logs/${task.process}/${task.process}.${workflow.sessionId}.log

    iqtree2 ${args} \
      -s ${msa} \
      -pre iqtree2/iqtree \
      -nt AUTO \
      -ntmax ${task.cpus} \
      ${outgroup} \
      | tee -a \$log_file

    if [ -f "iqtree2/iqtree.treefile" ]; then cp iqtree2/iqtree.treefile iqtree2/iqtree.treefile.nwk ; fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree2: \$(echo \$(iqtree2 -version 2>&1) | sed 's/^IQ-TREE multicore version //;s/ .*//')
    END_VERSIONS
  """
}
