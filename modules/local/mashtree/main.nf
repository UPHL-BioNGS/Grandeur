process MASHTREE {
  tag           "Phylogenetic analysis"
  label         "process_medium"
  container     'staphb/mashtree:1.4.6'
  
  input:
  file(assemblies)

  output:
  path "mashtree/*"                                                    , emit: tree
  tuple val("mashtree"), file("mashtree/mashtree.nwk"), optional: true , emit: newick
  path "logs/${task.process}/*.log"                                    , emit: log
  path "versions.yml"                                                  , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def prefix = task.ext.prefix ?: "mashtree"
  def args   = task.ext.args   ?: "--outmatrix mashtree/${prefix}.txt"
  def input  = assemblies.join(" ")
  """
    mkdir -p mashtree logs/${task.process}
    log_file=logs/${task.process}/${task.process}.${workflow.sessionId}.log

    mashtree ${args} \
      --numcpus ${task.cpus} \
      ${input} \
      --outtree mashtree/${prefix}.nwk \
      | tee -a \$log_file
      
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      mashtree: \$( echo \$( mashtree --version 2>&1 ) | sed 's/^.*Mashtree //' )
    END_VERSIONS
  """
}
