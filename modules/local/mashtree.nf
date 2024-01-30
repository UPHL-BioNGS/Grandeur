process mashtree {
  tag           "Phylogenetic analysis"
  label         "process_medium"
  label         "maxcpus"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/mashtree:1.4.6'
  time      '4h'
  //errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  
  input:
  file(contigs)

  output:
  path "mashtree/*"                                                    , emit: tree
  tuple val("mashtree"), file("mashtree/mashtree.nwk"), optional: true , emit: newick
  path "logs/${task.process}/*.log"                                    , emit: log
  path "versions.yml"                                                  , emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p mashtree logs/${task.process}
    log_file=logs/${task.process}/${task.process}.${workflow.sessionId}.log

    mashtree ${args} \
      --numcpus ${task.cpus} \
      ${contigs} \
      --outtree mashtree/mashtree.nwk \
      | tee -a \$log_file
      
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mashtree: \$( echo \$( mashtree --version 2>&1 ) | sed 's/^.*Mashtree //' )
    END_VERSIONS
  """
}
