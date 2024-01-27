process mashtree {
  tag           "Phylogenetic analysis"
  label         "process_medium"
  label         "maxcpus"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/mashtree:1.4.6'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA cpus 14
  //#UPHLICA memory 60.GB
  //#UPHLICA time '24h'
  
  input:
  file(contigs)

  output:
  path "mashtree/*"                                                    , emit: tree
  tuple val("mashtree"), file("mashtree/mashtree.nwk"), optional: true , emit: newick
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log", emit: log
  path  "versions.yml"                          , emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
      def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p mashtree logs/${task.process}
    log_file=logs/${task.process}/${task.process}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    mashtree --version >> $log_file
    echo "container : ${task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    mashtree ${params.mashtree_options} \
      --numcpus ${task.cpus} \
      ${contigs} \
      --outtree mashtree/mashtree.nwk \
      | tee -a $log_file
  """
}
