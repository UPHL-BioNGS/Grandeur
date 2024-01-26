process iqtree2 {
  tag           "Phylogenetic analysis"
  label         "process_high"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/iqtree2:2.2.2.6'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA cpus 14
  //#UPHLICA memory 60.GB
  //#UPHLICA time '24h'
  
  input:
  file(msa)

  output:
  path "iqtree2/iqtree*"                                               , emit: tree
  tuple val("iqtree"), file("iqtree2/iqtree.contree"), optional: true  , emit: newick
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log", emit: log
  path  "versions.yml"                          , emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
      def args = task.ext.args ?: '-t RANDOM -m GTR+F+I -bb 1000 -alrt 1000'
    def prefix = task.ext.prefix ?: "${meta.id}"
  '''
    mkdir -p iqtree2 logs/!{task.process}
    log_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    iqtree2 -v >> $log_file
    echo "container : !{task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    outgroup=''
    if [ -n "!{params.iqtree2_outgroup}" ] ; then outgroup="-o !{params.iqtree2_outgroup}" ; fi

    iqtree2 !{params.iqtree2_options} \
      -s !{msa} \
      -pre iqtree2/iqtree \
      -nt AUTO \
      -ntmax !{task.cpus} \
      $outgroup \
      | tee -a $log_file

    if [ -f "iqtree2/iqtree.treefile" ]; then cp iqtree2/iqtree.treefile iqtree2/iqtree.treefile.nwk ; fi
  '''
}
