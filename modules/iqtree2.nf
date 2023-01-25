process iqtree2 {
  tag           "Pylogenetic Analysis"
  label         "maxcpus"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/iqtree2:2.1.2'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-large'
  //#UPHLICA cpus 6
  //#UPHLICA memory 26.GB
  
  input:
  file(msa)

  output:
  path "iqtree2/iqtree*"                                               , emit: tree
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log", emit: log

  shell:
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
