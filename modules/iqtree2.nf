process iqtree2 {
  tag "Pylogenetic Analysis"
  label "maxcpus"

  when:
  params.phylogenetic_processes =~ /iqtree2/

  input:
  file(msa)

  output:
  path "iqtree2/iqtree*"                                                     , emit: tree
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.{log,err}", emit: log

  shell:
  '''
    mkdir -p iqtree2 logs/!{task.process}
    log_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    iqtree2 -v >> $log_file
    echo "container : !{task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    outgroup=''
    if [ -n "!{params.outgroup}" ] ; then outgroup="-o !{params.outgroup}" ; fi

    iqtree2 !{params.iqtree2_options} \
      -s !{msa} \
      -pre iqtree2/iqtree \
      -nt AUTO \
      -ntmax !{task.cpus} \
      $outgroup \
      2>> $err_file >> $log_file
  '''
}
