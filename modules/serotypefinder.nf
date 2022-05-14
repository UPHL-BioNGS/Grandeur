process serotypefinder {
  tag "${sample}"
  label "medcpus"

  when:
  ((params.fastq_processes =~ /serotypefinder/ && type == 'fastq' ) || (params.contig_processes =~ /serotypefinder/ && type == 'fasta' )) && flag =~ 'found'

  input:
  tuple val(sample), file(file), val(type), val(flag)

  output:
  path "${task.process}/${sample}/*"                                    , emit: files
  tuple val(sample), env(o_type)                                        , emit: ogroup
  tuple val(sample), env(h_type)                                        , emit: hgroup
  path "logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}" , emit:  log

  shell:
  if ( type == "fastq" )
  '''
    mkdir -p !{task.process}/!{sample} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    serotypefinder.py !{params.serotypefinder_options} \
      -i !{fasta} \
      -o !{task.process}/!{sample} \
      -x \
      2>> $err_file >> $log_file

    h_type=$(cut -f 3 !{task.process}/!{sample}/results_tab.tsv | grep ^H | sort | uniq | tr '\\n' ',' | sed 's/,$//g' )
    o_type=$(cut -f 3 !{task.process}/!{sample}/results_tab.tsv | grep ^O | sort | uniq | tr '\\n' ',' | sed 's/,$//g' )
    if [ -z "$h_type" ] ; then h_type="none" ; fi
    if [ -z "$o_type" ] ; then o_type="none" ; fi
  '''
  else if ( type == "fasta" )
  '''
    mkdir -p !{task.process}/!{sample} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    serotypefinder.py !{params.serotypefinder_options} \
      -i !{fasta} \
      -o !{task.process}/!{sample} \
      -x \
      2>> $err_file >> $log_file

    h_type=$(cut -f 3 !{task.process}/!{sample}/results_tab.tsv | grep ^H | sort | uniq | tr '\\n' ',' | sed 's/,$//g' )
    o_type=$(cut -f 3 !{task.process}/!{sample}/results_tab.tsv | grep ^O | sort | uniq | tr '\\n' ',' | sed 's/,$//g' )
    if [ -z "$h_type" ] ; then h_type="none" ; fi
    if [ -z "$o_type" ] ; then o_type="none" ; fi
  '''
}
