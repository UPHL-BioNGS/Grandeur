process serotypefinder {
  tag "${sample}"
  label "medcpus"

  when:
  params.contig_processes =~ /serotypefinder/ && flag =~ 'found'

  input:
  tuple val(sample), file(file), val(flag)

  output:
  path "serotypefinder/${sample}/*"                                    , emit: files
  tuple val(sample), env(o_type)                                        , emit: ogroup
  tuple val(sample), env(h_type)                                        , emit: hgroup
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log" , emit:  log

  shell:
  '''
    mkdir -p serotypefinder/!{sample} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    serotypefinder.py !{params.serotypefinder_options} \
      -i !{file} \
      -o serotypefinder/!{sample} \
      -x \
      | tee -a $log_file

    h_type=$(cut -f 3 serotypefinder/!{sample}/results_tab.tsv | grep ^H | sort | uniq | tr '\\n' ',' | sed 's/,$//g' )
    o_type=$(cut -f 3 serotypefinder/!{sample}/results_tab.tsv | grep ^O | sort | uniq | tr '\\n' ',' | sed 's/,$//g' )
    if [ -z "$h_type" ] ; then h_type="none" ; fi
    if [ -z "$o_type" ] ; then o_type="none" ; fi
  '''
}
