process plasmidfinder {
  tag "${sample}"

  when:
  params.fastq_processes =~ /plasmidfinder/ && task.process =~ /fastq_information/ || params.contig_processes =~ /plasmidfinder/ && task.process =~ /contig_information/

  input:
  tuple val(sample), file(file)

  output:
  path "plasmidfinder/${sample}/*"                                      , emit: files
  tuple val(sample), env(plasmids)                                      , emit: plasmids
  path "logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}" , emit: log

  shell:
  '''
    mkdir -p plasmidfinder/!{sample} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    echo "plasmidfinder.py: no version" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    plasmidfinder.py !{params.plasmidfinder_options} \
      -i !{file} \
      -o plasmidfinder/!{sample} \
      2>> $err_file >> $log_file

     plasmids=$(cat plasmidfinder/!{sample}/data.json | tr "," "\\n" | tr "{" "\\n" |  grep plasmid | awk '{print $2 }' | sort | uniq | grep \\" | tr "\\n" "," | sed 's/,$//g' | sed 's/\"//g' | sed 's/}//g' | sed 's/{//g')
     if [ -z "$plasmids" ] ; then plasmids="No hit found" ; fi
  '''
}
