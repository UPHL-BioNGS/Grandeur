process plasmidfinder {
  tag "${sample}"

  input:
  tuple val(sample), file(file)

  output:
  path "plasmidfinder/${sample}/*"                               , emit: files
  tuple val(sample), env(plasmids)                               , emit: plasmids
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log", emit: log

  shell:
  '''
    mkdir -p plasmidfinder/!{sample} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    echo "plasmidfinder.py: no version" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    plasmidfinder.py !{params.plasmidfinder_options} \
      -i !{file} \
      -o plasmidfinder/!{sample} \
      | tee -a $log_file

     plasmids=$(cat plasmidfinder/!{sample}/data.json | tr "," "\\n" | tr "{" "\\n" |  grep plasmid | awk '{print $2 }' | sort | uniq | grep \\" | tr "\\n" "," | sed 's/,$//g' | sed 's/\"//g' | sed 's/}//g' | sed 's/{//g')
     if [ -z "$plasmids" ] ; then plasmids="No hit found" ; fi
  '''
}
