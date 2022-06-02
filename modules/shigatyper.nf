process shigatyper {
  tag "${sample}"
  label "medcpus"

  when:
  params.fastq_processes =~ /shigatyper/ && flag =~ 'found'

  input:
  tuple val(sample), file(fastq), val(flag)

  output:
  path "shigatyper/${sample}_shigatyper.tsv"                           , emit: files
  path "shigatyper/${sample}-hits.tsv", optional: true                 , emit: hits
  tuple val(sample), env(predictions)                                  , emit: predictions
  tuple val(sample), env(lacy_cada)                                    , emit: cada
  path "logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}", emit: log

  shell:
  '''
    mkdir -p shigatyper logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    shigatyper --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    shigatyper !{params.shigatyper_options} \
      --R1 !{fastq[0]} \
      --R2 !{fastq[1]} \
      2>> $err_file \
      > shigatyper/!{sample}_shigatyper.tsv

    if [ -f "!{sample}.tsv" ]; then cp !{sample}.tsv shigatyper/!{sample}-hits.tsv ; fi

    predictions=$(grep -v "prediction" shigatyper/!{sample}_shigatyper.tsv | grep -wv "Hit" | cut -f 2 | tr '\\n' ',' | sed 's/,$//g' )
    lacy_cada="$(grep -ie "lac" -ie "cad" $err_file | head -n 1)"
    if [ -z "$predictions" ] ; then predictions='none' ; fi
    if [ -z "$lacy_cada" ] ; then lacy_cada='none' ; fi
  '''
}
