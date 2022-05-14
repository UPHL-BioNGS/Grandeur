process seqyclean {
  tag "${sample}"

  when:
  params.fastq_processes =~ /seqyclean/

  input:
  tuple val(sample), file(reads)

  output:
  tuple val(sample), file("seqyclean/${sample}_clean_PE{1,2}.fastq*")  , emit: clean_reads
  path "seqyclean/${sample}*"                                          , emit: collect
  path "seqyclean/${sample}_clean_SummaryStatistics.tsv"               , emit: for_multiqc
  path "logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}", emit: log
  tuple val(sample), env(perc_kept)                                    , emit: perc_kept
  tuple val(sample), env(kept)                                         , emit: pairskept

  shell:
  '''
    mkdir -p seqyclean logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    echo "seqyclean version: $(seqyclean -h | grep Version | head -n 1)" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    seqyclean !{params.seqyclean_options} \
      -c !{params.seqyclean_contaminant_file} \
      -1 !{reads[0]} \
      -2 !{reads[1]} \
      -o seqyclean/!{sample}_clean \
      -gz \
      2>> $err_file >> $log_file

    kept=$(cut -f 58 seqyclean/!{sample}_clean_SummaryStatistics.tsv | grep -v "Kept" | head -n 1)
    perc_kept=$(cut -f 59 seqyclean/!{sample}_clean_SummaryStatistics.tsv | grep -v "Kept" | head -n 1)
    if [ -z "$kept" ] ; then kept="0" ; fi
    if [ -z "$perc_kept" ] ; then perc_kept="0" ; fi
  '''
}
