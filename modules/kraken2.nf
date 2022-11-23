process kraken2_fastq {
  tag "${sample}"
  label "maxcpus"

  input:
  tuple val(sample), file(file), path(kraken2_db)

  output:
  path "kraken2/${sample}_kraken2_report.txt"                    , emit: for_multiqc
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log", emit: log
  tuple val(sample), env(top_hit)                                , emit: top_hit
  tuple val(sample), env(top_perc)                               , emit: top_perc
  tuple val(sample), env(top_reads)                              , emit: top_reads

  shell:
  '''
    mkdir -p kraken2 logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    date > $log_file
    echo "container : !{task.container}" >> $log_file
    kraken2 --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    kraken2 !{params.kraken2_options} \
      --paired \
      --classified-out cseqs#.fq \
      --threads !{task.cpus} \
      --db !{kraken2_db} \
      !{file} \
      --report kraken2/!{sample}_kraken2_report.txt \
      | tee -a $log_file

    top_hit=$(cat kraken2/!{sample}_kraken2_report.txt   | grep -w S | sort | tail -n 1 | awk '{print $6 " " $7}')
    top_perc=$(cat kraken2/!{sample}_kraken2_report.txt  | grep -w S | sort | tail -n 1 | awk '{print $1}')
    top_reads=$(cat kraken2/!{sample}_kraken2_report.txt | grep -w S | sort | tail -n 1 | awk '{print $2}')
    if [ -z "$top_hit" ] ; then top_hit="NA" ; fi
    if [ -z "$top_perc" ] ; then top_perc="0" ; fi
    if [ -z "$top_reads" ] ; then top_reads="0" ; fi
  '''
}

process kraken2_fasta {
  tag "${sample}"
  label "maxcpus"

  input:
  tuple val(sample), file(file), path(kraken2_db)

  output:
  path "kraken2/${sample}_kraken2_report_contigs.txt"            , emit: for_multiqc
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log", emit: log
  tuple val(sample), env(top_hit)                                , emit: top_hit
  tuple val(sample), env(top_perc)                               , emit: top_perc
  tuple val(sample), env(top_reads)                              , emit: top_reads

  shell:
    '''
    mkdir -p kraken2 logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    date > $log_file
    echo "container : !{task.container}" >> $log_file
    kraken2 --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    kraken2 !{params.kraken2_options} \
      --threads !{task.cpus} \
      --db !{kraken2_db} \
      !{file} \
      --report kraken2/!{sample}_kraken2_report_contigs.txt \
      | tee -a $log_file

    top_hit=$(cat kraken2/!{sample}_kraken2_report_contigs.txt   | grep -w S | sort | tail -n 1 | awk '{print $6 " " $7}')
    top_perc=$(cat kraken2/!{sample}_kraken2_report_contigs.txt  | grep -w S | sort | tail -n 1 | awk '{print $1}')
    top_reads=$(cat kraken2/!{sample}_kraken2_report_contigs.txt | grep -w S | sort | tail -n 1 | awk '{print $2}')
    if [ -z "$top_hit" ]   ; then top_hit="NA"  ; fi
    if [ -z "$top_perc" ]  ; then top_perc="0"  ; fi
    if [ -z "$top_reads" ] ; then top_reads="0" ; fi
  '''
}
