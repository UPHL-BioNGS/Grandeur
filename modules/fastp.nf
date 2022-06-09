process fastp {
  tag "${sample}"

  input:
  tuple val(sample), file(reads)

  output:
  tuple val(sample), file("fastp/${sample}_fastp_R{1,2}.fastq.gz"),      emit: fastq
  path "fastp/${sample}_fastp.html",                                     emit: html
  path "fastp/${sample}_fastp.json",                                     emit: fastp_files
  path "logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}", emit: log
  tuple val(sample), env(passed_reads),                                  emit: fastp_results

  shell:
  '''
    mkdir -p fastp logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    fastp --version >> $log_file 2>> $err_file
    echo "container : !{task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    fastp !{params.fastp_options} \
      -i !{reads[0]} \
      -I !{reads[1]} \
      -o fastp/!{sample}_fastp_R1.fastq.gz \
      -O fastp/!{sample}_fastp_R2.fastq.gz \
      -h fastp/!{sample}_fastp.html \
      -j fastp/!{sample}_fastp.json \
      2>> $err_file >> $log_file

    passed_reads=$(grep "reads passed filter" $err_file | tail -n 1 | cut -f 2 -d ":" | sed 's/ //g' )
    if [ -z "$passed_reads" ] ; then passed_reads="0" ; fi
  '''
}
