process bbduk{
  tag "${sample}"

  input:
  tuple val(sample), file(reads)

  output:
  tuple val(sample), file("bbduk/${sample}_rmphix_R{1,2}.fastq.gz"),     emit: fastq
  path "bbduk/*",                                                        emit: files
  path "bbduk/${sample}.phix.stats.txt",                                 emit: stats
  path "logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}", emit: log
  tuple val(sample), env(phix_reads),                                    emit: phix_reads

  shell:
  '''
    mkdir -p bbduk logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    bbduk.sh --version >> $log_file 2>> $err_file
    echo "container : !{task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    bbduk.sh !{params.bbduk_options} \
      in1=!{reads[0]} \
      in2=!{reads[1]} \
      out1=bbduk/!{sample}_rmphix_R1.fastq.gz \
      out2=bbduk/!{sample}_rmphix_R2.fastq.gz \
      outm=bbduk/!{sample}.matched_phix.fq \
      ref=/opt/bbmap/resources/phix174_ill.ref.fa.gz \
      stats=bbduk/!{sample}.phix.stats.txt \
      2>> $err_file | tee -a $log_file

    phix_reads=$(grep Matched bbduk/!{sample}.phix.stats.txt | cut -f 2)
  '''
}
