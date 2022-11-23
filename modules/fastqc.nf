process fastqc {
  tag "${sample}"
  
  input:
  tuple val(sample), file(raw)

  output:
  path "fastqc/*"                                                , emit: fastq_files
  path "fastqc/*_fastqc.zip"                                     , emit: for_multiqc
  tuple val(sample), env(raw_1)                                  , emit: fastqc_1_results
  tuple val(sample), env(raw_2)                                  , emit: fastqc_2_results
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log", emit: log_files

  shell:
  '''
    mkdir -p fastqc logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    fastqc --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    fastqc !{params.fastqc_options} \
      --outdir fastqc \
      --threads !{task.cpus} \
      !{raw} \
      | tee -a $log_file

    zipped_fastq=($(ls fastqc/*fastqc.zip) "")

    raw_1=$(unzip -p ${zipped_fastq[0]} */fastqc_data.txt | grep "Total Sequences" | awk '{ print $3 }' )
    raw_2=$(unzip -p fastqc/*fastqc.zip */fastqc_data.txt | grep "Total Sequences" | awk '{ print $3 }' )

    if [ -z "$raw_1" ] ; then raw_1="0" ; fi
    if [ -z "$raw_2" ] ; then raw_2="0" ; fi
  '''
}
