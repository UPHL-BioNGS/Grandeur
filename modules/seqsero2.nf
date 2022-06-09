process seqsero2_fastq {
  tag "${sample}"
  label "medcpus"

  when:
  params.fastq_processes =~ /seqsero2/ && flag =~ 'found'

  input:
  tuple val(sample), file(file), val(flag)

  output:
  tuple val(sample), env(antigenic_profile)                             , emit: profile
  tuple val(sample), env(serotype)                                      , emit: serotype
  tuple val(sample), env(contamination)                                 , emit: contamination
  path "seqsero2/${sample}/*"                                           , emit: files
  path "seqsero2/${sample}/SeqSero_result.tsv"                          , emit: collect
  path "logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}" , emit: log

  shell:
  '''
    mkdir -p seqsero2 logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    SeqSero2_package.py --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    SeqSero2_package.py !{params.seqsero2_options} \
      -t 2 \
      -i !{file} \
      -p !{task.cpus} \
      -d seqsero2/!{sample} \
      -n !{sample} \
      2>> $err_file >> $log_file

    serotype=$(cut -f 9 seqsero2/!{sample}/SeqSero_result.tsv | tail -n 1)
    contamination=$(cut -f 10 seqsero2/!{sample}/SeqSero_result.tsv | tail -n 1)
    antigenic_profile=$(cut -f 8 seqsero2/!{sample}/SeqSero_result.tsv | tail -n 1)
    enteritidis_check=$(grep "Enteritidis" seqsero2/!{sample}/SeqSero_result.tsv | head -n 1 )
    sdf_check=$(grep "Detected Sdf" seqsero2/!{sample}/SeqSero_result.tsv | head -n 1 )

    if [ -n "$sdf_check" ] && [ -n "$enteritidis_check" ]
    then
      serotype="$serotype (Sdf+)"
    elif [ -z "$sdf_check" ] && [ -n "$enteritidis_check" ]
    then
      serotype="$serotype (Sdf-)"
    fi

    if [ -z "$serotype" ] ; then serotype='NA' ; fi
    if [ -z "$contamination" ] ; then contamination='NA' ; fi
    if [ -z "$antigenic_profile" ] ; then antigenic_profile='NA' ; fi
  '''
}

process seqsero2_fasta {
  tag "${sample}"
  label "medcpus"

  when:
  params.contig_processes =~ /seqsero2/ && flag =~ 'found'

  input:
  tuple val(sample), file(file), val(flag)

  output:
  tuple val(sample), env(antigenic_profile)                             , emit: profile
  tuple val(sample), env(serotype)                                      , emit: serotype
  tuple val(sample), env(contamination)                                 , emit: contamination
  path "seqsero2/${sample}/*"                                           , emit: files
  path "seqsero2/${sample}/SeqSero_result.tsv"                          , emit: collect
  path "logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}" , emit: log

  shell:
  '''
    mkdir -p seqsero2 logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    SeqSero2_package.py --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    SeqSero2_package.py !{params.seqsero2_options} \
      -m k \
      -t 4 \
      -i !{file} \
      -p !{task.cpus} \
      -d seqsero2/!{sample} \
      -n !{sample} \
      2>> $err_file >> $log_file

    serotype=$(cut -f 9 seqsero2/!{sample}/SeqSero_result.tsv | tail -n 1)
    contamination=$(cut -f 10 seqsero2/!{sample}/SeqSero_result.tsv | tail -n 1)
    antigenic_profile=$(cut -f 8 seqsero2/!{sample}/SeqSero_result.tsv | tail -n 1)
    enteritidis_check=$(grep "Enteritidis" seqsero2/!{sample}/SeqSero_result.tsv | head -n 1 )
    sdf_check=$(grep "Detected Sdf" seqsero2/!{sample}/SeqSero_result.tsv | head -n 1 )

    if [ -n "$sdf_check" ] && [ -n "$enteritidis_check" ]
    then
      serotype="$serotype (Sdf+)"
    elif [ -z "$sdf_check" ] && [ -n "$enteritidis_check" ]
    then
      serotype="$serotype (Sdf-)"
    fi

    if [ -z "$serotype" ] ; then serotype='NA' ; fi
    if [ -z "$contamination" ] ; then contamination='NA' ; fi
    if [ -z "$antigenic_profile" ] ; then antigenic_profile='NA' ; fi
  '''
}
