process lyveset_shuffle {
  tag "${sample}"

  when:
  params.fastq_processes =~ /cg_pipeline/

  input:
  tuple val(sample), file(reads)

  output:
  tuple val(sample), file("shuffled/${sample}_shuffled.fastq.gz")      , emit: shuffled
  path "logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}", emit: log

  shell:
  '''
    mkdir -p shuffled logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    run_assembly_shuffleReads.pl -gz !{reads} 2>> $err_file > shuffled/!{sample}_shuffled.fastq.gz
  '''
}

process lyveset_cg_pipeline {
  tag "${sample}"
  label "medcpus"

  input:
  tuple val(sample), file(fastq), val(mash), val(genus), val(species), file(genome_file)

  output:
  path "cg_pipeline/${sample}_cg_pipeline_report.txt", optional: true  , emit: collect
  tuple val(sample), env(read_length)                                  , emit: read_length
  tuple val(sample), env(quality)                                      , emit: quality
  tuple val(sample), env(coverage)                                     , emit: coverage
  tuple val(sample), env(reference_genome_length)                      , emit: ref_genome_length
  path "logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}", emit: log

  shell:
  '''
    mkdir -p cg_pipeline logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    genome_length=''
    if [ "!{genus}" != "null" ] && [ "!{species}" != "null" ] ; then genome_length=$(grep !{genus} !{genome_file} | grep !{species} | grep -v "#" | head -n 1 | cut -f 2 -d ":" | cut -f 1 -d "," | awk '{ print $0 "e+06" }') ; fi
    if [ -z "$genome_length" ] && [ "!{mash}" != "null" ] ; then genome_length=$(echo !{mash} | xargs printf "%.0f" ) ; fi

    if [ -n "$genome_length" ]
    then
      run_assembly_readMetrics.pl !{fastq} \
        !{params.cg_pipeline_options} \
        --fast \
        --numcpus !{task.cpus} \
        -e $genome_length \
        2>> $err_file > cg_pipeline/!{sample}_cg_pipeline_report.txt

        read_length=$(cut -f 2 cg_pipeline/!{sample}_cg_pipeline_report.txt | tail -n 1 )
        quality=$(cut -f 6 cg_pipeline/!{sample}_cg_pipeline_report.txt | tail -n 1 )
        coverage=$(cut -f 9 cg_pipeline/!{sample}_cg_pipeline_report.txt | tail -n 1 )
    else
      genome_length='0'
      read_length='NA'
      quality='NA'
      coverage='NA'
      echo "Could not determine genome length of isolate, so could not run GC pipeline" | tee $log_file
    fi

    if [ -z "$read_length" ] ; then read_length='NA' ; fi
    if [ -z "$quality" ] ; then quality='NA' ; fi
    if [ -z "$coverage" ] ; then coverage='NA' ; fi
    reference_genome_length=$genome_length
  '''
}
