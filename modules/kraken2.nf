process kraken2_fastq {
  tag "${sample}"
  label "maxcpus"

  input:
  tuple val(sample), file(file), path(kraken2_db)

  output:
  path "kraken2/${sample}_kraken2_report_reads.txt"              , emit: for_multiqc
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log", emit: log
  path "kraken2/${sample}_reads_summary.csv"                     , emit: results

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
      --report kraken2/!{sample}_kraken2_report_reads.txt \
      | tee -a $log_file

    echo "Sample,Type,Percentage of fragments,Number of fragments,Number of fragments assigned directly to this taxon,Rank code,NCBI taxonomic ID number,Scientific name" > kraken2/!{sample}_reads_summary.csv
    cat kraken2/!{sample}_kraken2_report_reads.txt | grep -w S | \
      awk -v sample=!{sample} '{ if ($1 >= 5 ) print sample ",reads," $1 "," $2 "," $3 "," $4 "," $5 "," $6 "_" $7 }' | \
      sort >> kraken2/!{sample}_reads_summary.csv
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
  path "kraken2/${sample}_contigs_summary.csv"                   , emit: results

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

    echo "Sample,Type,Percentage of fragments,Number of fragments,Number of fragments assigned directly to this taxon,Rank code,NCBI taxonomic ID number,Scientific name" > kraken2/!{sample}_contigs_summary.csv
    cat kraken2/!{sample}_kraken2_report_contigs.txt | grep -w S | \
      awk -v sample=!{sample} '{ if ($1 >= 5 ) print sample ",contigs," $1 "," $2 "," $3 "," $4 "," $5 "," $6 "_" $7 }' | \
      sort >> kraken2/!{sample}_contigs_summary.csv
  '''
}
