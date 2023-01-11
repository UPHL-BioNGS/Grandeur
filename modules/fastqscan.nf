process fastqscan {
  tag           "${sample}"
  stageInMode   "copy"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/fastq-scan:1.0.1'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-large'
    
  input:
  tuple val(sample), file(fastq), val(size)

  output:
  path "fastq-scan/${sample}_fastqscan.json"                     , emit: json
  path "fastq-scan/${sample}_fastqscan.csv"                      , emit: collect
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log", emit: log_files

  shell:
  '''
    mkdir -p fastq-scan logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    fastq-scan -v >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    gzip --decompress !{fastq}

    cat *fastq | \
      fastq-scan !{params.fastqscan_options} \
      -g !{size} \
      | tee -a fastq-scan/!{sample}_fastqscan.json

    head -n 19 fastq-scan/!{sample}_fastqscan.json | tail -n 17 | awk '{print $1}' | sed 's/://g' | tr "\\n" "," | sed 's/,$/\\n/g' | awk '{print "sample," $0 }'                      >  fastq-scan/!{sample}_fastqscan.csv
    head -n 19 fastq-scan/!{sample}_fastqscan.json | tail -n 17 | awk '{print $2}' | sed 's/,//g' | tr "\\n" "," | sed 's/,$/\\n/g' | awk -v sample=!{sample} '{print sample "," $0 }' >> fastq-scan/!{sample}_fastqscan.csv
  '''
}