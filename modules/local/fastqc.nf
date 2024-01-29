process fastqc {
  tag           "$meta.id"
  label         "process_single"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/fastqc:0.12.1'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '10m'
    
  input:
  tuple val(meta), file(fastq)

  output:
  path "fastqc/*html"                                            , emit: fastq_files
  path "fastqc/*_fastqc.zip"                                     , emit: for_multiqc
  path "fastqc/${prefix}_summary.csv"                            , emit: collect
  path "logs/${task.process}/${prefix}.${workflow.sessionId}.log", emit: log_files
  path  "versions.yml"                          , emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p fastqc logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : ${task.container}" >> $log_file
    fastqc --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    fastqc ${params.fastqc_options} \
      --outdir fastqc \
      --threads ${task.cpus} \
      ${fastq} \
      --extract \
      | tee -a $log_file

    header=$(head -n 10 fastqc/*_fastqc/fastqc_data.txt | cut -f 1 | tr "\\n" ",")

    for data in fastqc/*_fastqc/fastqc_data.txt
    do
      if [ $ -f "fastqc/${prefix}_summary.csv" ]
      then
        head -n 10 $data | cut -f 1 | tr "\\n" "," | sed 's/,$/\\n/' | sed 's/#//g' | sed 's/>//g' | awk '{print "sample," $0}' > fastqc/${prefix}_summary.csv
      fi

      head -n 10 $data | cut -f 2 | tr "\\n" "," | sed 's/,$/\\n/' | awk -v sample=${prefix} '{print sample "," $0}' >> fastqc/${prefix}_summary.csv
    done
  """
}