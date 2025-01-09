process FASTQC {
  tag           "${meta.id}"
  label         "process_single"
  container     'staphb/fastqc:0.12.1'
    
  input:
  tuple val(meta), file(fastq)

  output:
  path "fastqc/*html"              , emit: fastq_files
  path "fastqc/*_fastqc.zip"       , emit: for_multiqc
  path "fastqc/*_summary.csv"      , emit: collect
  path "logs/${task.process}/*.log", emit: log_files
  path "versions.yml"              , emit: versions
  val meta                         , emit: meta

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  def reads  = fastq.join(" ")
  """
    mkdir -p fastqc logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    fastqc ${args} \
      --outdir fastqc \
      --threads ${task.cpus} \
      ${reads} \
      --extract \
      | tee -a \$log_file

    header=\$(head -n 10 fastqc/*_fastqc/fastqc_data.txt | cut -f 1 | tr "\\n" ",")

    for data in fastqc/*_fastqc/fastqc_data.txt
    do
      if [ ! -f "fastqc/${prefix}_summary.csv" ]
      then
        head -n 10 \$data | \
          cut -f 1 | \
          tr "\\n" "," | \
          sed 's/,\$/\\n/' | \
          sed 's/#//g' | \
          sed 's/>//g' | \
          awk '{print "sample," \$0}' \
          > fastqc/${prefix}_summary.csv
      fi

      head -n 10 \$data | \
        cut -f 2 | \
        tr "\\n" "," | \
        sed 's/,\$/\\n/' | \
        awk -v sample=${prefix} '{print sample "," \$0}' \
        >> fastqc/${prefix}_summary.csv
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
    END_VERSIONS
  """
}
