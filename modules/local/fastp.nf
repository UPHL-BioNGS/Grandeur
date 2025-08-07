process FASTP {
  tag           "${meta.id}"
  label         "process_low"
  container     'staphb/fastp:1.0.1'
  
  input:
  tuple val(meta), file(reads)

  output:
  tuple val(meta), file("fastp/*_fastp_R{1,2}.fastq.gz"), emit: fastq, optional: true
  path "fastp/*_fastp.html", emit: html, optional: true
  path "fastp/*_fastp.json", emit: fastp_files, optional: true
  path "logs/${task.process}/*.{log,err}", emit: log
  path  "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: '--detect_adapter_for_pe'
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p fastp logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log
    err_file=logs/${task.process}/${prefix}.${workflow.sessionId}.err

    fastp ${args} \
      -i ${reads[0]} \
      -I ${reads[1]} \
      -o fastp/${prefix}_fastp_R1.fastq.gz \
      -O fastp/${prefix}_fastp_R2.fastq.gz \
      -h fastp/${prefix}_fastp.html \
      -j fastp/${prefix}_fastp.json \
      2>> \$err_file | tee -a \$log_file

    passed_reads=\$(grep "reads passed filter" \$err_file | tail -n 1 | cut -f 2 -d ":" | sed 's/ //g' )

    # removes reads from proceeding further if there aren't enough
    if [ "\$passed_reads" -lt "${params.minimum_reads}" ]
    then
      mv fastp/${prefix}_fastp_R1.fastq.gz fastp/${prefix}_fastp_1.fastq.gz
      mv fastp/${prefix}_fastp_R2.fastq.gz fastp/${prefix}_fastp_2.fastq.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      fastp: \$(fastp --version 2>&1 | sed -e 's/fastp //g')
    END_VERSIONS
  """
}
