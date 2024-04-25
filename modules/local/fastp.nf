process fastp {
  tag           "${meta.id}"
  label         "process_low"
  publishDir    params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/fastp:0.23.4'
  time          '30m'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  
  input:
  tuple val(meta), file(reads)

  output:
  tuple val(meta), file("fastp/*_fastp_R{1,2}.fastq.gz"), emit: fastq, optional: true
  path "fastp/*_fastp.html",                              emit: html, optional: true
  path "fastp/*_fastp.json",                              emit: fastp_files, optional: true
  path "logs/${task.process}/*.{log,err}",                emit: log
  tuple val(meta), env(passed_reads),                     emit: fastp_results
  path  "versions.yml",                                   emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
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
    if [ -z "\$passed_reads" ] ; then passed_reads="0" ; fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      fastp: \$(fastp --version 2>&1 | sed -e 's/fastp //g')
    END_VERSIONS
  """
}
