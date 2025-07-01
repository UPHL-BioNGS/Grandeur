process CIRCULOCOV {
  tag           "${meta.id}"
  label         "process_medium"
  //stageInMode   "copy"
  container     'staphb/circulocov:0.1.20240104'
  
  input:
  tuple val(meta), file(contigs), file(fastqs)

  output:
  tuple val(meta), file("circulocov/*/*sr.bam*"), emit: bam
  tuple val(meta), file(contigs), file("circulocov/*/*sr.bam*"), emit: contig_bam
  path "circulocov/*/overall_summary.txt", emit: collect
  path "circulocov/*", emit: everything
  path "circulocov/fastq/*", emit: fastq, optional: true
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  def reads  = fastqs.join(" ")
  """
    mkdir -p circulocov/${prefix} logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    circulocov ${args} \
      --threads ${task.cpus} \
      --genome ${contigs} \
      --illumina ${reads} \
      --out circulocov/${prefix} \
      --sample ${prefix} \
      | tee -a \$log_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      circulocov: \$(circulocov -v)
    END_VERSIONS
  """
}
