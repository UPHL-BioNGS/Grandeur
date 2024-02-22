process circulocov {
  tag           "${meta.id}"
  label         "process_medium"
  stageInMode   "copy"
  publishDir    path: params.outdir, mode: 'copy'
  container     'quay.io/uphl/circulocov:0.1.20240104-2024-02-21'
  time          '30m'
  //errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  
  input:
  tuple val(meta), file(fastqs), file(contigs)

  output:
  tuple val(meta), file("circulocov/*/*sr.bam*"), emit: bam
  path "circulocov/*/overall_summary.txt", emit: collect
  path "circulocov/*", emit: everything
  path "circulocov/fastq/*", emit: fastq, optional: true
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
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
