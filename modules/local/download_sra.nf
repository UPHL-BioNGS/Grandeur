process download_sra {
  tag           "${SRR}"
  label         "process_single"
  publishDir    params.outdir, mode: 'copy'
  container     'quay.io/biocontainers/sra-tools:2.11.0--pl5321ha49a11a_3'
  time          '2h'
  //errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  
  input:
  val(SRR)

  output:
  tuple val(SRR), file("reads/${SRR}_{1,2}.fastq.gz"), emit: fastq
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p reads logs/${task.process}
    log_file=logs/${task.process}/${SRR}.${workflow.sessionId}.log

    fasterq-dump ${args} \
      ${SRR} \
      --split-files \
      --threads ${task.cpus} \
      --outdir reads | \
      tee -a \$log_file

    gzip reads/${SRR}_1.fastq
    gzip reads/${SRR}_2.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      sratools: \$(fasterq-dump --version 2>&1 | grep -Eo '[0-9.]+')
    END_VERSIONS
  """
}

