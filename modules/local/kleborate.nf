process kleborate {
  tag           "${meta.id}"
  label         "process_medium"
  publishDir    params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/kleborate:2.4.1'
  time          '10m'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}

  input:
  tuple val(meta), file(contig), file(script)

  output:
  path "kleborate/*_results.tsv"   , emit: collect, optional: true
  path "kleborate/*_results.txt"   , emit: result
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml"              , emit: versions

  when:
  (task.ext.when == null || task.ext.when)

  shell:
  def args   = task.ext.args   ?: '--all'
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p kleborate logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    kleborate ${args} \
      -o kleborate/${prefix}_results.txt \
      -a ${contig} \
      | tee -a \$log_file

    python3 ${script} kleborate/${prefix}_results.txt kleborate/${prefix}_results.tsv kleborate ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kleborate: \$( echo \$(kleborate --version | sed 's/Kleborate v//;'))
    END_VERSIONS
  """
}
