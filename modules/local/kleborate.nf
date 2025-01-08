process KLEBORATE {
  tag           "${meta.id}"
  label         "process_medium"
  container     'staphb/kleborate:3.1.2'

  input:
  tuple val(meta), file(contig), file(script)

  output:
  path "kleborate/*_results.tsv"   , emit: collect, optional: true
  path "kleborate/*_results.txt"   , emit: result, optional: true
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml"              , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: '-p kpsc --trim_headers'
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p kleborate logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    kleborate ${args} \
      -o kleborate/${prefix} \
      -a ${contig} \
      | tee -a \$log_file

    if [ -f "kleborate/${prefix}_results.txt" ]
    then
      python3 ${script} kleborate/${prefix}_results.txt kleborate/${prefix}_results.tsv kleborate ${prefix}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kleborate: \$( echo \$(kleborate --version | sed 's/Kleborate v//;'))
    END_VERSIONS
  """
}
