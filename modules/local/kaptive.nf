process KAPTIVE {
  tag           "${meta.id}"
  label         "process_medium"
  container     'staphb/kaptive:3.1.0'

  input:
  tuple val(meta), file(contigs)

  output:
  path "kaptive/${meta.id}_table.txt", emit: collect
  path "kaptive/*", emit: files
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p kaptive logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    kaptive\
      assembly \
      ${args} \
      /kaptive/reference_database/VibrioPara_Kaptivedb_K.gbk \
      ${contigs} \
      --threads ${task.cpus} \
      --out kaptive/${prefix}_VibrioPara_Kaptivedb_K.txt \
      | tee -a \$log_file

    kaptive \
      assembly \
      ${args} \
      /kaptive/reference_database/VibrioPara_Kaptivedb_O.gbk \
      ${contigs} \
      --threads ${task.cpus} \
      --out kaptive/${prefix}_VibrioPara_Kaptivedb_O.txt \
      | tee -a \$log_file

    grep -h "Other genes" kaptive/${prefix}* | head -n 1 > ${prefix}_table.txt
    grep -h ${prefix} kaptive/${prefix}* >> ${prefix}_table.txt
    mv ${prefix}_table.txt kaptive/${prefix}_table.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      kaptive: \$( echo \$(kaptive --version | sed 's/Kaptive v//;'))
    END_VERSIONS
  """
}