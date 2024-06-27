process kaptive {
  tag           "${meta.id}"
  label         "process_medium"
  stageInMode   "copy"
  publishDir    path: params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/kaptive:2.0.8'
  time          '30m'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}

  input:
  tuple val(meta), file(contigs)

  output:
  path "kaptive/${meta.id}_table.txt", emit: collect
  path "kaptive/*", emit: files
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  (task.ext.when == null || task.ext.when)

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p kaptive logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    kaptive.py ${args} \
      --k_refs /kaptive/reference_database/VibrioPara_Kaptivedb_K.gbk \
      --assembly ${contigs} \
      --threads ${task.cpus} \
      --out kaptive/${prefix}_VibrioPara_Kaptivedb_K \
      | tee -a \$log_file

    kaptive.py ${args} \
      --k_refs /kaptive/reference_database/VibrioPara_Kaptivedb_O.gbk \
      --assembly ${contigs} \
      --threads ${task.cpus} \
      --out kaptive/${prefix}_VibrioPara_Kaptivedb_O \
      | tee -a \$log_file

    grep -h "Other genes" kaptive/*table.txt | head -n 1 > ${prefix}_table.txt
    grep -h ${prefix} kaptive/*table.txt >> ${prefix}_table.txt
    mv ${prefix}_table.txt kaptive/${prefix}_table.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      kaptive.py: \$( echo \$(kaptive.py --version | sed 's/Kaptive v//;'))
    END_VERSIONS
  """
}