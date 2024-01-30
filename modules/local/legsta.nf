process legsta {
  tag           "${meta.id}"
  label         "process_medium"
  stageInMode   "copy"
  publishDir    path: params.outdir, mode: 'copy'
  container     'staphb/legsta:0.5.1'
  maxForks      10
  time          '30m'
  //errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}

  when:
  (task.ext.when == null || task.ext.when) && flag =~ 'found'

  input:
  tuple val(meta), file(contigs), val(flag)

  output:
  path "legsta/*_legsta.csv"       , emit: collect
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml"              , emit: versions

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p legsta logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log
    
    cat ${contigs} | awk '{(\$1=\$1); print \$0}' > ${prefix}.fasta.tmp
    mv ${prefix}.fasta.tmp ${prefix}.fasta

    legsta ${args} \
      ${prefix}.fasta \
      --csv \
      | tee -a \$log_file \
      | awk -v sample=${prefix} '{print sample "," \$0 }' \
      | sed '0,/${prefix}/{s/${prefix}/sample/}' > legsta/${prefix}_legsta.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        legsta: \$(echo \$(legsta --version 2>&1) | sed 's/^.*legsta //; s/ .*\$//;')
    END_VERSIONS
  """
}