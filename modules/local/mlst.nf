process MLST {
  tag           "${meta.id}"
  label         "process_medium"
  container     'staphb/mlst:2.23.0-2025-01-01'

  input:
  tuple val(meta), file(contig), file(script)

  output:
  path "mlst/*_mlst.tsv", emit: collect
  path "versions.yml"   , emit: versions
  val meta              , emit: meta

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p mlst 

    mlst ${args} \
      --threads ${task.cpus} \
      ${contig} | \
      tr ' ' '_' \
      > ${prefix}_mlst.txt

    python3 ${script} ${prefix}_mlst.txt mlst/${prefix}_mlst.tsv mlst ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mlst: \$( echo \$(mlst --version 2>&1) | sed 's/mlst //' )
    END_VERSIONS
  """
}
