process mlst {
  tag           "${meta.id}"
  label         "process_medium"
  publishDir    params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/mlst:2.23.0-2024-03-11'
  maxForks      10
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '10m'

  input:
  tuple val(meta), file(contig), file(script)

  output:
  path "mlst/*_mlst.tsv", emit: collect
  path "versions.yml"   , emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p mlst 

    mlst ${args} \
      --threads ${task.cpus} \
      ${contig} \
      > ${prefix}_mlst.txt

    python3 ${script} ${prefix}_mlst.txt mlst/${prefix}_mlst.tsv mlst ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mlst: \$( echo \$(mlst --version 2>&1) | sed 's/mlst //' )
    END_VERSIONS
  """
}
