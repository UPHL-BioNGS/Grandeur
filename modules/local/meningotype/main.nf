process MENINGOTYPE {
  tag           "${meta.id}"
  label         "process_medium"
  container     'staphb/meningotype:0.8.5'

  input:
  tuple val(meta), file(contigs)

  output:
  path "meningotype/*.tsv", emit: files
  path "versions.yml", emit: versions
  val meta, emit: meta

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: '--finetype'
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p meningotype

    meningotype \
      ${args} \
      ${contigs} \
      > meningotype/${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        meningotype: \$(echo \$(meningotype --version 2>&1 | grep meningotype | awk '{print \$NF}'))
    END_VERSIONS
  """
}