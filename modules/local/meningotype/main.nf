process MENINGOTYPE {
  tag           "${meta.id}"
  label         "process_medium"
  container     'staphb/meningotype:0.8.5'

  input:
  tuple val(meta), file(contigs)

  output:
  tuple val(meta), file("meningotype/*.tsv"), emit: files
  path "*meningotype.tsv", emit: summary
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

    # ensure prefix is in summary file
    head -n 1  meningotype/${prefix}.tsv | awk '{print "sample\\t"    \$0 }' >  ${prefix}_meningotype.tsv
    tail -n +2 meningotype/${prefix}.tsv | awk '{print "${prefix}\\t" \$0 }' >> ${prefix}_meningotype.tsv


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        meningotype: \$(echo \$(meningotype --version 2>&1 | grep meningotype | awk '{print \$NF}'))
    END_VERSIONS
  """
}