process SNPDISTS {
  tag           "${aln}"
  label         "process_medium"
  container     'staphb/snp-dists:1.2.0'

  input:
  file(aln)

  output:
  path "snp-dists/snpdists*.txt", emit: snp_matrix, optional: true
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: '-c'
  def prefix = task.ext.prefix ?: "snpdists_${aln.baseName}"
  """
    mkdir -p snp-dists

    snp-dists ${args} \
      ${aln} \
      > snp-dists/${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpdists: \$(snp-dists -v 2>&1 | sed 's/snp-dists //;')
    END_VERSIONS
  """
}
