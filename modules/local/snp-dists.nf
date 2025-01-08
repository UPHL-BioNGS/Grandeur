process SNPDISTS {
  tag           "SNP matrix"
  label         "process_medium"
  container     'staphb/snp-dists:0.8.2'

  
  input:
  file(contigs)

  output:
  path "snp-dists/snp_matrix.txt", emit: snp_matrix
  path "versions.yml"            , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: '-c'
  """
    mkdir -p snp-dists

    snp-dists ${args} \
      ${contigs} \
      > snp-dists/snp_matrix.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpdists: \$(snp-dists -v 2>&1 | sed 's/snp-dists //;')
    END_VERSIONS
  """
}
