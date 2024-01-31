process snp_dists {
  tag           "SNP matrix"
  label         "process_medium"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/snp-dists:0.8.2'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '24h'
  
  input:
  file(contigs)

  output:
  path "snp-dists/snp_matrix.txt", emit: snp_matrix
  path "versions.yml"            , emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
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
