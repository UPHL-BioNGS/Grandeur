process DATASETS_SUMMARY {
  tag           "${taxon}"
  label         "process_single"
  container     'staphb/ncbi-datasets:18.18.0'

  input:
  tuple val(taxon), file(script)

  output:
  path "datasets/*_genomes.csv", emit: genomes, optional: true
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args     = task.ext.args  ?: '--reference --mag exclude --limit 5'
  def fields   = task.ext.fields ?: 'accession,assminfo-refseq-category,assminfo-level,organism-name,assmstats-total-ungapped-len'
  def gen_spec = taxon.toString().replace('_', ' ')
  def prefix   = task.ext.prefix ?: taxon
  """
    mkdir -p datasets logs/!{task.process}

    echo ${fields} > datasets/${prefix}_genomes.csv

    datasets summary genome taxon "${gen_spec}" ${args} --as-json-lines | \
      dataformat tsv genome --fields ${fields} | \
      awk '{if (\$NF < 15000000 ) print \$0}' | \
      sort | uniq | tr "\\t" "," \
      >> datasets/${prefix}_genomes.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        datasets: \$(datasets --version | awk '{print \$NF}' )
        dataformat: \$(dataformat version | awk '{print \$NF}' )
    END_VERSIONS
  """
}
