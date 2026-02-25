process NAMES {
  tag           "${meta.id}"
  label         "process_single"
  container     'staphb/pandas:2.3.0'
  
  input:
  tuple val(meta), file(input)

  output:
  path "summary/*_names.csv", emit: collect

  when:
  task.ext.when == null || task.ext.when

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def files  = input.join(" ")
  """
  mkdir -p summary

  echo "sample,file,version" > summary/${prefix}_names.csv
  echo "${prefix},${files},${workflow.manifest.version}" >> summary/${prefix}_names.csv
  """
}
