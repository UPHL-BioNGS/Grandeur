process REFERENCES {
  tag       "Preparing references"
  label     "process_single"
  container 'staphb/grandeur_ref:4.6'

  output:
  path "ref/*", emit: fastas

  when:
  task.ext.when == null || task.ext.when

  script:
  """
  mkdir ref

  cp /ref/* ref/.
  """
}
