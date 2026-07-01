process JSON_CONVERT {
  tag       "${meta.id}"
  label     "process_single"
  container 'staphb/pandas:3.0.3'

  input:
  tuple val(meta), val(analysis), file(json), file(script)

  output:
  path "${analysis}/*_${analysis}*", emit: collect, optional: true

  when:
  task.ext.when == null || task.ext.when

  script:
  """
  mkdir -p ${analysis}

  python3 ${script} ${json} ${analysis}

  mv *${analysis}*tsv ${analysis}/.
  """
}
