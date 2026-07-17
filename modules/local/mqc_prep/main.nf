process MQC_PREP {
  tag           "prepping files"
  label         "process_single"
  container     'staphb/pandas:3.0.3'
  
  input:
  file(input)
  file(script)

  output:
  path "*mqc*", emit: for_multiqc, optional: true

  when:
  task.ext.when == null || task.ext.when

  script:
  """
  python3 ${script}
  """
}
