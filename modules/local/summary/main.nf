process SUMMARY {
  tag           "Creating summary files"
  container     'staphb/pandas:3.0.1'
  label         "process_single"

  input:
  file(input)

  output:
  path "grandeur_summary.tsv"                 , emit: summary_tsv, optional: true
  path "grandeur_summary.txt"                 , emit: summary_txt, optional: true
  path "summary/grandeur_extended_summary.tsv", emit: extended_tsv, optional: true
  path "summary/grandeur_extended_summary.txt", emit: extended_txt, optional: true

  when:
  task.ext.when == null || task.ext.when

  script:
  """
    mkdir -p summary

    python3 summary.py
  """
}
