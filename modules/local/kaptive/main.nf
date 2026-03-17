process KAPTIVE {
  tag           "${meta.id}"
  label         "process_medium"
  container     'staphb/kaptive:3.1.0'

  input:
  tuple val(meta), file(contigs)

  output:
  tuple val(meta), file("kaptive/*.txt"), emit: files, optional: true
  path "kaptive/*_kaptive.tsv", emit: collect, optional: true
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  mkdir -p kaptive logs/${task.process}
  log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

  header="sample"
  result="${prefix}"

  for ref in \$(ls /kaptive/reference_database/*.gbk | cut -f 4 -d "/" | sed 's/.gbk//g')
  do
    echo "Running kaptive against \${ref}.gbk " | tee -a \$log_file
    kaptive \
      ${args} \
      assembly \
      /kaptive/reference_database/\${ref}.gbk \
      ${contigs} \
      --threads 6 \
      --out kaptive/${prefix}_\${ref}.txt  \
      | tee -a \$log_file

    new_header=\$(head -n 1 kaptive/${prefix}_\${ref}.txt | sed "s/\t/\t\${ref}_/g")
    line_count=\$(wc -l < kaptive/${prefix}_\${ref}.txt)
    if [ "\$line_count" -gt 1 ]
    then
      new_result=\$(tail -n +2 kaptive/${prefix}_\${ref}.txt)
    else
      new_result="\\t\\t\\t\\t\\t\\t\\t\\t\\t\\t\\t\\t\\t\\t\\t\\t\\t\\t"
    fi

    header="\$header\\tGBK\\t\$new_header"
    result="\$result\\t\${ref}.gbk\\t\$new_result"
  done

  # combining all files into one
  echo -e "\$header" >  kaptive/${prefix}_kaptive.tsv
  echo -e "\$result" >> kaptive/${prefix}_kaptive.tsv

  cat <<-END_VERSIONS > versions.yml
  "GRANDEUR:SUBTYPING:KAPTIVE":
    kaptive: \$( echo \$(kaptive --version | sed 's/Kaptive v//;'))
  END_VERSIONS
  """
}