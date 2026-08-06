process KLEBORATE {
  tag           "${meta.id}"
  label         "process_low"
  container     'staphb/kleborate:3.2.4-micromamba'

  input:
  tuple val(meta), file(contig), file(script)

  output:
  path "kleborate/*_kleborate.tsv", emit: collect
  tuple val(meta), file("kleborate/*_output.txt"), emit: result, optional: true
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: '--trim_headers'
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  mkdir -p kleborate logs/${task.process}
  log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

  cat ${contig} > input_${prefix}.fasta

  header="sample"
  result="${prefix}"

  for module in \$(kleborate -o kleborate --list_modules | grep ":" | cut -f 1 -d ":" | sed 's/\\x1b\\[[0-9;]*m//g' )
  do
    echo "Running kleborate with module \$module" | tee -a \$log_file

    kleborate ${args} \
      -m \$module \
      -o kleborate \
      -a input_${prefix}.fasta \
      | tee -a \$log_file

    new_header=\$(head -n 1 kleborate/\${module}_output.txt)
    new_result=\$(tail -n 1 kleborate/\${module}_output.txt)
    header="\$header\\t\$new_header"
    result="\$result\\t\$new_result"
  done

  # combining all files into one
  echo -e "\$header" >  kleborate/${prefix}_kleborate.tsv
  echo -e "\$result" >> kleborate/${prefix}_kleborate.tsv

  cat <<-END_VERSIONS > versions.yml
  "GRANDEUR:SUBTYPING:KLEBORATE":
    kleborate: \$( echo \$(kleborate --version | sed 's/Kleborate v//;'))
  END_VERSIONS
  """
}
