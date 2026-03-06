process SEQSERO2 {
  tag           "${meta.id}"
  label         "process_medium"
  container     'staphb/seqsero2:1.3.2'


  input:
  tuple val(meta), file(file)

  output:
  tuple val(meta), file("seqsero2/*/*"), emit: files, optional: true
  path "seqsero2/*_seqsero_result.tsv", emit: collect, optional: true
  path "logs/${task.process}/*.log", emit: log
  path  "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args     ?: '-m a -b mem'
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p seqsero2 logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    SeqSero2_package.py ${args} \
      -m k \
      -t 4 \
      -i ${file} \
      -p ${task.cpus} \
      -d seqsero2/${prefix} \
      -n ${prefix} \
      | tee -a \$log_file

    if [ -f "seqsero2/${prefix}/SeqSero_result.tsv" ]
    then
      head -n 1 seqsero2/${prefix}/SeqSero_result.tsv | sed 's/Sample name/sample/g' > seqsero2/${prefix}_seqsero_result.tsv

      enteritidis_check=\$(grep "Enteritidis" seqsero2/${prefix}/SeqSero_result.tsv | head -n 1)
      sdf_check=\$(grep "Detected Sdf" seqsero2/${prefix}/SeqSero_result.tsv | head -n 1 )

      if [ -n "\$enteritidis_check" ] && [ -n "\$sdf_check" ]
      then
        tail -n 1 seqsero2/${prefix}/SeqSero_result.tsv | awk -F "\\t" -v OFS='\t' '{(\$9 = \$9 " (Sdf+)") ; print \$0}' >> seqsero2/${prefix}_seqsero_result.tsv
      elif [ -n "\$enteritidis_check" ] && [ -z "\$sdf_check" ]
      then
        tail -n 1 seqsero2/${prefix}/SeqSero_result.tsv | awk -F "\\t" -v OFS='\t' '{(\$9 = \$9 " (Sdf-)") ; print \$0}' >> seqsero2/${prefix}_seqsero_result.tsv
      else
        tail -n 1 seqsero2/${prefix}/SeqSero_result.tsv >> seqsero2/${prefix}_seqsero_result.tsv
      fi

    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqsero2: \$( echo \$( SeqSero2_package.py --version 2>&1) | sed 's/^.*SeqSero2_package.py //' )
    END_VERSIONS
  """
}
