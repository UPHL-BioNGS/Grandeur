process SEQSERO2S {
  tag           "${meta.id}"
  label         "process_medium"
  container     'staphb/seqsero2s:1.1.4'


  input:
  tuple val(meta), file(file)

  output:
  tuple val(meta), file("seqsero2s/*/*"), emit: files, optional: true
  path "seqsero2s/*_seqsero2s_result.tsv", emit: collect, optional: true
  path "logs/${task.process}/*.log", emit: log
  path  "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args     ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p seqsero2s logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    SeqSero2S.py \
      -i ${file} \
      -t 4 \
      -m k \
      -d seqsero2s/${prefix} \
      -n ${prefix} \
      -p ${task.cpus}  \
      | tee -a \$log_file


    if [ -f "seqsero2s/${prefix}/SeqSero_result.tsv" ]
    then
      head -n 1 seqsero2s/${prefix}/SeqSero_result.tsv | sed 's/Sample name/sample/g' > seqsero2s/${prefix}_seqsero2s_result.tsv

      enteritidis_check=\$(grep "Enteritidis" seqsero2s/${prefix}/SeqSero_result.tsv | head -n 1)
      sdf_check=\$(grep "Detected Sdf" seqsero2s/${prefix}/SeqSero_result.tsv | head -n 1 )

      if [ -n "\$enteritidis_check" ] && [ -n "\$sdf_check" ]
      then
        tail -n 1 seqsero2s/${prefix}/SeqSero_result.tsv | awk -F "\\t" -v OFS='\t' '{(\$9 = \$9 " (Sdf+)") ; print \$0}' >> seqsero2s/${prefix}_seqsero2s_result.tsv
      elif [ -n "\$enteritidis_check" ] && [ -z "\$sdf_check" ]
      then
        tail -n 1 seqsero2s/${prefix}/SeqSero_result.tsv | awk -F "\\t" -v OFS='\t' '{(\$9 = \$9 " (Sdf-)") ; print \$0}' >> seqsero2s/${prefix}_seqsero2s_result.tsv
      else
        tail -n 1 seqsero2s/${prefix}/SeqSero_result.tsv >> seqsero2s/${prefix}_seqsero2s_result.tsv
      fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqsero2S: \$( SeqSero2S.py --version | awk '{print \$NF}' )
    END_VERSIONS
  """
}
