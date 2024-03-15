process seqsero2 {
  tag           "${meta.id}"
  label         "process_medium"
  publishDir    params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/seqsero2:1.2.1'
  time          '10m'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}

  input:
  tuple val(meta), file(file), val(flag)

  output:
  path "seqsero2/*/*"                                    , emit: files
  path "seqsero2/*/SeqSero_result.tsv"                   , emit: collect
  path "logs/${task.process}/*.log", emit: log
  path  "versions.yml"                          , emit: versions

  when:
  (task.ext.when == null || task.ext.when) && flag =~ 'found'

  shell:
  def args = task.ext.args     ?: '-m a -b mem'
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p seqsero2 logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > \$log_file
    echo "container : ${task.container}" >> \$log_file
    SeqSero2_package.py --version >> \$log_file
    echo "Nextflow command : " >> \$log_file
    cat .command.sh >> \$log_file

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
      enteritidis_check=\$(grep "Enteritidis" seqsero2/${prefix}/SeqSero_result.tsv | head -n 1)
      sdf_check=\$(grep "Detected Sdf" seqsero2/${prefix}/SeqSero_result.tsv | head -n 1 )

      if [ -n "\$enteritidis_check" ] && [ -n "\$sdf_check" ]
      then
        head -n 1 seqsero2/${prefix}/SeqSero_result.tsv > SeqSero_result.tsv.tmp
        tail -n 1 seqsero2/${prefix}/SeqSero_result.tsv | awk -F "\\t" -v OFS='\t' '{(\$9 = \$9 " (Sdf+)") ; print \$0}' >> SeqSero_result.tsv.tmp
        mv SeqSero_result.tsv.tmp seqsero2/${prefix}/SeqSero_result.tsv
      elif [ -n "\$enteritidis_check" ] && [ -z "\$sdf_check" ]
      then
        head -n 1 seqsero2/${prefix}/SeqSero_result.tsv > SeqSero_result.tsv.tmp
        tail -n 1 seqsero2/${prefix}/SeqSero_result.tsv | awk -F "\\t" -v OFS='\t' '{(\$9 = \$9 " (Sdf-)") ; print \$0}' >> SeqSero_result.tsv.tmp
        mv SeqSero_result.tsv.tmp seqsero2/${prefix}/SeqSero_result.tsv
      fi

      cat seqsero2/${prefix}/SeqSero_result.tsv | sed 's/Sample name/sample/g' > seqsero2/${prefix}/SeqSero_result.tsv.tmp
      mv seqsero2/${prefix}/SeqSero_result.tsv.tmp seqsero2/${prefix}/SeqSero_result.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqsero2: \$( echo \$( SeqSero2_package.py --version 2>&1) | sed 's/^.*SeqSero2_package.py //' )
    END_VERSIONS
  """
}
