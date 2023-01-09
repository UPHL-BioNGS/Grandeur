process seqsero2 {
  tag           "${sample}"
  label         "medcpus"
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/seqsero2:1.2.1'
  maxForks      10
  
  when:
  flag =~ 'found'

  input:
  tuple val(sample), file(file), val(flag)

  output:
  path "seqsero2/${sample}/*"                                    , emit: files
  path "seqsero2/${sample}/SeqSero_result.tsv"                   , emit: collect
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log", emit: log

  shell:
  '''
    mkdir -p seqsero2 logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    SeqSero2_package.py --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    SeqSero2_package.py !{params.seqsero2_options} \
      -m k \
      -t 4 \
      -i !{file} \
      -p !{task.cpus} \
      -d seqsero2/!{sample} \
      -n !{sample} \
      | tee -a $log_file

    enteritidis_check=$(grep "Enteritidis" seqsero2/!{sample}/SeqSero_result.tsv | head -n 1)
    sdf_check=$(grep "Detected Sdf" seqsero2/!{sample}/SeqSero_result.tsv | head -n 1 )

    if [ -n "$enteritidis_check" ] && [ -n "$sdf_check" ]
    then
      head -n 1 seqsero2/!{sample}/SeqSero_result.tsv > SeqSero_result.tsv.tmp
      tail -n 1 seqsero2/!{sample}/SeqSero_result.tsv | awk -F "\\t" -v OFS='\t' '{($9 = $9 " (Sdf+)") ; print $0}' >> SeqSero_result.tsv.tmp
      mv SeqSero_result.tsv.tmp seqsero2/!{sample}/SeqSero_result.tsv
    elif [ -n "$enteritidis_check" ] && [ -z "$sdf_check" ]
    then
      head -n 1 seqsero2/!{sample}/SeqSero_result.tsv > SeqSero_result.tsv.tmp
      tail -n 1 seqsero2/!{sample}/SeqSero_result.tsv | awk -F "\\t" -v OFS='\t' '{($9 = $9 " (Sdf-)") ; print $0}' >> SeqSero_result.tsv.tmp
      mv SeqSero_result.tsv.tmp seqsero2/!{sample}/SeqSero_result.tsv
    fi

    cat seqsero2/!{sample}/SeqSero_result.tsv | sed 's/Sample name/sample/g' > seqsero2/!{sample}/SeqSero_result.tsv.tmp
    mv seqsero2/!{sample}/SeqSero_result.tsv.tmp seqsero2/!{sample}/SeqSero_result.tsv
  '''
}
