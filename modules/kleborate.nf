process kleborate {
  tag "${sample}"
  label "medcpus"

  when:
  params.contig_processes =~ /kleborate/ && flag =~ 'found'

  input:
  tuple val(sample), file(contig), val(flag)

  output:
  tuple val(sample), env(kleborate_score)                               , emit: score
  tuple val(sample), env(kleborate_mlst)                                , emit: mlst
  path "kleborate/${sample}_results.txt"                          , emit: collect
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log" , emit: log

  shell:
  '''
    mkdir -p kleborate logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    kleborate --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    kleborate !{params.kleborate_options} \
      -o kleborate/!{sample}_results.txt \
      -a !{contig} \
      | tee -a $log_file

    virulence_column=$(head -n 1 kleborate/!{sample}_results.txt | tr '\\t' '\\n' | grep -n virulence_score | cut -f 1 -d ":" )
    mlst_column=$(head -n 1 kleborate/!{sample}_results.txt | tr '\\t' '\\n' | grep -n ST | cut -f 1 -d ":" )
    if [ -n "$virulence_column" ] ; then kleborate_score=$(cut -f $virulence_column kleborate/!{sample}_results.txt | tail -n 1 ) ; fi
    if [ -n "$mlst_column" ] ; then kleborate_mlst=$(cut -f $mlst_column kleborate/!{sample}_results.txt | tail -n 1 ) ; fi
    if [ -z "$kleborate_score" ] ; then kleborate_score='NA' ; fi
    if [ -z "$kleborate_mlst" ] ; then kleborate_mlst='NA' ; fi
  '''
}
