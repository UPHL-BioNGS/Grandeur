process fastani {
  tag "${sample}"
  label "medcpus"

  when:
  params.contig_processes =~ /fastani/

  input:
  tuple val(sample), file(contigs), file(genomes)

  output:
  path "fastani/${sample}.txt", optional: true                          , emit: collect
  path "logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}" , emit: log
  tuple val(sample), env(top_ref)                                       , emit: ref
  tuple val(sample), env(ani_score)                                     , emit: ani
  tuple val(sample), env(fragment)                                      , emit: fragment
  tuple val(sample), env(total)                                         , emit: total

  shell:
  '''
    mkdir -p fastani logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    echo "fastANI version: " >> $log_file
    fastANI --version 2>> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    tar -xvf !{genomes}
    ls genomes/*fna > reference_list.txt

    fastANI !{params.fastani_options} \
      --threads !{task.cpus} \
      -q !{contigs} \
      --rl reference_list.txt \
      -o fastani/!{sample}.txt \
      2>> $err_file | tee -a $log_file

    top_ref=$(sort   -k3,3n -k 4,4n fastani/!{sample}.txt | tail -n 1 | cut -f 2 )
    ani_score=$(sort -k3,3n -k 4,4n fastani/!{sample}.txt | tail -n 1 | cut -f 3 )
    fragment=$(sort  -k3,3n -k 4,4n fastani/!{sample}.txt | tail -n 1 | cut -f 4 )
    total=$(sort     -k3,3n -k 4,4n fastani/!{sample}.txt | tail -n 1 | cut -f 5 )

    if [ -z "$top_ref"   ] ; then top_ref='not identified' ; fi
    if [ -z "$ani_score" ] ; then ani_score='0'            ; fi
    if [ -z "$fragment"  ] ; then fragment='0'             ; fi
    if [ -z "$total"     ] ; then total='0'                ; fi
  '''
}
