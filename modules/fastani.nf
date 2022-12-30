process fastani {
  tag "${sample}"
  label "medcpus"

  input:
  tuple val(sample), file(contigs), file(genomes)

  output:
  path "fastani/${sample}.txt"                   , emit: results
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log", emit: log

  shell:
  '''
    mkdir -p fastani logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    echo "fastANI version: " >> $log_file
    fastANI --version 2>> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    mkdir db

    cp !{genomes}/* 2>/dev/null db/.

    ls db/{*fna,*fa,*fasta} 2>/dev/null > reference_list.txt

    fastANI \
      -q !{contigs} \
      --rl reference_list.txt \
      -o test

    echo "WTF?"

      exit 1

    fastANI !{params.fastani_options} \
      --threads !{task.cpus} \
      -q !{contigs} \
      --rl reference_list.txt \
      -o fastani/!{sample}.txt \
      | tee -a $log_file
  '''
}

//total=$(sort     -k3,3n -k 4,4n fastani/!{sample}.txt | tail -n 1 | cut -f 5 )
