process blastn {
  tag "${sample}"
  label "medcpus"

  input:
  tuple val(sample), file(contig), path(blastdb)

  output:
  tuple val(sample), file("blastn/${sample}.tsv")                       , emit: blastn
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log" , emit: log

  shell:
  '''
    mkdir -p blastn logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    blastn -version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    blastn -query !{contig} \
      -out blastn/!{sample}.tsv \
      -num_threads !{task.cpus} \
      -db !{blastdb}/!{params.local_db_type} \
      -outfmt '6 qseqid staxids bitscore std' \
      -max_target_seqs 10 \
      -max_hsps 1 \
      -evalue 1e-25 \
      | tee -a $log_file
  '''
}
