process blastn {
  tag           "$meta.id"
  label         "process_medium"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/blast:2.15.0'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize' , value: 'himem-medium'
  //#UPHLICA cpus 8
  //#UPHLICA memory 100.GB
  //#UPHLICA time '2h'

  input:
  tuple val(sample), file(contig), path(blastdb)

  output:
  tuple val(sample), file("blastn/${sample}.tsv")                , emit: blastn
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log", emit: log
  path  "versions.yml"                          , emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
      def args = task.ext.args ?: '-max_target_seqs 10 -max_hsps 1 -evalue 1e-25'
    def prefix = task.ext.prefix ?: "${meta.id}"
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
      -db !{blastdb}/!{params.blast_db_type} \
      -outfmt '6 qseqid staxids bitscore std' \
      !{params.blastn_options} \
      | tee -a $log_file
  '''
}
