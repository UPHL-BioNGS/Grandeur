process BLASTN {
  tag           "${meta.id}"
  label         "process_medium"
  container     'staphb/blast:2.17.0'

  input:
  tuple val(meta), file(contig), path(blastdb)

  output:
  tuple val(meta), file("blastn/*.tsv"), emit: blastn
  path "logs/${task.process}/*.${workflow.sessionId}.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: '-max_target_seqs 10 -max_hsps 1 -evalue 1e-25'
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p blastn logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    blastn -query ${contig} \
      -out blastn/${prefix}.tsv \
      -num_threads ${task.cpus} \
      -db ${blastdb}/${params.blast_db_type} \
      -outfmt '6 qseqid staxids bitscore std' \
      ${args} \
      | tee -a \$log_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
  """
}
