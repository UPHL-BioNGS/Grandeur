process roary {
  tag           "Core Genome Alignment"
  label         "process_high"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/roary:3.13.0'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '10h'
  
  input:
  file(contigs)

  output:
  path "roary/*"                                                                       , emit: roary_files
  path "roary/fixed_input_files/*"                                                     , emit: roary_input_files
  tuple path("roary/core_gene_alignment.aln"), path("roary/gene_presence_absence.Rtab"), emit: core_gene_alignment
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log"                , emit: log_files
  path "versions.yml"                          , emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p logs/${task.process}
    log_file=logs/${task.process}/${task.process}.${workflow.sessionId}.log

    roary ${args} \
      -p ${task.cpus} \
      -f roary \
      -e -n \
      *.gff \
      | tee -a \$log_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        roary: \$( roary --version )
    END_VERSIONS
  """
}
