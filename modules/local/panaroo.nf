process panaroo {
  tag           "Core Genome Alignment"
  label         "process_high"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/panaroo:1.3.4'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '10h'
  
  input:
  file(gff)

  output:
  path "panaroo/*"                                                                         , emit: files
  tuple path("panaroo/core_gene_alignment.aln"), path("panaroo/gene_presence_absence.Rtab"), emit: core_gene_alignment
  path "logs/${task.process}/*.log"                                                        , emit: log_files
  path "versions.yml"                                                                      , emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args       = task.ext.args   ?: '--clean-mode strict --remove-invalid-genes'
  def prefix     = task.ext.prefix ?: "panaroo"
  def assemblies = gff.join(' ')
  """
    mkdir -p logs/${task.process}
    log_file=logs/${task.process}/${task.process}.${workflow.sessionId}.log

    panaroo ${args} \
      -t ${task.cpus} \
      -o ${prefix} \
      -i ${assemblies} \
      -a core \
      | tee -a \$log_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        panaroo: \$(echo \$(panaroo --version 2>&1) | sed 's/^.*panaroo //' )
    END_VERSIONS
  """
}
