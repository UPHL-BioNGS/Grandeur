process PANAROO {
  tag           "Core Genome Alignment"
  label         "process_high"
  container     'staphb/panaroo:1.5.0'
  
  input:
  file(gff)

  output:
  path "panaroo/*"                                                                         , emit: files
  tuple path("panaroo/core_gene_alignment.aln"), path("panaroo/gene_presence_absence.Rtab"), emit: core_gene_alignment
  path "logs/${task.process}/*.log"                                                        , emit: log_files
  path "versions.yml"                                                                      , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args       = task.ext.args   ?: '--clean-mode strict --remove-invalid-genes --alignment core'
  def prefix     = task.ext.prefix ?: 'panaroo'
  def assemblies = gff.join(' ')
  """
    mkdir -p input logs/${task.process}
    log_file=logs/${task.process}/${task.process}.${workflow.sessionId}.log

    for assembly in ${assemblies}
    do
      cp \$assembly input/.
    done

    ls input/* > input_genomes.txt

    panaroo ${args} \
      -t ${task.cpus} \
      -o ${prefix} \
      -i input_genomes.txt \
      | tee -a \$log_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        panaroo: \$(echo \$(panaroo --version 2>&1) | sed 's/^.*panaroo //' )
    END_VERSIONS
  """
}
