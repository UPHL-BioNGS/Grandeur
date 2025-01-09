process PROKKA {
  tag           "${meta.id}"
  label         "process_high"
  container     'staphb/prokka:1.14.6'

  input:
  tuple val(meta), file(contigs), val(organism)

  output:
  path "prokka/*/*"                 , emit: prokka_files
  path "prokka/*/*.txt"             , emit: for_multiqc
  path "gff/*.gff"                  , emit: gff, optional: true
  path "logs/${task.process}/*.log" , emit: log
  path "versions.yml"               , emit: versions
  val meta                          , emit: meta

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: '--mincontiglen 500 --compliant --locustag locus_tag --centre STAPHB'
  def prefix = task.ext.prefix ?: "${meta.id}"
  def gen_sp = organism ? "--genus ${organism[0]} --species ${organism[1]}" : ""
  """
    mkdir -p prokka gff logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    prokka ${args} \
      --cpu ${task.cpus} \
      --outdir prokka/${prefix} \
      --prefix ${prefix} \
      ${gen_sp} \
      --force ${contigs} \
      | tee -a \$log_file

    cp prokka/${prefix}/${prefix}.gff gff/${prefix}.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prokka: \$(echo \$(prokka --version 2>&1) | sed 's/^.*prokka //')
    END_VERSIONS
  """
}
