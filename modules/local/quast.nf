process QUAST {
  tag           "${meta.id}"
  label         "process_medium"
  container     'staphb/quast:5.3.0'

  input:
  tuple val(meta), file(contigs), file(reads)

  output:
  path "quast/*"                                                   , emit: files
  path "quast/*_quast_report.tsv"                  , optional: true, emit: for_multiqc
  tuple val(meta), file("quast/*_quast_report.tsv"), optional: true, emit: results
  path "quast/*/quast_transposed_report.tsv"             , optional: true, emit: collect
  path "quast/*/quast_transposed_report_contig.tsv"      , optional: true, emit: collect_contig
  path "logs/${task.process}/*.log"                                , emit: log
  path "versions.yml"                                              , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  def fastq  = reads[1] ? "--pe1 ${reads[0]} --pe2 ${reads[1]}" : ""
  def fin    = reads[1] ? "quast/${prefix}/quast_transposed_report.tsv" : "quast/${prefix}/quast_transposed_report_contig.tsv"
  """
    mkdir -p ${task.process} logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    quast.py ${args} \
      ${contigs} \
      --output-dir quast/${prefix} \
      --threads ${task.cpus} \
      ${fastq} \
      | tee -a \$log_file

    if [ -f "quast/${prefix}/report.tsv" ] ; then cp quast/${prefix}/report.tsv quast/${prefix}_quast_report.tsv ; fi

    if [ -f "quast/${prefix}/transposed_report.tsv" ]
    then
      head -n 1 quast/${prefix}/transposed_report.tsv | awk '{print "sample\\t" \$0 }' > quast/${prefix}/transposed_report.tsv.tmp
      tail -n 1 quast/${prefix}/transposed_report.tsv | awk -v sample=${prefix} '{print sample "\\t" \$0}' >> quast/${prefix}/transposed_report.tsv.tmp
      mv quast/${prefix}/transposed_report.tsv.tmp ${fin}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: \$(quast.py --version 2>&1 | sed 's/^.*QUAST v//; s/ .*\$//')
    END_VERSIONS
  """
}
