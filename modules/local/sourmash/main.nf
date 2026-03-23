process SOURMASH {
  tag        "${meta.id}"
  label      "process_medium"
  container  'quay.io/biocontainers/sourmash:4.8.4--pyhdfd78af_0'

  input:
  tuple val(meta), file(reads), file(reference)

  output:
  path "sourmash/*.search.csv",                    emit: search
  tuple val(meta), file("sourmash/*.summary.csv"), emit: results
  path "logs/${task.process}/*.log",               emit: log
  path "versions.yml",                             emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
    def args        = task.ext.search_args ?: ""
    def prefix      = task.ext.prefix      ?: "${meta.id}"
    def input_files = (reads instanceof List) ? reads.join(" ") : reads

    """
    mkdir -p sourmash logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    # Build sourmash signature
    sourmash sketch ${moltype} \\
      -k ${ksize} \\
      --scaled ${scaled} \\
      -p ${task.cpus} \\
      -o ${prefix}.sig \\
      ${input_files} \\
      2> sourmash_${prefix}.err | tee -a \$log_file

    # Search against user-provided database
    sourmash search \\
      ${prefix}.sig \\
      ${reference} \\
      ${args_search} \\
      --csv sourmash/${prefix}.search.csv \\
      -p ${task.cpus} \\
      | tee -a \$log_file

    # Build standardized summary file
    echo "sample,match,name,similarity,containment" > sourmash/${prefix}.summary.csv

    tail -n +2 sourmash/${prefix}.search.csv | \\
      awk -F',' -v sample=${prefix} \\
      '{print sample "," \$1 "," \$2 "," \$3 "," \$4}' \\
      >> sourmash/${prefix}.summary.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      sourmash: \$( sourmash --version )
    END_VERSIONS
  """
}
