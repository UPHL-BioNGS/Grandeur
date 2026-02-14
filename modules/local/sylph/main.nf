process SYLPH{
    tag           "${meta.id}"
    label         "process_medium"
    container     'staphb/skani:latest'

    input:
    tuple val(meta), file(contigs), file(references)

    output:
    tuple val(meta), file("skani/*.tsv"), emit: results
    path "skani/*",                       emit: files
    path "logs/${task.process}/*.log",    emit: log
    path "versions.yml",                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Join the list of reference files into a space-separated string
    def refs   = references.join(" ") 
    """
    mkdir -p skani logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    skani dist \
      ${args} \
      --threads ${task.cpus} \
      --query ${contigs} \
      --ref ${refs} \
      --output skani/${prefix}_skani.tsv \
      | tee -a \$log_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        skani: \$(skani --version | sed 's/skani //')
    END_VERSIONS
    """
}