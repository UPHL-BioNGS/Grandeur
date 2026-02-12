process SHIGAPASS {
    tag           "${meta.id}"
    label         "process_medium"
    container     'staphb/shigapass:latest'

    input:
    tuple val(meta), file(contigs)

    output:
    path "shigapass/*_summary.csv", emit: summary
    path "shigapass/*", emit: files
    path "logs/${task.process}/*.log", emit: log
    path "versions.yml", emit: versions
    val meta, emit: meta

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p shigapass logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    shigapass \
      --contigs ${contigs} \
      --out_dir shigapass \
      --threads ${task.cpus} \
      | tee -a \$log_file

    # Rename output to include sample name if not already present or to standardize
    mv shigapass/summary.csv shigapass/${prefix}_shigapass_summary.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shigapass: \$(shigapass --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """
}