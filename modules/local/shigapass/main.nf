process SHIGAPASS {
    tag           "${meta.id}"
    label         "process_medium"
    container     'staphb/shigapass:latest'

    input:
    tuple val(meta), file(contigs)

    output:
    path "shigapass/*_summary.csv", emit: summary, optional: true
    path "shigapass/*", emit: all_files, optional: true
    path "logs/${task.process}/*.log", emit: log
    path "versions.yml", emit: versions
    val meta, emit: meta

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p shigapass logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    echo "${contigs}" > list.txt

    ShigaPass.sh \
        ${args} \
        -l list.txt \
        -o shigapass/${prefix} \
        -p \${DB_PATH} \
        -t ${task.cpus} \
        | tee -a \$log_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shigapass: \$(ShigaPass.sh -v | awk '{print \$NF}')
    END_VERSIONS

    # add sample to column name
    exit 1
    """
}