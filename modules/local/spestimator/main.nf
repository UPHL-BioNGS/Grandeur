process SPESTIMATOR {
    tag           "${meta.id}"
    label         "process_medium"
    container     'staphb/spestimator:0.3.0.233'

    input:
    tuple val(meta), file(contigs)

    output:
    tuple val(meta), file("spestimator/*.csv"), emit: results, optional: true
    path "logs/${task.process}/*.log", emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p spestimator logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    spestimator \
        ${args} \
        --input ${contigs} \
        --output spestimator/${prefix}_spestimator.csv \
        --threads ${task.cpus} \
        | tee -a \$log_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spestimator: \$(spestimator --version | awk '{print \$2}')
    END_VERSIONS
    """
}