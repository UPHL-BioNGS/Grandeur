process SPESTIMATOR {
    tag           "${meta.id}"
    label         "process_medium"
    container     'staphb/spestimator:latest'

    input:
    tuple val(meta), file(contigs)

    output:
    path "spestimator/*.tsv",          emit: results
    path "spestimator/*",              emit: files
    path "logs/${task.process}/*.log", emit: log
    path "versions.yml",               emit: versions
    val meta,                          emit: meta

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
      --assembly ${contigs} \
      --output spestimator/${prefix}_spestimator.tsv \
      --threads ${task.cpus} \
      | tee -a \$log_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spestimator: \$(spestimator --version | awk '{print \$2}')
    END_VERSIONS
    """
}