process SKANI_SKETCH {
    tag           "Sketching reference genomes with SKANI"
    label         "process_medium"
    container     'staphb/skani:0.3.1'
    stageInMode   'copy'

    input:
    path genomes

    output:
    path "skani_db" , emit: db, optional: true
    path "logs/${task.process}/*.log", emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "skani_db"
    """
    mkdir -p decompressed_genomes logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    skani sketch ${args} \
        --separate-sketches \
        ${genomes}/* \
        -o ${prefix} \
        -t ${task.cpus} \
        | tee -a \$log_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        skani: \$(skani --version | awk '{print \$NF}')
    END_VERSIONS
    """
}