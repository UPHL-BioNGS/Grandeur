process CHECKM2 {
    tag           "${meta.id}"
    label         "process_high"
    container     'staphb/checkm2:latest'

    input:
    tuple val(meta), file(contigs), path(db)

    output:
    tuple val(meta), file ("checkm2/*/quality_report.tsv"), emit: results, optional: true
    path "checkm2/*_quality_report.tsv", emit: report, optional: true
    path "checkm2/*/*",                  emit: files
    path "logs/${task.process}/*.log",   emit: log
    path "versions.yml",                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p checkm2 logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    checkm2 predict \
        ${args} \
        --threads ${task.cpus} \
        --input ${contigs} \
        --output-directory checkm2/${prefix} \
        --database_path ${db} \
        | tee -a \$log_file

    # for summary file
    head -n 1  checkm2/${prefix}/quality_report.tsv | awk '{print "sample\\t"    \$0 }' >  checkm2/${prefix}_quality_report.tsv
    tail -n +2 checkm2/${prefix}/quality_report.tsv | awk '{print "${prefix}\\t" \$0 }' >> checkm2/${prefix}_quality_report.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm2: \$(checkm2 --version)
    END_VERSIONS
    """
}