process SHIGAPASS {
    tag           "${meta.id}"
    label         "process_medium"
    container     'staphb/shigapass:1.5.0'

    input:
    tuple val(meta), file(contigs)

    output:
    path "shigapass/*_shigapass.tsv", emit: summary, optional: true
    tuple val(meta), file("shigapass/*/*"), emit: all_files, optional: true
    path "logs/${task.process}/*.log", emit: log
    path "versions.yml", emit: versions

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

    # ensure prefix is in summary file
    head -n 1  shigapass/${prefix}/ShigaPass_summary.csv | sed 's/;/\\t/g' | awk '{print "sample\\t"    \$0 }' >  shigapass/${prefix}_shigapass.tsv
    tail -n +2 shigapass/${prefix}/ShigaPass_summary.csv | sed 's/;/\\t/g' | awk '{print "${prefix}\\t" \$0 }' >> shigapass/${prefix}_shigapass.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shigapass: \$(ShigaPass.sh -v | awk '{print \$NF}')
    END_VERSIONS
    """
}