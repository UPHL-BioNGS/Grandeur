process ENA_DOWNLOAD {
    tag           "${SRR}"
    label         "process_single"
    container     'staphb/enabrowsertools:1.7.1'
    
    input:
    val(SRR)

    output:
    tuple val(SRR), file("*/*{1,2}.fastq.gz"), emit: fastq, optional: true
    path "logs/*/*.log", emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${SRR}"
    """
    mkdir -p reads logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log
        
    enaDataGet \
        ${args} \
        -f fastq \
        ${SRR} \
        | tee -a \$log_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        enaDataGet: \$( enaDataGet -v | awk '{print \$NF}' )
    END_VERSIONS
    """
}