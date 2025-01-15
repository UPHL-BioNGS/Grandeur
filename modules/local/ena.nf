process ENA_DOWNLOAD {
    tag           "${SRR}"
    label         "process_single"
    container     'staphb/viridian:1.3.0'
    
    input:
    val(SRR)

    output:
    tuple val(SRR), file("reads/${SRR}_{1,2}.fastq.gz"), emit: fastq
    path "logs/${task.process}/*.log", emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p reads logs/${task.process}
    log_file=logs/${task.process}/${SRR}.${workflow.sessionId}.log
        
    enaDataGet -f fastq ${SRR}

    mv */*fastq.gz reads/.

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        enaDataGet: \$( enaDataGet -v | awk '{print \$NF}' )
    END_VERSIONS
    """
}