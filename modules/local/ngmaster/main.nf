process NGMASTER {
    tag           "${meta.id}"
    label         "process_medium"
    container     'staphb/ngmaster:1.1.1'

    input:
    tuple val(meta), file(contigs)

    output:
    tuple val(meta), file("ngmaster/*"), emit: files
    path "*ngmaster.csv", emit: collect
    path "logs/${task.process}/*.log", emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: '--csv'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ngmaster logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    cat ${contigs} > input_${prefix}.fasta

    ngmaster \
        ${args} \
        input_${prefix}.fasta \
        > ngmaster/${prefix}_ngmaster.csv \
        2>> \$log_file

        # ensure prefix is in summary file
    head -n 1  ngmaster/${prefix}_ngmaster.csv | awk '{print "sample,"    \$0 }' >  ${prefix}_ngmaster.csv
    tail -n +2 ngmaster/${prefix}_ngmaster.csv | awk '{print "${prefix}," \$0 }' >> ${prefix}_ngmaster.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ngmaster: \$(echo \$(ngmaster --version 2>&1 | awk '{print \$2}'))
    END_VERSIONS
    """
}