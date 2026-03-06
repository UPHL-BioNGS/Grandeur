process NGMASTER {
    tag           "${meta.id}"
    label         "process_medium"
    container     'staphb/ngmaster:0.5.8-2023-05'

    input:
    tuple val(meta), file(contigs)

    output:
    path "ngmaster/*.tsv", emit: collect, optional: true
    path "logs/${task.process}/*.log", emit: log
    path "versions.yml", emit: versions
    val meta, emit: meta

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: '--csv'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ngmaster logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    ngmaster \
      ${args} \
      ${contigs} \
      > ngmaster/${prefix}_ngmaster.tsv \
      2>> \$log_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ngmaster: \$(echo \$(ngmaster --version 2>&1 | awk '{print \$2}'))
    END_VERSIONS

    # add sample to column name
    exit 1
    """
}