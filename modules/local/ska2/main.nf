process SKA2 {
    tag           "SKA alignment"
    label         "process_medium"
    container     'staphb/ska2:0.5.1'

    input:
    file(contigs)

    output:
    path "ska/*_alignment.aln",        emit: aln, optional: true
    path "ska/*",                      emit: files, optional: true
    path "logs/${task.process}/*.log", emit: log
    path "versions.yml",               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: '--filter no-filter'
    def build_args = task.ext.build_args ?: '-k 31'
    def prefix = task.ext.prefix ?: "ska"
    """
    mkdir -p ska logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    ska build ${build_args} \
        -f ${contigs} \
        -o ska/${prefix}_index
    
    ska align ${args} \
        -o ska/${prefix}_alignment.aln \
        ska/${prefix}_index.skf

    ska --version

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        skani: \$(skani --version | sed 's/skani //')
    END_VERSIONS
    """
}