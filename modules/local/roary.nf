process roary {
    tag           "Core Genome Alignment"
    label         "process_high"
    publishDir    params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    container     'staphb/roary:3.13.0'
    errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
    time          '10h'
    
    input:
    file(contigs)

    output:
    path "roary/*"                                                                       , emit: files
    path "roary/fixed_input_files/*"                                                     , emit: roary_input_files
    tuple path("roary/core_gene_alignment.aln"), path("roary/gene_presence_absence.Rtab"), emit: core_gene_alignment, optional: true
    path "logs/${task.process}/${task.process}.${workflow.sessionId}.log"                , emit: log_files
    path "versions.yml"                                                                  , emit: versions

    shell:
    def args       = task.ext.args   ?: ''
    def prefix     = task.ext.prefix ?: 'roary'
    """
        mkdir -p logs/${task.process}
        log_file=logs/${task.process}/${task.process}.${workflow.sessionId}.log

        roary ${args} \
        -p ${task.cpus} \
        -f roary \
        -e -n \
        *.gff \
        | tee -a \$log_file

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            roary: \$( roary --version )
        END_VERSIONS
    """
}