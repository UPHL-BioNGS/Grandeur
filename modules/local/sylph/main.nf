process SYLPH {
    tag           "${meta.id}"
    label         "process_medium"
    container     'staphb/sylph:0.9.0'

    input:
    tuple val(meta), file(reads), path(db)

    output:
    tuple val(meta), file("sylph/*sylph.tsv"), emit: tsv, optional: true
    path "sylph/*sylph_results.tsv", emit: results, optional: true
    path "download/*.txt", emit: for_download, optional: true
    path "logs/${task.process}/*.log", emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args        ?: "" 
    def is_fastq     = (reads instanceof List) || reads.name.toString().matches('.*\\.(fastq|fq)(\\.gz)?$')
    def sylph_input  = is_fastq ? "-1 ${reads[0]} -2 ${reads[1]}" : "-r ${reads}"   
    def prefix       = task.ext.prefix      ?: "${meta.id}"
    """
    mkdir -p sylph download logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    sylph profile \
        ${db} \
        ${sylph_input} \
        -t ${task.cpus} \
        ${args} \
        -o sylph/${prefix}_sylph.tsv \
        | tee -a \$log_file

    head -n 1  sylph/${prefix}_sylph.tsv | awk '{print "sample\\t"    \$0 }' >   sylph/${prefix}_sylph_results.tsv
    tail -n +2 sylph/${prefix}_sylph.tsv | awk '{print "${prefix}\\t" \$0 }' >>  sylph/${prefix}_sylph_results.tsv

    # extract the genome accessions to download if set
    cut -f 2 sylph/${prefix}_sylph.tsv \
        | tail -n +2  \
        | grep "G" \
        | rev \
        | cut -f 1 -d / \
        | rev \
        | sed "s/_genomic.fna.gz//g" \
        > download/download_${prefix}.txt

    [ -s download/download_${prefix}.txt ] || rm download/download_${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph: \$(sylph --version | awk '{print \$NF}')
    END_VERSIONS
    """
}