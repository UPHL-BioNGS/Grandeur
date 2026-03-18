process KSNP4 {
    tag           "kSNP4 Tree"
    label         "process_high"
    container     'staphb/ksnp4:4.1'

    input:
    file(fastas)

    output:
    path "ksnp4/*",                        emit: files
    path "ksnp4/tree.parsimony.tre",       emit: newick, optional: true
    path "ksnp4/tree.ML.tre",              emit: tree_ml, optional: true
    path "ksnp4/SNPs_all_matrix.fasta",    emit: snp_matrix, optional: true
    path "logs/${task.process}/*.log",     emit: log
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: '-k 19 -core -ML -NJ'
    def prefix = task.ext.prefix ?: "ksnp4_run"
    """
    mkdir -p ksnp4 logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    # kSNP4 requires an input file listing the paths to the genomes    
    for f in *.fa *.fasta *.fna; do
        if [ -f "\$f" ]; then
            # Clean filename to create ID
            ID=\$(basename "\$f" | sed 's/.fasta//;s/.fa//;s/.fna//')
            echo -e "\$PWD/\$f\t\$ID" >> ksnp4_input_list.txt
        fi
    done
    
    kSNP4 \
        ${args} \
        -in ksnp4_input_list.txt \
        -outdir ksnp4 \
        -CPU ${task.cpus} \
        | tee -a \$log_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kSNP4: \$(kSNP4 -version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """
}