process SKANI_DIST {
    tag           "${meta.id}"
    label         "process_medium"
    container     'staphb/skani:0.3.1'

    input:
    tuple val(meta), file(contigs)
    path db_folder

    output:
    tuple val(meta), file(contigs), file("skani/*.tsv"), emit: hits, optional: true
    tuple val(meta), file("skani/*.tsv"), emit: results, optional: true
    tuple val(meta), file("Mycobacteri/*"), emit: myco, optional: true
    tuple val(meta), file("gas/*"), emit: gas, optional: true
    tuple val(meta), file("{Klebsiella,Enterobacter,Serratia}/*"), emit: kleb, optional: true
    tuple val(meta), file("Legionella/*"), emit: legionella, optional: true
    tuple val(meta), file("Streptococcus/*"), emit: strep, optional: true
    tuple val(meta), file("Acinetobacter/*"), emit: acinetobacter, optional: true
    tuple val(meta), file("Salmonella/*"), emit: salmonella, optional: true
    tuple val(meta), file("Escherichia/*"), emit: ecoli, optional: true
    tuple val(meta), file("Vibrio/*"), emit: vibrio, optional: true
    tuple val(meta), file("Neisseriac/*"), emit: gc, optional: true
    path "skani/*txt", emit: skani, optional: true
    path "top_hit/*", emit: top_hit, optional: true
    path "logs/${task.process}/*.log", emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: '--short-header -s 90'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p skani top_hit logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    skani dist ${args} \
        -q ${contigs} \
        -r ${db_folder}/* \
        -t ${task.cpus} \
        -o skani/${prefix}_skani.tsv \
        | tee -a \$log_file

    # ensure prefix is in summary file
    head -n 1 skani/${prefix}_skani.tsv | awk '{print "sample\\t" \$0}' > skani/${prefix}_skani.txt
    tail -n +2 skani/${prefix}_skani.tsv | awk '{print "${prefix}\\t" \$0}' >> skani/${prefix}_skani.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        skani: \$(skani --version | awk '{print \$NF}')
    END_VERSIONS

    # getting the values of all the top hits
    tail -n +2 skani/${prefix}_skani.tsv | sort -rk 3 | head -n 1 | cut -f 1 > top_hit/${prefix}_top_hit.tsv

    # using the skani results to separate contigs for subtyping

    for organism in Salmonella Escherichia Klebsiella Enterobacter Serratia Legionella Vibrio Acinetobacter Mycobacteri Neisseria
    do
        if grep -q "^\$organism" skani/${prefix}_skani.tsv
        then
            echo "Contigs from ${prefix} matched to \${organism} in SKANI results. Copying contigs to \${organism}/ directory." >> \$log_file
            mkdir -p \${organism} 
            cp ${contigs} \${organism}/
        fi
    done

    if grep -q "^Shigella" skani/${prefix}_skani.tsv
    then
        echo "Contigs from ${prefix} matched to Shigella in SKANI results. Copying contigs to Escherichia/ directory." >> \$log_file
        mkdir -p Escherichia
        cp ${contigs} Escherichia/
    fi

    if grep -q "^Streptococcus" skani/${prefix}_skani.tsv
    then
        if grep "pyogenes\\|dysgalactiae\\|anginosus" skani/${prefix}_skani.tsv | grep -q "^Streptococcus"
        then
            echo "Contigs from ${prefix} matched to Streptococcus pyogenes/dysgalactiae/anginosus in SKANI results. Copying contigs to gas/ directory." >> \$log_file
            mkdir -p gas
            cp ${contigs} gas/
        fi
        if grep "pneumoniae" skani/${prefix}_skani.tsv | grep -q "^Streptococcus"
        then
            echo "Contigs from ${prefix} matched to Streptococcus pneumoniae in SKANI results. Copying contigs to Streptococcus/ directory." >> \$log_file
            mkdir -p Streptococcus
            cp ${contigs} Streptococcus/
        fi
    fi
    """
}