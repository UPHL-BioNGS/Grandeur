process AMRFINDER {

    tag       "${meta.id}"
    label     "process_high"
    container 'staphb/ncbi-amrfinderplus:4.0.22-2025-03-25.1'

    input:
    tuple val(meta), file(contigs), val(genus), val(species)

    output:
    path "amrfinder/*_amrfinder.txt", emit: collect, optional: true
    val meta, emit: meta
    path "logs/*/*.log", emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: '--plus'
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    set -euo pipefail

    mkdir -p amrfinder logs/${task.process}
    log_file="logs/${task.process}/${prefix}.${workflow.sessionId}.log"

    echo "Running AMRFinder on ${contigs} for ${genus} ${species}" >> "\${log_file}"

    # Fetch list of supported organisms from AMRFinder
    organism_list=\$(amrfinder -l | tr " " "\\n")

    # Try to match genus and species
    organism_match=\$( (echo "\$organism_list" | grep -i "${genus}" | grep -i "${species}" | sed 's/,//g' | head -n 1) || true )

    # Fallback to just genus
    if [ -z "\$organism_match" ]; then
        organism_match=\$( (echo "\$organism_list" | grep -i "${genus}" | sed 's/,//g' | head -n 1) || true )
    fi

    # Special fallback for Shigella
    if [ -z "\$organism_match" ] && [ "${genus}" = "Shigella" ]; then
        organism_match="Escherichia"
        echo "[WARN] No match found for ${genus} ${species}; using Escherichia as fallback" >> "\${log_file}"
    fi
    
    # Build organism flag
    if [ -n "\$organism_match" ]; then
        organism_flag="--organism=\$organism_match"
        echo "[INFO] Using organism match: \$organism_match" >> "\${log_file}"
    else
        organism_flag=""
        echo "[WARN] No organism match found; running AMRFinder without --organism" >> "\${log_file}"
    fi

    # Log and run the command
    echo "[INFO] Running AMRFinder..." >> "\${log_file}"

    amrfinder ${args} \\
        --nucleotide ${contigs} \\
        --threads ${task.cpus} \\
        --name ${prefix} \\
        --output amrfinder/${prefix}_amrfinder.txt \\
        \${organism_flag} | tee -a "\${log_file}"

    echo "[INFO] AMRFinder finished." >> "\${log_file}"

    # Capture version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinderplus: \$(amrfinder --version)
        amrfinderplus-database: \$(amrfinder --version)
    END_VERSIONS
    """
}
