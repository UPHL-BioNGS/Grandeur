process AMRFINDER {
  tag           "${meta.id}"
  label         "process_high"
  container     'staphb/ncbi-amrfinderplus:4.0.23-2025-07-16.1'

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
    mkdir -p amrfinder logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    organism=\$(amrfinder -l | tr " " "\\n" | grep -i ${genus} | grep -i ${species} | sed 's/,//g' | head -n 1 )
    if [ -z "\$organism" ] ; then organism=\$(amrfinder -l | tr " " "\\n" | grep -i ${genus} | sed 's/,//g' | head -n 1 ) ; fi
    if [ -n "\$organism" ]
    then
      organism_check="--organism \$organism"
      echo "Top organism result of ${genus} ${species} matched with \$organism" >> \$log_file
    elif [ "${genus}" == "Shigella" ]
    then
      organism_check="--organism Escherichia"
      echo "--organism Escherichia with be used because of top organism result of ${genus}" >> \$log_file
    else
      organism_check=''
      echo "Top organism result of ${genus} ${species} did not match any of the organisms" >> \$log_file
    fi

    amrfinder ${args} \
      --nucleotide ${contigs} \
      --threads ${task.cpus} \
      --name ${prefix} \
      --output amrfinder/${prefix}_amrfinder.txt \
      \$organism_check \
      | tee -a \$log_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinderplus: \$(amrfinder --version)
        amrfinderplus-database: \$(echo \$(echo \$(amrfinder --database amrfinderdb --database_version 2> stdout) | rev | cut -f 1 -d ' ' | rev))
    END_VERSIONS
  """
}
