process MASH_SCREEN {
  tag        "${meta.id}"
  label      "process_medium"
  container  'staphb/mash:2.3'

  input:
  tuple val(meta), file(reads), file(reference)

  output:
  tuple val(meta), file("mash/*.mashscreen.txt"), emit: screen, optional: true
  path "mash/*_mashscreen_summary.csv",           emit: screen_results, optional: true
  path "versions.yml",                            emit: versions
  // prints to stdout, so no logs

  when:
  task.ext.when == null || task.ext.when

  script:
    def args        = task.ext.args        ?: "-v 0 -i 0.9 -w" 
    def better_ref  = reference.toString().contains("input") ? "/db/RefSeqSketchesDefaults.msh" : reference
    def org_cmd     = reference.toString().contains("input") ? "echo \$line | awk '{print \$5}' | cut -f 8 -d - | cut -f 1,2 -d _ | cut -f 1 -d ." : "echo \$line | awk '{print \$5}' | cut -f 1,2 -d _ | cut -f 1 -d ."
    def prefix      = task.ext.prefix      ?: "${meta.id}"

    """
    mkdir -p mash

    # mash screen
    mash screen ${args} \
      -p ${task.cpus} \
      ${better_ref} \
      ${reads[0]} \
      > mash/${prefix}.mashscreen.txt

    if [[ -s "mash/${prefix}.mashscreen.txt" ]]
    then 
      echo "sample,identity,shared-hashes,median-multiplicity,p-value,query-ID,organism" > mash/${prefix}_mashscreen_summary.csv

      while read line
      do
        organism2=\$(${org_cmd})
        echo \$line | \
          sed 's/,//g' | \
          awk -v sample=${prefix} \
          -v org=\$organism2 \
          '{print sample "," \$1 "," \$2 "," \$3 "," \$4 "," \$5 "," org}' >> mash/${prefix}_mashscreen_summary.csv
      done < mash/${prefix}.mashscreen.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      mash: \$( mash --version )
    END_VERSIONS
  """
}