process MASH_SKETCH {
  tag           "${meta.id}"
  label         "process_medium"
  container     'staphb/mash:2.3'

  input:
  tuple val(meta), file(files)

  output:
  tuple val(meta), file("mash/*.msh"), emit: msh
  tuple val(meta), file("mash/*.err"), optional: true, emit: err
  path "mash/*estimates.csv", optional: true, emit: summary
  path "logs/${task.process}/*.log", emit: log
  path  "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: ""
  def mode   = files[1] ? "-m 2 $args" : "$args"
  def prefix = task.ext.prefix ?: "${meta.id}"
  def input  = files.join(" ")
  """
  mkdir -p mash logs/${task.process}
  log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log
  err_file=mash/${prefix}.${workflow.sessionId}.err

  cat ${input} | \
    mash sketch ${mode} \
    -o mash/${prefix} - \
    2>> \$err_file | tee -a \$log_file

  genome_size=\$(grep "Estimated genome size:" \$err_file | awk '{print \$NF}')
  coverage=\$(grep "Estimated coverage:"  \$err_file | awk '{print \$NF}')

  echo "sample,mash_estimated_genome_size,mash_estimated_coverage" > mash/${prefix}_mash_estimates.csv
  echo "${prefix},\$genome_size,\$coverage" >> mash/${prefix}_mash_estimates.csv

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    mash: \$( mash --version )
  END_VERSIONS
  """
}

process MASH_DIST {
  tag           "${meta.id}"
  label         "process_medium"
  container     'staphb/mash:2.3'

  input:
  tuple val(meta), file(msh), file(reference)

  output:
  path "mash/*.mashdist.txt",                       emit: mashdist
  tuple val(meta), file("mash/*.summary.mash.csv"), emit: results
  path "logs/${task.process}/*.log",                emit: log
  path  "versions.yml",                             emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: "-v 0 -d 0.5"
  def prefix = task.ext.prefix ?: "${meta.id}"
  if ( reference =~ "input" ) {
    """
    mkdir -p mash logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    mash dist ${args} \
      -p ${task.cpus} \
      /db/RefSeqSketchesDefaults.msh \
      ${msh} \
      > mash/${prefix}.mashdist.txt.tmp

    if [ ! -s "mash/${prefix}.mashdist.txt.tmp" ]
    then
      echo "No mash dist results with args ${args}. Trying again without them."  | tee -a \$log_file
      mash dist \
        -p ${task.cpus} \
        /db/RefSeqSketchesDefaults.msh \
        ${msh} \
        > mash/${prefix}.mashdist.txt.tmp
    fi

    sort -gk3 mash/${prefix}.mashdist.txt.tmp | head -n ${params.mash_max_hits} > mash/${prefix}.mashdist.txt

    echo "sample,reference,query,mash-distance,P-value,matching-hashes,organism" > mash/${prefix}.summary.mash.csv

    while read line
    do
      organism=\$(echo \$line | cut -f 8 -d "-" | cut -f 1,2 -d "_" | cut -f 1 -d ".")
      echo \$line | sed 's/,//g' | awk -v sample=${prefix} -v org=\$organism '{print sample "," \$1 "," \$2 "," \$3 "," \$4 "," \$5 "," org}' >> mash/${prefix}.summary.mash.csv
    done < mash/${prefix}.mashdist.txt

    wc -l mash/${prefix}.mashdist.txt >> \$log_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      mash: \$( mash --version )
    END_VERSIONS
    """
  } else {
    """
    mkdir -p mash logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    mash dist ${args} \
      -p ${task.cpus} \
      ${reference} \
      ${msh} \
      > mash/${prefix}.mashdist.txt.tmp

    if [ ! -s "mash/${prefix}.mashdist.txt.tmp" ]
    then
      echo "No mash dist results with args ${args}. Trying again without them." | tee -a \$log_file
      mash dist \
        -p ${task.cpus} \
        ${reference} \
        ${msh} \
        > mash/${prefix}.mashdist.txt.tmp
    fi

    sort -gk3 mash/${prefix}.mashdist.txt.tmp | head -n ${params.mash_max_hits} > mash/${prefix}.mashdist.txt

    echo "sample,reference,query,mash-distance,P-value,matching-hashes,organism" > mash/${prefix}.summary.mash.csv

    while read line
    do
      organism=\$(echo \$line | cut -f 3,4 -d "_" | cut -f 1 -d ".")
      echo \$line | sed 's/,//g' | awk -v sample=${prefix} -v org=\$organism '{print sample "," \$1 "," \$2 "," \$3 "," \$4 "," \$5 "," org}' >> mash/${prefix}.summary.mash.csv
    done < mash/${prefix}.mashdist.txt

    wc -l mash/${prefix}.mashdist.txt >> \$log_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      mash: \$( mash --version )
    END_VERSIONS
    """
  }
}
