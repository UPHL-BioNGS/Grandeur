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
  def args   = task.ext.args   ?: "-v 0 -d 0.25"
  def prefix = task.ext.prefix ?: "${meta.id}"
  if ( reference =~ "input" ) {
    """
    mkdir -p mash logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    mash dist ${args} \
      -p ${task.cpus} \
      /db/RefSeqSketchesDefaults.msh \
      ${msh} | \
      sort -gk3 \
      > mash/${prefix}.mashdist.txt

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
      ${msh} | \
      sort -gk3 \
      > mash/${prefix}.mashdist.txt

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
