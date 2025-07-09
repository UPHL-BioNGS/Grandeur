process FASTANI {
  tag           "${meta.id}"
  label         "process_medium"
  container     'staphb/fastani:1.34'
  
  input:
  tuple val(meta), file(contigs), file(genomes)

  output:
  tuple val(meta), file("fastani/*_fastani.csv"), emit: results, optional: true
  tuple val(meta), env("top_hit"), path("top_hit/*"), emit: top_hit, optional: true
  path "fastani/*_fastani_len.csv", emit: top_len, optional: true
  path "fastani/*", emit: everything
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  def ends   = genomes.collect { it.Name[-6..-1] }.flatten().unique().join(' *')
  """
    mkdir -p fastani logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log
    ls *${ends} | grep -v ${contigs} | sort > reference_list.txt

    fastANI ${args} \
      --threads ${task.cpus} \
      -q ${contigs} \
      --rl reference_list.txt \
      -o fastani/${prefix}.txt \
      | tee -a \$log_file

    echo "sample,query,reference,ANI estimate,total query sequence fragments,fragments aligned as orthologous matches" > fastani/${prefix}_fastani.csv
    cat fastani/${prefix}.txt | sed 's/,//g' | tr "\\t" "," | awk -v sample=${prefix} '{ print sample "," \$0 }' >> fastani/${prefix}_fastani.csv

    top_hit=\$(head -n 2 fastani/${prefix}_fastani.csv | tail -n 1 | cut -f 3 -d , )
    if [ -f "\$top_hit" ]
    then
      mkdir -p top_hit
      cp \$top_hit top_hit/.
      gzip -d top_hit/*.gz || ls top_hit
      chmod 664 top_hit/*

      top_len=\$(grep -v ">" top_hit/* | wc -c)

      echo "sample,top_hit,top_len" > fastani/${prefix}_fastani_len.csv
      echo "${prefix},\$top_hit,\$top_len" >> fastani/${prefix}_fastani_len.csv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      fastani: \$(fastANI --version 2>&1 | sed 's/version//;')
    END_VERSIONS
  """
}
