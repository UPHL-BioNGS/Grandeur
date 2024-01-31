process fastani {
  tag           "${meta.id}"
  label         "process_medium"
  stageInMode   "copy"
  publishDir    path: params.outdir, mode: 'copy', pattern: 'logs/*/*log'
  publishDir    path: params.outdir, mode: 'copy', pattern: 'fastani/*' 
  container     'staphb/fastani:1.34'
  time          '10m'
  //errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  
  input:
  tuple val(meta), file(contigs), file(genomes)

  output:
  tuple val(meta), file("fastani/*_fastani.csv"), emit: results, optional: true
  tuple val(meta), env(top_hit), path("top_hit/*"), emit: top_hit, optional: true
  path "fastani/*", emit: everything
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  def ref    = genomes.join(",")
  """
    mkdir -p fastani logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    echo ${ref} | tr "," "\\n" | sort > reference_list.txt

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
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      fastani: \$(fastANI --version 2>&1 | sed 's/version//;')
    END_VERSIONS
  """
}
