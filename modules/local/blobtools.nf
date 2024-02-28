process blobtools_create {
  tag           "${meta.id}"
  label         "process_medium"
  publishDir    params.outdir, mode: 'copy'
  container     'chrishah/blobtools:v1.1.1'
  time          '45m'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  
  input:
  tuple val(meta), file(contig), file(blastn), file(bam)

  output:
  tuple val(meta), file("blobtools/*.blobDB.json"), emit: json
  path "blobtools/*.bam.cov"                      , emit: files
  path "logs/${task.process}/*.log"               , emit: log
  path "versions.yml"                             , emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p blobtools logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    blobtools create ${args} \
      -o blobtools/${prefix} \
      -i ${contig} \
      -b ${bam[0]} \
      -t ${blastn} \
      | tee -a \$log_file
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtools: \$(blobtools -v)
    END_VERSIONS
  """
}

process blobtools_view {
  tag           "${meta.id}"
  publishDir    params.outdir, mode: 'copy'
  container     'chrishah/blobtools:v1.1.1'
  time          '10m'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  
  input:
  tuple val(meta), file(json)

  output:
  tuple val(meta), file("blobtools/*.blobDB.table.txt"), emit: file
  path "logs/${task.process}/*.log"                    , emit: log
  path "versions.yml"                                  , emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p blobtools logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    blobtools view ${args} \
      -i ${json} \
      -o blobtools/ \
      | tee -a \$log_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtools: \$(blobtools -v)
    END_VERSIONS
  """
}

process blobtools_plot {
  tag           "${meta.id}"
  publishDir    params.outdir, mode: 'copy'
  container     'chrishah/blobtools:v1.1.1'
  time          '10m'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  
  input:
  tuple val(meta), file(json)

  output:
  path "blobtools/*"                                , emit: files
  tuple val(meta), file("blobtools/*_blobtools.txt"), emit: results
  path "blobtools/*_summary.txt"                    , emit: collect
  path "logs/${task.process}/*.log"                 , emit: log
  path "versions.yml"                               , emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: '--format png -r species'
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p blobtools logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    blobtools plot ${args} \
      -i ${json} \
      -o blobtools/ \
      | tee -a \$log_file

    grep "^# " blobtools/${prefix}*.stats.txt | \
      sed 's/# //g' | \
      awk '{print "sample\t" \$0 }' > blobtools/${prefix}_summary.txt

    grep -v "^#" blobtools/${prefix}*.stats.txt | \
      sed 's/%//g' | \
      tr " " "_" | \
      awk -v sample=${prefix} '{if (\$13 >= 5.0 ) print sample "\\t" \$0}' | \
      tr " " "\\t" | \
      sort -k 14rn,14 >> blobtools/${prefix}_summary.txt

    grep -vw all blobtools/${prefix}_summary.txt > blobtools/${prefix}_blobtools.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtools: \$(blobtools -v)
    END_VERSIONS
  """
}