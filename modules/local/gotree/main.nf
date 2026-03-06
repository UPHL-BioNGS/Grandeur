process GOTREE {
  tag           "${analysis}"
  label         "process_medium"
  container     'staphb/gotree:0.5.1'

  
  input:
  tuple val(analysis), file(newick)

  output:
  path "gotree/*.png", emit: for_multiqc
  path "gotree/*", emit: results
  path "gotree/*_stats_all.tsv", emit: stats
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${analysis}"
  """
    mkdir -p gotree logs/${task.process}
    log_file=logs/${task.process}/${analysis}.${task.process}.${workflow.sessionId}.log

    gotree reroot midpoint -i ${newick} | \
      gotree draw png ${args} -o gotree/${prefix}.png | \
      tee -a \$log_file
    gotree stats -i ${newick} -o gotree/${prefix}_stats.tsv

    head -n 1  gotree/${prefix}_stats.tsv | awk '{print "sample\\t"    \$0 }' >  gotree/${prefix}_stats_all.tsv
    tail -n +2 gotree/${prefix}_stats.tsv | awk '{print "${prefix}\\t" \$0 }' >> gotree/${prefix}_stats_all.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      gotree: \$(gotree version | awk '{print \$NF}')
    END_VERSIONS
  """
}