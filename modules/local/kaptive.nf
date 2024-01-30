process kaptive {
  tag           "${meta.id}"
  label         "process_medium"
  stageInMode   "copy"
  publishDir    path: params.outdir, mode: 'copy'
  container     'staphb/kaptive:2.0.5'
  time          '30m'
  //errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}

  input:
  tuple val(meta), file(contigs), val(flag)

  output:
  path "kaptive/*_kaptive.csv", emit: collect
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  (task.ext.when == null || task.ext.when) && flag =~ 'found'

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    echo "1"
    mkdir -p kaptive logs/${task.process}
    echo "2"
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log
    
    echo "9"
    for reference in /kaptive/reference_database/*
    do
      echo \$reference    
        kaptive ${args} \
            --k_refs \$reference \
            --assembly ${contigs} \
            | tee -a \$log_file 

    done
    exit 1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kaptive: \$( echo \$(kleborate --version | sed 's/Kleborate v//;'))
    END_VERSIONS
  """
}