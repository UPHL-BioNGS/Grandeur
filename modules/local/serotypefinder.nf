process serotypefinder {
  tag           "${meta.id}"
  label         "process_medium"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/serotypefinder:2.0.1'
  maxForks      10
  time          '10m'
  //errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  
  when:
  flag =~ 'found'

  input:
  tuple val(meta), file(file), val(flag), file(script)

  output:
  path "serotypefinder/*/*"                                     , emit: files
  path "serotypefinder/*_serotypefinder.tsv"                    , emit: collect, optional: true
  path "logs/${task.process}/*.log"       , emit: log
  path "versions.yml"                          , emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
      def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p serotypefinder/${prefix} logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : ${task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    serotypefinder.py ${params.serotypefinder_options} \
      -i ${file} \
      -o serotypefinder/${prefix} \
      -x \
      | tee -a $log_file

    cp serotypefinder/${prefix}/results_tab.tsv serotypefinder/${prefix}_serotypefinder.tsv

    python3 ${script} serotypefinder/${prefix}/results_tab.tsv serotypefinder/${prefix}_serotypefinder.tsv serotypefinder ${prefix}
  """
}
