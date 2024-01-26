process serotypefinder {
  tag           "${sample}"
  label         "process_medium"
  label         "medcpus"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/serotypefinder:2.0.1'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '10m'
  
  when:
  flag =~ 'found'

  input:
  tuple val(sample), file(file), val(flag), file(script)

  output:
  path "serotypefinder/${sample}/*"                                     , emit: files
  path "serotypefinder/${sample}_serotypefinder.tsv"                    , emit: collect, optional: true
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"       , emit: log
  path  "versions.yml"                          , emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
      def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
  '''
    mkdir -p serotypefinder/!{sample} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    serotypefinder.py !{params.serotypefinder_options} \
      -i !{file} \
      -o serotypefinder/!{sample} \
      -x \
      | tee -a $log_file

    cp serotypefinder/!{sample}/results_tab.tsv serotypefinder/!{sample}_serotypefinder.tsv

    python3 !{script} serotypefinder/!{sample}/results_tab.tsv serotypefinder/!{sample}_serotypefinder.tsv serotypefinder !{sample}
  '''
}
