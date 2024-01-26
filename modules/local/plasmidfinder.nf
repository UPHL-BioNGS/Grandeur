process plasmidfinder {
  tag           "$meta.id"
  label         "process_medium"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/plasmidfinder:2.1.6'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '10m'

  input:
  tuple val(sample), file(file), file(script)

  output:
  path "plasmidfinder/${sample}/*"                                    , emit: files
  path "plasmidfinder/${sample}_plasmidfinder.tsv"                    , emit: collect, optional: true
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"     , emit: log
  path  "versions.yml"                          , emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
      def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
  '''
    mkdir -p plasmidfinder/!{sample} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    echo "plasmidfinder.py: no version" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    plasmidfinder.py !{params.plasmidfinder_options} \
      -i !{file} \
      -o plasmidfinder/!{sample} \
      --extented_output \
      | tee -a $log_file

    python3 !{script} plasmidfinder/!{sample}/results_tab.tsv plasmidfinder/!{sample}_plasmidfinder.tsv plasmidfinder !{sample}
  '''
}
